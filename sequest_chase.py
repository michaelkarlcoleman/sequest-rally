#!/usr/bin/env python

'''This program invokes individual SEQUEST (unify_sequest) searches.  It
listens for a connection from sequest-rally, which sends all of the
information necessary to do the search, and returns the search results on the
same connection.  It will continue to process further search requests on the
connection, timing out if none are received within the specified interval.
Then new connections will be accepted, one at a time.

On startup, the program will attempt to listen on a sequence of ports,
depending on the values of --port and --port-count.  If these are 12345 and 4,
for example, it will try ports 12345, 12346, 12347, and 12348, using the first
available.  Since only one process is allowed to listen on each port, this
functions as a form of locking that will prevent more than --port-count
occurrences of this program from running on the host simultaneously.

This program (or unify_sequest) will access ./sequest.params, the sequence
database it mentions, and the ms2 files via the filesystem.
'''

from __future__ import with_statement

__copyright__ = '''
    sequest-rally, a cluster driver for SEQUEST, derived from greylag
    Copyright (C) 2006-2009  Stowers Institute for Medical Research
'''


import cPickle as pickle
import errno
import fcntl
import itertools
import logging; from logging import debug, info, warning
from operator import itemgetter
import optparse
import os
from pprint import pprint, pformat
import re
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import time

from sequest_rally_lib import *


# gc possibly harms performance here, so disable it.  gc only matters for
# cycles, which we (hope we) don't create.  See the gc module docs.
import gc; gc.disable()


def error(s, *args):
    "fatal error --> exit with error"
    # if we're unit testing, just throw an exception
    if __name__ != "__main__":
        raise Exception((s + " (fatal error)") % args)
    logging.error(s, *args)
    sys.exit(1)

class ChaseException(Exception): pass

def chase_error(s, *args):
    "error --> disconnect from client"
    logging.error(s, *args)
    raise ChaseException((s + " (disconnecting)") % args)


def check_file_readable(fn):
    try:
        with open(fn) as f:
            pass
    except OSError, e:
        chase_error("could not read '%s' (%s)" % (fn, e))


def set_parameters(state, arg):
    pwd = arg['pwd']
    assert pwd.startswith('/'), "pwd must be absolute"
    try:
        os.chdir(pwd)
    except OSError, e:
        chase_error("could not cd to '%s' (%s)" % (pwd, e))

    assert arg['command'] in ('unify_sequest', 'unify_sequest-1')
    state['original_directory'] = pwd

    absolute_command = subprocess.Popen("type -P %s" % arg['command'],
                                        stdout=subprocess.PIPE,
                                        shell=True).communicate()[0].strip()
    debug("sequest executable = %s", absolute_command)
    state['command'] = absolute_command

    # Set up a local work directory so that we can use custom versions of
    # sequest.params to search partial sequence databases
    assert 'work_directory' not in state, "directory leak"
    tmp_prefix = os.path.basename(sys.argv[0]) + '-'
    state['work_directory'] = tempfile.mkdtemp(prefix=tmp_prefix)

    # double-check that we're looking at the same sequence database as the
    # master
    db = arg['sequence_db']
    check_file_readable(db)
    # NB: SEQUEST sigsegv's on long (>80?) database filenames
    state['db_fn'] = os.path.join(state['work_directory'], 'db.fasta')
    shutil.copyfile(db, state['db_fn'])

    checksum = arg['sequence_db_checksum']
    checksum1 = file_sha1(state['db_fn'])
    if checksum != checksum1:
        chase_error("database checksum does not match [%s]",
                    (checksum1, checksum))

    # NB: must open as binary because we're calculating seek/tell offsets
    dbfn_contents = file_contents(state['db_fn'], binary=True)
    # list of offsets for each clump in sequence db
    # sentinel at end so we can handle half-open ranges correctly
    state['sequence_db_offsets'] = (arg['sequence_db_clump_index']
                                    + [ len(dbfn_contents) ])
    del dbfn_contents

    CONFIG = 'sequest.params'
    check_file_readable(CONFIG)
    sequest_params = file_contents(CONFIG)
    state['db_partial_fn'] = state['db_fn'] + '.partial'
    wsp, n = re.subn(r'(?m)^\s*database_name\s*=.*$',
                     'database_name = %s' % state['db_partial_fn'],
                     sequest_params)
    if n < 1:
        chase_error("did not find proper 'database_name' line in '%s'", CONFIG)
    with open(os.path.join(state['work_directory'], CONFIG), 'w') as wspf:
        wspf.write(wsp)

    state['prior_search_specs'] = None  # i.e., entire sequence database
    shutil.copyfile(state['db_fn'], state['db_partial_fn'])


# FIX: won't deal well with huge command output (not supposed to happen)

last_search_time = 0

def single_search(options, state, piece):
    """Search piece (part of all of one spectrum)"""
    debug("piece: %s", piece)

    sp_spec = piece[0]

    # each db_spec is None (meaning "entire db"), or (start, end), where start
    # and end are indices into state['sequence_db_offsets'].  (end is
    # exclusive)
    ### each db_spec is None (meaning "entire db"), or (part_count, part_no)
    db_spec = piece[1]

    # create new partial sequence database, if necessary
    if db_spec != state['prior_search_specs']:
        with open(state['db_partial_fn'], 'w') as pdbf:
            if db_spec is None:
                s = file_contents(state['db_fn'])
            else:
                offset_s = state['sequence_db_offsets'][db_spec[0]]
                offset_e = state['sequence_db_offsets'][db_spec[1]]
                s = file_contents(state['db_fn'], offset_s, offset_e-offset_s)
            assert s[0] == '>', "invalid db subpart"
            assert s[-1] == '\n', "invalid db subpart (end)"
            pdbf.write(s)
            debug("db subpart length is %s", len(s))
            state['prior_search_specs'] = db_spec

    spfn = sp_spec[2].split(':', 1)[0]

    # symlink ms2 file into work directory, if necessary
    work_spfn = os.path.join(state['work_directory'], spfn)
    if not os.path.exists(work_spfn):
        check_file_readable(spfn)
        os.symlink(os.path.join(state['original_directory'], spfn), work_spfn)

    # Try pretty hard to kill child on exception/signals.  (SIGIO usually
    # drops the child, but not always.)
    running = False
    try:
        global last_search_time
        limit_time = last_search_time + 1.0/options.rate_limit
        start_time = time.time()
        if start_time < limit_time:
            debug("sleeping %.3fs, per rate limit", limit_time - start_time)
            time.sleep(limit_time - start_time)
            start_time = time.time()
        last_search_time = start_time
        sub = subprocess.Popen([state['command'], sp_spec[2]],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               cwd=state['work_directory'])
        running = True
        out_s, err_s = sub.communicate()
        running = False
        elapsed_time = time.time() - start_time
    finally:
        if running:
            try:
                debug("killing (SIGKILL) child (pid = %s)", sub.pid)
                os.kill(sub.pid, signal.SIGKILL)
            except OSError, e:
                debug("failed to kill child (already dead?) (%s)", e)

    # tell the OS it can discard these pages, because the aggregate size of
    # the spectrum files may be huge, and we want to avoid pushing more
    # important stuff out of memory
    with open(work_spfn) as f:
        posix_fadvise(f.fileno(), 0, 0, POSIX_FADV_DONTNEED)

    if len(err_s) > 0:
        chase_error("got error '%s' from '%s %s'"
                    % (err_s.strip(), state['command'], sp_spec[2]))
    if sub.returncode != 0:
        chase_error("got error return code '%s' (%s) from '%s %s'"
                    % (sub.returncode, out_s.strip(), state['command'],
                       sp_spec[2]))

    # use our search time, because SEQUEST's is buggy
    TIME_FIELD = 5
    if out_s.startswith('\nS\t'):
        parts = out_s.split('\t', TIME_FIELD)
        if len(parts) == TIME_FIELD+1:
            parts[TIME_FIELD-1] = '%.*f' % (ELAPSED_TIME_PRECISION,
                                            elapsed_time)
            out_s = '\t'.join(parts)

    return (sp_spec, db_spec, out_s, elapsed_time)


def perform_search(options, state, arg):
    """Search the pieces (arg), batching searches for the same spectra into
    calls to single_search, for efficiency.
    """
    assert len(arg) > 0

    info("searching %s pieces", len(arg))

    results = [ single_search(options, state, piece) for piece in arg ]

    info("returning results")
    debug("results: %s", results)
    return results


# inspired by Beazley's generator talk
def listen_for_connections(port, port_count):
    """Listen on a port from the range [port, port+port_count) and yield a
    series of connections."""
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)

    for p in range(port, port+port_count):
        try:
            s.bind(('', p))             # socket.INADDR_ANY?
            s.listen(20)                # works as a poor man's queue
            info("listening on port %s", p)
            break
        except socket.error, e:
            if e[0] == errno.EADDRINUSE:
                continue
            raise
    else:
        error("could not listen, all specified ports in use")

    while True:
        client_socket, client_addr = s.accept()
        info("received connection from %s", client_addr)

        # try to better notice reboots, net failures, etc
        try:
            client_socket.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE,
                client_socket.getsockopt(socket.SOL_SOCKET,
                                         socket.SO_KEEPALIVE) | 1)
        except socket.error:
            pass

        # get ready for SIGIO handling on this socket
        fcntl.fcntl(client_socket.fileno(), fcntl.F_SETOWN, -os.getpgrp())

        yield client_socket


# FIX: give cleaner diagnostics on protocol failure?


def reset_state(state):
    state.clear()


def set_socket_abort(s, enabled):
    """Enable or disable SIGIO for this socket."""
    O_ASYNC = 020000
    flags = fcntl.fcntl(s.fileno(), fcntl.F_GETFL)
    if enabled:
        flags |= O_ASYNC
    else:
        flags &= ~O_ASYNC
    fcntl.fcntl(s.fileno(), fcntl.F_SETFL, flags)


def handle_command(client_socket, options, state, command, arg):
    """Handle command/arg from client, returning (response, response_arg) or
    None.  State that persists across commands is stored in 'state'.
    """
    if command == 'parameters':
        reset_state(state)
        set_parameters(state, arg)
        return None
    elif command == 'search':
        try:
            # With abort enabled, we'll receive a SIGIO on socket activity,
            # as will our possible child unify_sequest (which will die
            # immediately, as desired).
            set_socket_abort(client_socket, True)

            # FIX: revamp to return individual spectrum results as they are
            # known?
            results = perform_search(options, state, arg)
        finally:
            set_socket_abort(client_socket, False)
        return 'found', results
    else:
        assert False, "unknown command '%s'" % command


def signal_handler(signum, frame):
    debug("received signal %s", signum)
    raise socket.error("received signal %s" % signum)


def reap_child_processes():
    """Wait for any children, to avoid leaving zombie processes."""
    while True:
        pid = 0
        try:
            pid, status = os.waitpid(-1, os.WNOHANG)
        except OSError, e:
            if e.errno != errno.ECHILD:
                raise
        if pid == 0:
            break
        debug("child %s exited with status %s", pid, status)


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options]",
                                   description=__doc__, version=VERSION)
    pa = parser.add_option
    DEFAULT_PORT = 20078
    pa("--port", dest="port", type="int", default=DEFAULT_PORT,
       help="first listener port [default=%s]" % DEFAULT_PORT)
    DEFAULT_PORT_COUNT = 4
    pa("--port-count", dest="port_count", type="int",
       default=DEFAULT_PORT_COUNT,
       help="number of ports to try [default=%s]" % DEFAULT_PORT_COUNT,
       metavar="COUNT")
    DEFAULT_TIMEOUT = 300
    pa("--timeout", dest="timeout", type="int",
       default=DEFAULT_TIMEOUT,
       help=("inactivity timeout, in whole seconds [default=%ss]"
             % DEFAULT_TIMEOUT), metavar="SECONDS")
    DEFAULT_SESSION_LIMIT = 3600
    pa("--session-limit", dest="session_limit", type="float",
       default=DEFAULT_SESSION_LIMIT,
       help="maximum session length, in seconds [default=%ss]; after this,"
       " chaser will disconnect, providing simple sharing between multiple"
       " sequest-rally sessions" % DEFAULT_SESSION_LIMIT, metavar="SECONDS")
    DEFAULT_RATE_LIMIT = 12.0
    pa("--rate-limit", dest="rate_limit", type="float",
       default=DEFAULT_RATE_LIMIT,
       help="maximum spectrum search rate, in spectra per second"
       " [default=%ss]; used to avoid swamping I/O channels"
       % DEFAULT_RATE_LIMIT, metavar="SPECTRA-PER-SECOND")
    pa("-q", "--quiet", action="store_true", dest="quiet", help="no warnings")
    pa("-p", "--show-progress", action="store_true", dest="show_progress",
       help="show running progress")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("-l", "--logfile", dest="logfile",
       help="log to FILE instead of stderr", metavar="FILE")
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    pa("--debug", action="store_true", dest="debug",
       help="output debugging info")
    pa("--profile", action="store_true", dest="profile",
       help="dump Python profiling output to './sequest-chase.prof.<pid>'")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)

    if (len(args) > 0 or options.timeout <= 0 or options.session_limit <= 0
        or options.rate_limit <= 0 or options.port_count < 1
        or options.port <= 0 or options.port + options.port_count >= 65536):
        parser.print_help()
        sys.exit(1)

    set_logging(options)

    signal.signal(signal.SIGIO, signal_handler)
    signal.signal(signal.SIGALRM, signal_handler)

    info("starting on %s", socket.gethostname())

    state = {}
    for client in listen_for_connections(options.port, options.port_count):
        try:
            client_f_to = client.makefile('w')
            client_f_from = client.makefile('r')

            print >> client_f_to, "sequest-chase %s ready" % VERSION
            client_f_to.flush()

            reset_state(state)
            session_limit = time.time() + options.session_limit

            while True:
                if time.time() > session_limit:
                    info("session limit reached, disconnecting")
                    break
                signal.alarm(options.timeout)
                command_line = client_f_from.readline()

                if not command_line:
                    info("client closed connection")
                    break
                # Allow chasers to be remotely killed.
                # FIX: need a more secure method.
                if command_line.startswith('DIE'):
                    error("DIE request received")
                command, arg_length = command_line.split()
                arg = pickle.loads(client_f_from.read(int(arg_length)))
                signal.alarm(0)

                try:
                    r = handle_command(client, options, state, command, arg)
                except ChaseException, e:
                    r = ('error', e)
                if not r:
                    continue
                response, response_arg = r

                signal.alarm(options.timeout)
                if response == 'error':
                    print >> client_f_to, 'error', response_arg
                    break
                assert response == 'found'
                p_response_arg = pickle.dumps(response_arg)
                print >> client_f_to, 'found', len(p_response_arg)
                client_f_to.write(p_response_arg)
                client_f_to.flush()
        except socket.error, e:
            info("closing connection on error [%s]", e)
        except EOFError, e:
            info("closing connection on EOF [%s]", e)
        except Exception, e:
            try: print >> client_f_to, ('error [%s "%s"]'
                                        % (sys.exc_info()[0], e))
            except: pass                # FIX: use more specific except clauses
            raise
        finally:
            try: client_f_to.close()
            except: pass
            try: client_f_from.close()
            except: pass
            try: client.close()
            except: pass
            signal.alarm(0)
            try: shutil.rmtree(state['work_directory'], ignore_errors=True)
            except: pass
            reap_child_processes()
            reset_state(state)

    info("exiting")


if __name__ == '__main__':
#     try:
#         import psyco
#         psyco.full()
#         ###psyco.bind(formula_mass)
#         warning('psyco enabled')
#     except ImportError:
#         pass

    try:
        if '--profile' in sys.argv:
            assert False, "fix chdir problem first"

            import cProfile
            import pstats
            report_fn = "sequest-chase.prof.%s" % os.getpid()
            data_fn = report_fn + ".tmp"
            prof = cProfile.run('main()', data_fn)
            with open(report_fn, 'w') as report_f:
                try:
                    stats = pstats.Stats(data_fn, stream=report_f)
                    stats.strip_dirs()
                    stats.sort_stats('cumulative')
                    stats.print_stats(50)
                    stats.sort_stats('time')
                    stats.print_stats(50)
                    print "# profile report written to '%s'" % report_fn
                finally:
                    try:
                        os.remove(data_fn)
                    except:
                        pass
        else:
            main()
    except SystemExit:
        raise
    except:
        logging.exception("unhandled exception")
        sys.exit(1)
    finally:
        logging.shutdown()
