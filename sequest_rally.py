#!/usr/bin/env python

'''This program coordinates the search of mass spectra by sending subsearches
to running instances of sequest-chase, which are expected to be listening on
the specified hosts/ports.

'''

from __future__ import with_statement

__copyright__ = '''
    sequest-rally, a cluster driver for SEQUEST, derived from greylag
    Copyright (C) 2006-2009  Stowers Institute for Medical Research
'''


import asynchat
import asyncore
from collections import defaultdict
import ConfigParser
import cPickle as pickle
import errno
import exceptions
import heapq
import itertools
import logging; from logging import debug, info, warning
import math
import optparse
import os
from pprint import pprint, pformat
import random
import re
import shutil
import signal
import smtplib
import socket
import subprocess
import sys
import tempfile
import time

import sqt_header

from sequest_rally_lib import *


# gc possibly harms performance here, so disable it.  gc only matters for
# cycles, which we (hope we) don't create.  See the gc module docs.
import gc; gc.disable()

START_TIME = time.time()

# FIX: put this in a config file
EMAIL_DOMAIN = 'stowers.org'
ADMIN_USER = 'mkc'
EMAIL_ADDRESSES = ['%s@%s' % (os.getenv('USER', 'nobody'), EMAIL_DOMAIN),
                   '%s@%s' % (ADMIN_USER, EMAIL_DOMAIN) ]

# if True, try to send email message to user on termination
mail_notification = False


def send_mail(successful, s=''):
    if not mail_notification:
        return
    from_ = 'nobody@%s' % EMAIL_DOMAIN

    debug("attempting to mail status to %s", EMAIL_ADDRESSES)

    body_prefix = ("After %.1f minutes, your command\n\n\t%s\n\n"
                   "on %s in directory\n\n\t%s\n\n"
                   % (((time.time() - START_TIME) / 60), ' '.join(sys.argv),
                      socket.gethostname(), os.getenv('PWD', '[unknown]')))

    if successful:
        status = "succeeded"
        body = ("%shas completed successfully.\n\n%s\n"
                % (body_prefix, s))
    else:
        status = "failed"
        body = ("%shas failed with error\n\n%s\n\n"
                % (body_prefix, s))

    server = smtplib.SMTP('localhost')
    message = ('From: %s\r\n'
               'To: %s\r\n'
               'Subject: %s\r\n'
               'Precedence: bulk'
               '\r\n'
               '%s\r\n'
               % (from_, ', '.join(EMAIL_ADDRESSES),
                  "your sequest-rally job %s" % status, body))
    try:
        server.sendmail(from_, EMAIL_ADDRESSES, message)
        server.quit()
    except smtplib.SMTPException, e:
        debug("sendmail failed: %s", e)


def error(s, *args):
    "fatal error"
    # if we're unit testing, just throw an exception
    if __name__ != "__main__":
        raise Exception((s + " (fatal error)") % args)
    logging.error(s, *args)
    send_mail(False, s % args)
    sys.exit(1)

# fatal for rally, but chase just disconnects
chase_error = error


# A lockfile is used to try to avoid collisions between multiple searches in
# the same directory.
LOCKFILENAME = '.lock'
lockfile_created = False


def read_fasta_file(f):
    """Yield (locusname, defline, sequence) tuples as read from the given
    FASTA file (uppercasing sequence)."""

    locusname, defline = None, None
    seqs = []

    for line in f:
        line = line.strip()
        if line[:1] == '>':
            if defline != None:
                yield (locusname, defline, ''.join(seqs))
            defline = line[1:]
            locusname_rest = defline.split(None, 1)
            locusname = locusname_rest[0] if locusname_rest else ''
            seqs = []
        else:
            if defline == None:
                chase_error("bad format: line precedes initial defline"
                            " in '%s'", (f.name if hasattr(f, 'name')
                                         else 'unknown FASTA file'))
            seqs.append(line.upper())
    if defline != None:
        yield (locusname, defline, ''.join(seqs))


def read_fasta_files(filenames):
    """Yield (locusname, defline, sequence, filename) tuples as read from
    FASTA files (uppercasing sequence).  An error is given if locusname is
    empty or not unique across all sequences.

    >>> list(read_fasta_files(['/dev/null']))
    []
    """

    loci_seen = set()

    for filename in filenames:
        with open(filename) as f:
            for locusname, defline, sequence in read_fasta_file(f):
                if not locusname:
                    chase_error("empty locus name not allowed in '%s'",
                                filename)
                if locusname in loci_seen:
                    chase_error("locus name '%s' is not unique in the search"
                                " database(s) in '%s'", locusname, filename)
                loci_seen.add(locusname)

                yield (locusname, defline, sequence, filename)


def effective_sequence_length(sequence):
    """Estimate the effective length of sequence, which accounts for ends and
    invalid residues (which are treated as ends).  So, for example, a sequence
    of length 25 will have somewhat fewer than half the number of peptides
    that a sequence of length 50 does.

    >>> effective_sequence_length('')
    0
    >>> effective_sequence_length(30*'A' + 'X' + 8 * 'V')
    15
    """

    AA_SEQUENCE = re.compile(r'[ARNDCQEGHILKMFPSTWYV]+')

    # If we subtract this from the length of all runs, we hope the result will
    # be proportional to the number of peptides in the run
    # (just a guess for now)
    END_CORRECTION = 15

    return sum(max(0, len(m.group()) - END_CORRECTION)
               for m in AA_SEQUENCE.finditer(sequence))


def create_sequence_clump_index(sequence_db_fn, max_clumps):
    """Given a sequence database filename, return the pair (list of clump
    start indices (database file offsets), effective database size).  These
    indicate how the database is broken into nearly equal sized contiguous
    sequence clumps.
    """
    sequence_sizes = [ effective_sequence_length(sequence)
                       for locusname, defline, sequence, filename
                       in read_fasta_files([sequence_db_fn]) ]
    effective_db_size = sum(sequence_sizes)
    clump_limit = max(max(sequence_sizes), effective_db_size / max_clumps)
    #debug("clump seq limit/sizes: %s %s", clump_limit, sequence_sizes)
    index = [0]
    current_size = 0
    for n, s in enumerate(sequence_sizes):
        current_size += s
        if current_size > clump_limit:
            current_size = 0
            index.append(n)
    # FIX: could try to equalize index tail (e.g., [100, 50, 50, 1])

    sequence_offsets = [ m.start() for m
                         in re.finditer(r'^>',
                                        file_contents(sequence_db_fn,
                                                      binary=True),
                                        re.MULTILINE) ]
    assert len(sequence_offsets) == len(sequence_sizes)

    return [ sequence_offsets[i] for i in index ], effective_db_size


reserved_fds = []

def reserve_fds():
    """Open some placeholder files.  This kludge keeps holds some low-numbered
    file descriptors available after connections to the chasers have been
    opened.  This is necessary so that we can use the subprocess module, which
    uses select, which can only use low-numbered file descriptors.  (This
    Python issue will be fixed in a future version.)
    """
    global reserved_fds
    assert not reserved_fds
    reserved_fds = [ open('/dev/null') for i in range(20) ]

def release_fds():
    "Close the placeholder files."
    global reserved_fds
    for f in reserved_fds:
        f.close()
    reserved_fds = []


def generate_file_lines_in_order(sort_args, fn, error_tag, force=False):
    """Yield lines from a sorted version of the specified file.  The file is
    first sorted with GNU sort, using the specified arguments.  The error tag
    is used in any diagnostics (to distinguish the source of the problem).

    >>> list(generate_file_lines_in_order(['-k1,1'], '/dev/null', 'test'))
    []
    """

    outfd, outfn = tempfile.mkstemp(prefix='sequest-rally.sort.tmp-')
    try:
        try:
            release_fds()
            command = ['sort'] + sort_args + ['-o', outfn, fn]
            p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            outs, errs = p.communicate()
            os.close(outfd)                 # just wanted unique fn
            if outs or errs or p.returncode != 0:
                error("sort failed: %s", error_tag)
            del p               # try to force resource cleanup (?)
        except ValueError, e:
            error("sort failed: %s", e)
        finally:
            reserve_fds()

        with open(outfn) as outf:
            for line in outf:
                yield line
    finally:
        try:
            os.unlink(outfn)
        except OSError, e:
            pass


def generate_descriptors_from_ms2_block(spectrum_fns):
    """Read ms2 files, yielding tuples (file_index, index, descriptor), in
    decreasing mass order, where descriptor is
    'example.ms2:234:host:00002.00002', where 234 is the byte offset of the
    start of this spectrum in example.ms2.  The indexes are 0-based,
    correspond to the original file order, and 'index' begins at zero for each
    file.

    >>> list(generate_descriptors_from_ms2_block(['/dev/null']))
    []
    """

    massfile = tempfile.NamedTemporaryFile(prefix='sequest-rally.mass.tmp-')
    for file_index, spfn in enumerate(spectrum_fns):
        # binary because we need accurate offsets
        contents = file_contents(spfn, binary=True)

        for index, m in enumerate(re.finditer(r'^:([^.]+\.[^.]+)\.[0-9]+\s+([0-9.]+)\s',
                                              contents, re.MULTILINE)):
            print >> massfile, m.group(2), file_index, index, \
                  '%s:%s:host:%s' % (spfn, m.start(), m.group(1))
    massfile.flush()

    for line in generate_file_lines_in_order(['--stable', '-k1,1nr'],
                                             massfile.name, "ms2 masses"):
        fs = line.rstrip().split(None, 3)
        debug("generating sp having mass = %s", fs[0])
        yield (int(fs[1]), int(fs[2]), fs[3])

    massfile.close()


def generate_descriptors_from_ms2_files(spectrum_fns, block_size=1000000000):
    """Read ms2 files, yielding tuples (file_index, index, descriptor,
    block_number), in decreasing mass order within blocks, where descriptor is
    'example.ms2:234:host:00002.00002', where 234 is the byte offset of the
    start of this spectrum in example.ms2.  The indexes are 0-based,
    correspond to the original file order, and 'index' begins at zero for each
    file.

    The spectrum files are partitioned into groups called blocks, so that each
    block is smaller than block_size (bytes) or has only one file.  These are
    numbered by block_number.

    >>> list(generate_descriptors_from_ms2_files(['/dev/null']))
    []
    """
    blocks = limitsplit(block_size,
                        ((os.path.getsize(fn), fn) for fn in spectrum_fns))
    fns_done = 0
    index = 0
    for block_no, b in enumerate(blocks):
        info("generating spectrum block %s", block_no)
        block_fns = [ p[1] for p in b ]
        for fn_index, index, desc \
                in generate_descriptors_from_ms2_block(block_fns):
            yield (fn_index+fns_done, index, desc, block_no)
        fns_done += len(b)


def check_and_count_spectra_in_ms2(spectrum_fn):
    """Return the number of spectra in the file.  Also, check the file format
    for general correctness and issue a fatal error if any problem is found.

    >>> check_and_count_spectra_in_ms2('/dev/null')
    0
    """

    # NB: We stipulate here that each spectrum may include just one colon,
    # an assumption that spectrum counting (below) depends on.
    ms2_spectrum_re = re.compile(r'''
          (?:
           :[^:]*\n
           [1-9]\d*\.\d+\ [1-9]\d*\r?\n
          )+
          (?:
           \d+\.\d+\ [0-9.]+\r?\n
          )+
                                  ''', re.VERBOSE)

    with open(spectrum_fn) as f:
        contents = f.read()

    prev_end = 0
    count = 0

    for m in re.finditer(ms2_spectrum_re, contents):
        if m.start() != prev_end:
            error("bad ms2 file '%s' beginning near offset %s of %s",
                  spectrum_fn, prev_end, len(contents))
        prev_end = m.end()
        count += 1

    if len(contents) != prev_end:
            error("bad ms2 file '%s' beginning near offset %s of %s",
                  spectrum_fn, prev_end, len(contents))

    # each match above may contain one or more actual spectra; explicitly
    # count them here
    return contents.count(':')


def count_ms2_spectra(spectrum_fns):
    """Return count of ms2 spectra, for ETA calculation, etc.

    >>> count_ms2_spectra(['/dev/null'])
    0
    """
    return sum(check_and_count_spectra_in_ms2(fn) for fn in spectrum_fns)


def get_sequest_param_value(param_name, params_contents):
    """Return stripped value of this sequest.params parameter, warning on
    duplicates.
    """
    v = re.findall(r'^%s\s*=\s*(.+)$' % param_name, params_contents,
                   re.MULTILINE)
    if len(v) > 1:
        error("multiple '%s' lines in 'sequest.params'")
    if not v:
        return '(none)'
    return re.split(r';|#', v[0])[0].strip()


def summarize_search_parameters(params):
    """Give the contents of sequest.params, return a summary line.
    Also check validity of mods.
    """

    mods = get_sequest_param_value('diff_search_options', params)
    m = re.match(r'^(([0-9]+(\.[0-9]+)?\s+[A-WYZ]+)\s+|(0(\.0+)?\s+X\s+)){3}$',
                 mods + ' ')
    if not m:
        error("invalid diff_search_options in sequest.params")

    ps = [ get_sequest_param_value(s, params)
           for s in ('diff_search_options', 'max_num_differential_AA_per_mod',
                     'enzyme_number', 'peptide_mass_tolerance',
                     'ppm_peptide_mass_tolerance') ]
    return 'mods %s depth %s enz %s tol_da %s tol_ppm %s' % tuple(ps)


def predsplit(predicate, iterable):
    """Yield sublists of iterable for which predicate is true of the initial
    item.
    >>> list(predsplit(lambda x: x % 2 == 0, [2, 3, 4]))
    [[2, 3], [4]]
    >>> list(predsplit(lambda x: x % 2 == 0, [1, 2, 5, 6, 6, 7, 7]))
    [[1], [2, 5], [6], [6, 7, 7]]
    >>> list(predsplit(lambda x: x % 2 == 0, []))
    []
    """
    group = []
    for i in iterable:
        if predicate(i):
            if group:
                yield group
            group = [i]
        else:
            group.append(i)
    if group:
        yield group


def pairpredsplit(predicate, iterable):
    """Yield sublists of iterable for which predicate is true of the final and
    initial item.
    >>> list(pairpredsplit(lambda x, y: x > y, [2, 3, 4, 1, 2]))
    [[2, 3, 4], [1, 2]]
    >>> list(pairpredsplit(lambda x, y: x < y, [2, 3, 4, 1, 2]))
    [[2], [3], [4, 1], [2]]
    >>> list(pairpredsplit(lambda x, y: x > y, [2]))
    [[2]]
    >>> list(pairpredsplit(lambda x, y: x > y, []))
    []
    """

    group = []
    for i in iterable:
        if group:
            if predicate(group[-1], i):
                yield group
                group = []
        group.append(i)
    if group:
        yield group


def limitsplit(limit, iterable):
    """Yield lists that partition iterable, such that each part is as long as
    possible, subject to the condition that the sum over the first elements is
    less than limit.  Each list contains at least one element, even if that
    element is over the limit.

    >>> list(limitsplit(100, [(75, 'a'), (50, 'b', 'bb'), (125, 'c')]))
    [[(75, 'a')], [(50, 'b', 'bb')], [(125, 'c')]]
    >>> list(limitsplit(100, [(75, 'a'), (10, 'b', 'bb'), (125, 'c')]))
    [[(75, 'a'), (10, 'b', 'bb')], [(125, 'c')]]
    >>> list(limitsplit(100, [(75, 'a'), (125, 'c'), (50, 'b', 'bb')]))
    [[(75, 'a')], [(125, 'c')], [(50, 'b', 'bb')]]
    >>> list(limitsplit(100, [(125, 'c'), (50, 'b', 'bb')]))
    [[(125, 'c')], [(50, 'b', 'bb')]]
    >>> list(limitsplit(100, []))
    []
    """

    s, p = 0, []
    for x in iterable:
        v = x[0]
        if s + v > limit:
            if p:
                yield p
            if v > limit:
                yield [x]
            else:
                s, p = v, [x]
                continue
            s, p = 0, []
        else:
            s += v
            p.append(x)
    if p:
        yield p


def parse_sqt_results(r):
    """Given an SQT results string, return
            (lines, S_fields, [(M_fields, [L_fields,...]), ...])
    where lines is the newline-split of the input and *_fields have been
    tab-split into fields, but not stripped.

    """

    lines = [ l for l in r.split('\n') if l ]
    S_fields = lines[0].split('\t')
    assert len(S_fields) == 10, "bad SQT format (wrong S field count)"
    assert S_fields[:1] == ['S'], ("bad SQT format (bad prefix) [%s]"
                                   % S_fields)

    ML_result = []
    for ML_lines in predsplit(lambda x: x[0][:1] == 'M',
                              [ l.split('\t') for l in lines[1:] ]):
        assert ML_lines[0][0] == 'M', "bad SQT format (M line not found)"
        ML_result.append((ML_lines[0], ML_lines[1:]))
    return (lines, S_fields, ML_result)


def ssum(s1, s2, places=0):
    """
    >>> ssum('123.111', '5.222', places=0)
    '128'
    >>> ssum('123.111', '5.222', places=3)
    '128.333'
    """
    return '%.*f' % (places, float(s1) + float(s2))


def M_fields_sort_key(M_fields):
    """Returns (Xcorr as float, peptide without flanks)
    >>> M_fields_sort_key((['M', '  1', ' 54', '3143.86878', '0.0000',
    ...                    ' 2.250', '  2.884', ' 22', '100',
    ...                    '  T.LDESVFRINMLIAMLINKRLRHLKFA.S', 'U'], [['L',]]))
    (-2.25, 'LDESVFRINMLIAMLINKRLRHLKFA')
    """
    return (-float(M_fields[0][5]), M_fields[0][9].split('.')[1])


def uniqify_list_of_lists(lol):
    """Return a sorted, uniqified list of the given list's sublists.
    >>> uniqify_list_of_lists([[1,2], [3,4], [1,2]])
    [[1, 2], [3, 4]]
    """
    return sorted(list(t) for t in set(tuple(l) for l in lol))


def flatten(gen_of_gens):
    """Yield items from a generator of generators.  (Only the top level is
    flattened.)

    >>> list(flatten([(1,2,3), (4, (5, 'b'), 6)]))
    [1, 2, 3, 4, (5, 'b'), 6]
    """
    for gen in gen_of_gens:
        for item in gen:
            yield item


def generate_merged_ML_groups(ML_fields):
    """For identical adjacent M lines, merge L lines.  Also, update score
    ranks.
    """

    def M_not_equal(m0, m1):
        "not identical, except that peptide flanks may differ"
        assert m0[0][0] == m1[0][0] == 'M'
        return (m0[0][3] != m1[0][3]    # mass
                or m0[0][5] != m1[0][5] # score
                or m0[0][9].strip()[2:-2] != m1[0][9].strip()[2:-2]) # PEPTIDE

    prev_score = None
    rank = 0
    for same in pairpredsplit(M_not_equal, ML_fields):
        # To merge M lines, choose the one with the smallest SpRank, as a best
        # guess.  (FIX: why SpRank?  better most tryptic instead?)
        best_M = min((ml[0] for ml in same), key=lambda x: int(x[2]))
        best_M_score = best_M[5]
        if best_M_score != prev_score:
            prev_score = best_M_score
            rank += 1
        # correct rank (replicating bogus SEQUEST whitespace, yuck)
        best_M = [ best_M[0], '  %s' % rank ] + best_M[2:]
        yield (best_M,
               uniqify_list_of_lists(flatten(x[1] for x in same)))


def ML_update_delta(best_score, M_fields):
    """Update the delta in this M line, using the given best_score.
    >>> ML_update_delta(2.0, ['M', '2', '2', '', '0.0000', '1.5000', '', '',
    ...                       '', '', ''])
    ['M', '2', '2', '', '0.2500', '1.5000', '', '', '', '', '']
    """
    if best_score > 0:                  # typical case
        delta = (best_score - float(M_fields[5])) / best_score
    else:
        delta = 0.0
    old_delta = M_fields[4].rstrip()
    width = len(old_delta) - old_delta.index('.') - 1
    return M_fields[:4] + ['%.*f' % (width, delta)] + M_fields[5:]


# FIX: get 'keep' from sequest.params
def merge_sqt_results(r0, r1, keep=5):
    """Given two search results (in SQT format), merge them, returning the
    best of both (according to the count to keep).
    """
    result = []

    # fs = list of tab-split fields; ML_fsx -> [ (M, [L,...]), ... ]
    lines0, S_fs0, ML_fs0 = parse_sqt_results(r0)
    lines1, S_fs1, ML_fs1 = parse_sqt_results(r1)
    assert S_fs0[:4] == S_fs1[:4], "results not from same spectra"

    S_fs = (S_fs0[:4] + [ ssum(S_fs0[4], S_fs1[4],
                               places=ELAPSED_TIME_PRECISION), 'host' ]
            + S_fs0[6:9] + [ ssum(S_fs0[9], S_fs1[9]) ])
    result.append('\t'.join(S_fs))

    # output best N M lines
    # note that sometimes there are no M lines at all
    ML_fs = sorted(ML_fs0 + ML_fs1, key=M_fields_sort_key)
    if ML_fs:
        ML_fs = list(generate_merged_ML_groups(ML_fs))
        best_score = float(ML_fs[0][0][5])
        best_ML_fs = ML_fs[:keep]

        for mlfs in best_ML_fs:
            result.append('\t'.join(ML_update_delta(best_score, mlfs[0])))
            for lfs in mlfs[1]:
                result.append('\t'.join(lfs))

    return '\n' + '\n'.join(result)


def write_sqt_files(spectrum_fns, result_files, result_file_indices):
    """Given result files and corresponding indices (one for each spectrum
    file), write results to sqt files, with spectra in the correct order.

    Each result file contains result strings, appended in arbitrary order,
    while each result index contains lines '<index> <offset> <length>', where
    offset and length indicate the position of the results for the
    spectrum_split in result_file.
    """
    for f in result_files + result_file_indices:
        f.flush()               # probably necessary

    # Writing these files directly to their final destination may be quite a
    # bit slower, because we're appending to all of them in parallel, so write
    # them to a temp directory (on a fast local disk) first, and then copy
    # them (in serial) into place.

    sqt_fns = [ os.path.splitext(os.path.basename(sfn))[0] + '.sqt'
                for sfn in spectrum_fns ]

    sqt_fs = [ tempfile.TemporaryFile(prefix='sequest-rally.tmp-sqt-')
               for fn in sqt_fns ]

    header = sqt_header.header('sequest.params')

    for result_file, result_file_index, sqt_f in \
            zip(result_files, result_file_indices, sqt_fs):
        sqt_f.write(header)

        # FIX: replace by itertools.groupby?
        for sp_results in pairpredsplit(lambda x, y: x[0] != y[0],
            ([ int(v) for v in line.split() ] for line in
             generate_file_lines_in_order(['-k1,1n', '-k2,2n'],
                                          result_file_index.name,
                                          "final results"))):
            r = file_block_at_offset(result_file, sp_results[0][1],
                                     sp_results[0][2])
            for sp_result in sp_results[1:]:
                r1 = file_block_at_offset(result_file, sp_result[1], sp_result[2])
                r = merge_sqt_results(r, r1)
            sqt_f.write(r)
        sqt_f.write('\n')

    for f in result_files + result_file_indices:
        f.flush()               # reclaim disk space

    # write final result files, sync them, and force errors, if any
    info("writing result files")
    for sqtf, sqt_fn in zip(sqt_fs, sqt_fns):
        sqtf.flush()
        sqtf.seek(0)
        f = open(sqt_fn, 'w')
        shutil.copyfileobj(sqtf, f)
        sqtf.close()
        f.flush()
        os.fsync(f.fileno())
        f.close()

    dfd = os.open('.', os.O_RDONLY)
    os.fsync(dfd)
    os.close(dfd)


def check_fingerprint(fingerprint):
    """Error if any of the fingerprinted files have changed.
    """
    for fp in fingerprint:
        fn = fp[0]
        current_fp = get_file_fingerprint(fn)
        if current_fp != fp:
            error("search input file '%s' has changed! [%s != %s]",
                  fn, current_fp, fp)


class chase_client(asynchat.async_chat):

    # heapq of (deathtime, host, port) for connections that have died or
    # timed-out.  deathtime is the time of death.  It's up to external code to
    # recreate after a suitable wait, if desired.
    dead_clients = []

    # set of (host, port) for connected clients
    connected_clients = set()

    # total number of chasers
    client_count = None

    # cached so that we don't have thousands of identical copies
    pickled_parameters = None

    # yields spectra to search
    spectrum_generator = None

    # if not None, this is the remaining portion of the previously yielded
    # spectrum that has not yet been sent to a chaser: (sp, start, end)
    spectrum_fragment = None

    # (host, port) -> (submittime, [ piece, ... ])
    # tracks the pieces each chase_client is currently searching
    #   where a piece is either an entire spectrum (spectrum, None)
    #                        or a (spectrum, (start, end))
    #   where start and end are indices into the database clump index
    #       (end is exclusive)
    pieces_in_progress = {}

    # heapq of (lastsubmittime, submitcount, piece) for pieces submitted (at
    # least once) and (possibly) not yet completed.  lastsubmittime is
    # provided so that we can submit the least recently submitted piece next.
    # submitcount is the number of times this piece has been submitted for
    # search.  Not populated until spectrum_generator has been exhausted (to
    # avoid keeping lots of spectra in memory).
    piece_queue = None

    # set of pieces already searched (only those searched since we populated
    # piece_queue)
    searched_pieces = None

    # number of clumps in sequence database
    db_clump_count = None

    # number of slivers for which results not yet received
    slivers_to_go = None

    # Adaptively try to adjust batch size so that chasers take about this long
    # (in seconds) to return with results.
    # FIX: this should depend on #hosts?
    batch_duration_goal = 120

    # This is incremented each time we start working our way through the
    # remaining pieces from the start of the list, or when we start a new
    # spectrum block.
    batch_generation = 0

    # This is the highest spectrum block id seen.  It's tracked so that we can
    # reset sliver_batch_size when we start a new block, which is necessary
    # because the new block's spectra may require much more time to search
    # than that of the old.
    spectrum_block_seen = 0

    # Sum of chaser time spent calculating results that ultimately were not
    # used.  (It's not really wasted, of course.)
    wasted_time = 0

    # If True, do not split the search of spectra.  This keeps the results
    # deterministic, at the possible price of more skew.
    do_not_fragment_spectra = False

    # files on which to write pickled results, one for each spectrum file
    result_files = None
    # indices for result_file, containing lines
    # <spectrum index> <offset> <length>
    result_file_indices = None

    @classmethod
    def set_parameters(self, parameters):
        assert self.pickled_parameters == None, "no change once set"
        self.pickled_parameters = pickle.dumps(parameters)

    @classmethod
    def set_result_files(self, result_files, result_file_indices):
        assert self.result_files == None, "no change once set"
        self.result_files = result_files
        self.result_file_indices = result_file_indices

    @classmethod
    def set_spectra(self, spectrum_generator, db_clump_count, slivers_to_go):
        """Fix the list of spectra to be searched, which must be done before
        any clients are created (and cannot be changed later).
        """
        assert self.spectrum_generator == None, "no change once set"
        self.spectrum_generator = spectrum_generator
        self.db_clump_count = db_clump_count
        self.searched_pieces = set()
        #self.searched_piece_counts = [0] * locus_count
        self.slivers_to_go = slivers_to_go


    def __init__(self, host, port):
        asynchat.async_chat.__init__(self)
        assert self.pickled_parameters
        assert self.spectrum_generator
        self.host = host
        self.port = port

        # adaptive count of slivers to send to our chaser each time; starts at
        # 1 and adjusts upwards to try to meet batch_duration_goal; reset to 1
        # at the beginning of each new generation
        self.sliver_batch_size = 1
        self.batch_generation_seen = self.batch_generation

        self.ibuffer = []
        self.set_terminator('\n')

        self.banner = None

        self.create_socket(socket.AF_INET, socket.SOCK_STREAM)
        self.set_keep_alive()
        debug("commence connect to %s:%s", host, port)
        self.connect((host, port))

    # override asynchat version to use a reasonable buffer size
    def push(self, data):
        self.producer_fifo.push(
            asynchat.simple_producer(data, buffer_size=self.ac_out_buffer_size))
        self.initiate_send()

    def __repr__(self):
        return "<%s connected to %s:%s>" % (self.__class__.__name__,
                                            self.host, self.port)


    def set_keep_alive(self):
        # try to better notice reboots, net failures, etc
        try:
            self.socket.setsockopt(
                socket.SOL_SOCKET, socket.SO_KEEPALIVE,
                self.socket.getsockopt(socket.SOL_SOCKET,
                                       socket.SO_KEEPALIVE) | 1
                )
        except socket.error:
            pass


    def _send_command(self, command, pickled_argument):
        self.push(command + ' ')
        self.push("%s\n" % len(pickled_argument))
        self.push(pickled_argument)


    # FIX: do better?
    def _adjust_batch_size(self, results):
        """adjust sliver_batch_size according to result search times,
        batch_duration_goal, slivers_to_go, and the number of chasers
        """
        total_time = sum(r[3] for r in results)
        slivers = sum((self.db_clump_count if r[1] is None
                       else r[1][1] - r[1][0]) for r in results)
        rate = slivers / total_time
        time_left = self.slivers_to_go / rate
        active_clients = len(self.connected_clients)
        bdg = self.batch_duration_goal
        if time_left / (active_clients + 0.01) < bdg:
            bdg /= 10.0
        self.sliver_batch_size = max(1, int(rate * bdg))
        debug("adjust: %s", (total_time, slivers, rate, time_left,
                             active_clients, bdg, self.sliver_batch_size))


    def _check_block(self, block):
        """If this spectrum block is newer than anything previously seen, bump
        batch_generation to restart batch size adjustment from scratch.
        """
        debug("check_block: %s %s [%s:%s]", block, self.spectrum_block_seen,
              self.host, self.port)
        if block > self.spectrum_block_seen:
            self.__class__.spectrum_block_seen = block
            self.__class__.batch_generation += 1
            debug("check_block: bumping generation to %s [%s:%s]",
                  self.batch_generation, self.host, self.port)


    def _submit_search(self):
        "send a search to client"

        if self.slivers_to_go <= 0:
            self.close()
            return

        submit_pieces = []

        if self.batch_generation > self.batch_generation_seen:
            self.sliver_batch_size = 1
            self.batch_generation_seen = self.batch_generation
        need = self.sliver_batch_size

        # first, check whether we were working on something and died
        if self.pieces_in_progress is not None:
            prev = self.pieces_in_progress.get((self.host, self.port))
            if prev:
                submit_pieces = prev[1]
                need = 0

        # get some never-searched spectra/slivers, if there are any
        # first use up the fragment, if any
        if need > 0 and self.spectrum_fragment:
            sp, f_sliver = self.spectrum_fragment
            f_start, f_end = f_sliver
            if f_end - f_start <= need:
                submit_pieces.append(self.spectrum_fragment)
                self.__class__.spectrum_fragment = None
                need -= f_end - f_start
            else:
                f_start_1 = f_start + need
                submit_pieces.append((sp, (f_start, f_start_1)))
                self.__class__.spectrum_fragment = (sp, (f_start_1, f_end))
                need = 0

        if need > 0 and self.spectrum_generator != 'done':
            # If we need more than a half spectrum, grab a whole number of
            # spectra, for efficiency.  Otherwise, grab one spectrum and split
            # it.
            if self.do_not_fragment_spectra or need >= self.db_clump_count / 2:
                more = [ (sp, None) for sp
                         in itertools.islice(self.spectrum_generator,
                                             max(1, (need /
                                                     self.db_clump_count))) ]
                if more:
                    # FIX: this only affects subsequent submissions after a
                    # new block is detected, so the first might be overly
                    # large, wasting one chaser in the worst case.  Probably
                    # not much effect in practice, though.
                    self._check_block(more[-1][0][3])
                    submit_pieces.extend(more)
                    need = 0
            else:
                more = list(itertools.islice(self.spectrum_generator, 1))
                if more:
                    self._check_block(more[-1][3])
                    assert 0 < need < self.db_clump_count
                    submit_pieces.append((more[0], (0, need)))
                    assert not self.spectrum_fragment
                    self.__class__.spectrum_fragment = (more[0],
                                                        (need,
                                                         self.db_clump_count))
                    need = 0

        if submit_pieces:
            self.pieces_in_progress[(self.host, self.port)] = (time.time(),
                                                               submit_pieces)
        else:
            # spectrum generator exhausted
            if self.piece_queue == None:
                # move pieces_in_progress to piece_queue
                # note that heap entries are single spectra (or parts thereof)
                self.__class__.piece_queue = []
                for t, subsp in self.pieces_in_progress.values():
                    for sp in subsp:
                        heapq.heappush(self.piece_queue, (t, 1, sp))
                self.__class__.pieces_in_progress = None
                self.__class__.piece_generator = 'done'
                self.__class__.batch_generation += 1
                info("queued final %s pieces", len(self.piece_queue))

            # generate submit_pieces from piece_queue (without wrapping around)
            for _i in xrange(len(self.piece_queue)):
                if not self.piece_queue or need <= 0:
                    break

                _submittime, subcount, piece = heapq.heappop(self.piece_queue)
                if piece in self.searched_pieces:
                    continue            # already have search result
                heapq.heappush(self.piece_queue, (time.time(), subcount+1,
                                                  piece))
                submit_pieces.append(piece)
                if piece[1] == None:
                    need -= self.db_clump_count
                else:
                    need -= piece[1][1] - piece[1][0]

        assert submit_pieces
        debug("submitting %s slivers (%s pieces) to %s:%s (%s sl to go) %s",
              sum((self.db_clump_count if p[1] is None else p[1][1]-p[1][0])
                  for p in submit_pieces),
              len(submit_pieces), self.host, self.port, self.slivers_to_go,
              submit_pieces)
        self._send_command('search', pickle.dumps(submit_pieces))


    def _remember_result(self, result):
        """Write result to temp file, iff these pieces have not already been
        searched.  The variable slivers_to_go must be maintained exactly,
        because it is used to detect completion.
        """
        if self.pieces_in_progress is not None:
            pipkey = (self.host, self.port)
            if pipkey in self.pieces_in_progress:
                del self.pieces_in_progress[pipkey]

        sp_spec, clump_spec, out_s, elapsed_time = result
        piece = (sp_spec, clump_spec)

        debug("result %s[%s:%s] %s",
              ('(seen) ' if piece in self.searched_pieces else ''), self.host,
              self.port, result)
        if piece in self.searched_pieces:
            self.__class__.wasted_time += elapsed_time
            return                     # already have all results for piece
        self.searched_pieces.add(piece)

        if clump_spec is None:
            self.__class__.slivers_to_go -= self.db_clump_count
        else:
            self.__class__.slivers_to_go -= clump_spec[1] - clump_spec[0]
        assert self.slivers_to_go >= 0

        # remember result since it has not already been seen
        which_file = sp_spec[0]
        result_file = self.result_files[which_file]
        result_file_index = self.result_file_indices[which_file]
        offset = result_file.tell()
        result_file.write(out_s)
        length = result_file.tell() - offset
        print >> result_file_index, sp_spec[1], offset, length


    def _receive_response(self):
        r = ''.join(self.ibuffer)
        self.ibuffer = []
        return r


    def handle_connect(self):
        debug("connecting to %s:%s", self.host, self.port)

    def handle_expt(self):
        pass

    def handle_error(self):
        #self.handle_close()
        exc, why, _traceback = sys.exc_info()
        self.connected = False
        if (self.host, self.port) in self.connected_clients:
            self.connected_clients.remove((self.host, self.port))
        if exc == exceptions.KeyboardInterrupt:
            error("received keyboard interrupt")
        if exc == exceptions.SystemExit:
            error("received SIGTERM (SystemExit)")
        if exc == socket.error:
            if why[0] == errno.ECONNREFUSED:
                debug("no chaser at %s:%s (connection refused)", self.host,
                      self.port)
            else:
                debug("network error on connection %s:%s (%s %s)", self.host,
                      self.port, exc, why)
                debug("  traceback: %s", asyncore.compact_traceback()[3])
        else:
            info("unexpected error on connection %s:%s (%s %s)", self.host,
                 self.port, exc, why)
            info("  traceback: %s", asyncore.compact_traceback()[3])


    def handle_close(self):
        # FIX
        if (self.host, self.port) in self.connected_clients:
            self.connected_clients.remove((self.host, self.port))
        self.dead_clients.append((time.time(), self.host, self.port))
        asynchat.async_chat.handle_close(self)


    def collect_incoming_data(self, data):
        self.ibuffer.append(data)

    def found_terminator(self):
        if self.banner == None:
            self.banner = self._receive_response()
            banner_words = self.banner.split()
            if banner_words[0] != 'sequest-chase':
                warning("%s:%s is not a sequest-chase", self.host, self.port)
                self.handle_close()
                return
            assert banner_words[2:3] == ['ready']
            # now we can start talking
            self.connected_clients.add((self.host, self.port))
            self._send_command('parameters', self.pickled_parameters)
            self._submit_search()
            return

        if self.get_terminator() == '\n':
            reply = self._receive_response()
            reply_words = reply.split()
            if reply_words[0] == 'error':
                warning("%s:%s gave '%s'", self.host, self.port, reply)
                self.handle_close()
                return

            assert reply_words[0] == 'found' and len(reply_words) == 2
            self.set_terminator(int(reply_words[1]))
            return

        results = pickle.loads(self._receive_response())

        #debug("received: %s" % results)
        # FIX: check that multiply searched spectra had same results
        for result in results:
            self._remember_result(result)

        self._adjust_batch_size(results)

        self.set_terminator('\n')
        self._submit_search()


def signal_handler(signum, frame):
    error("received signal %s, exiting", signum)


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <ms2-file>...",
                                   description=__doc__, version=VERSION)
    pa = parser.add_option
    DEFAULT_HOSTFILE = "/etc/sequest-rally/hosts"
    pa("--hostfile", dest="hostfile", default=DEFAULT_HOSTFILE,
       help="file listing host:port locations where sequest-chase workers are"
       " listening.  [default '%s']" % DEFAULT_HOSTFILE, metavar="FILE")
    pa("-u", "--host-utilization", dest="host_utilization", type="float",
       default=1.0, help="proportion of available chasers to use"
       " [default=1.0, meaning all chasers]", metavar="PROPORTION")
    pa("--old-sequest", action="store_true", dest="old_sequest",
       help="use the old 'unify_sequest'; not suitable for Orbitrap searches")
    pa("--no-fragment", action="store_true", dest="no_fragment",
       help="do not fragment search of single spectra; this keeps the results"
       " deterministic, at the possible price of more skew")
    DEFAULT_BLOCK_SIZE = 1024 ** 3 / 2
    pa("--block-size", dest="block_size", type="int",
       default=DEFAULT_BLOCK_SIZE,
       help=("size of spectrum files subject to simultaneous reads [default=%s]"
             % DEFAULT_BLOCK_SIZE), metavar="BYTES")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("-q", "--quiet", action="store_true", dest="quiet", help="no warnings")
    pa("-l", "--logfile", dest="logfile",
       help="log to FILE instead of stderr", metavar="FILE")
    pa("--no-mail", action="store_true", dest="no_mail",
       help="do not send email on termination")
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    pa("--debug", action="store_true", dest="debug",
       help="output debugging info")
    pa("--profile", action="store_true", dest="profile",
       help="dump Python profiling output to './sequest-rally.prof.<pid>'")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    spectrum_fns = args

    if (any(not f.endswith('.ms2') for f in spectrum_fns)
        or not (0.0 < options.host_utilization <= 1.0)
        or options.block_size < 1):
        parser.print_help()
        sys.exit(1)

    if any(f.startswith(('tmp.', 'phs.', 'mth.')) for f in spectrum_fns):
        sys.exit("error: spectrum filename clashes with SEQUEST tempfile")

    if (options.logfile is not None
        and (options.logfile.endswith(('.ms2', '.sqt', '.params'))
             or options.logfile.startswith('-'))):
        sys.exit("error: cannot use '%s' as a log file" % options.logfile)

    # check for lockfile
    lock_contents = "%s on %s" % (os.getpid(), socket.gethostname())
    try:
        os.symlink(lock_contents, LOCKFILENAME)
    except OSError, e:
        pass

    lock_contents_read = 'unreadable'
    try:
        lock_contents_read = os.readlink(LOCKFILENAME)
    except OSError, e:
        pass
    if lock_contents_read != lock_contents:
        sys.exit("this directory locked (pid = %s); remove '%s' if not"
                 % (lock_contents_read, LOCKFILENAME))
    global lockfile_created
    lockfile_created = True

    set_logging(options, also_to_system_log=True)

    if not options.no_mail:
        global mail_notification
        mail_notification = True

    # try to notice/log/mail on SIGTERM (and maybe other signals)
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGHUP, signal.SIG_IGN)

    info("%s (%s) starting on %s", os.path.basename(sys.argv[0]), VERSION,
         socket.gethostname())
    info("command is %s", sys.argv)
    cwd = os.getcwd()
    info("cwd is '%s'", cwd)

    # try to check early that we can open result files
    # (this file is opened in binary mode)
    temp_result_files = [
        tempfile.TemporaryFile(prefix='sequest-rally.%s.tmp-' % fn)
        for fn in spectrum_fns ]
    temp_result_file_indices = [
        tempfile.NamedTemporaryFile(prefix='sequest-rally.%s-idx.tmp-' % fn)
        for fn in spectrum_fns ]

    # read chaser (host, port) listener list
    if not os.path.exists(options.hostfile):
        error("'%s' does not exist (not a cluster?)", options.hostfile)
    try:
        with open(options.hostfile) as hostf:
            hosts = [ l.split(':', 1) for l in hostf ]
            hosts = [ (host, int(port)) for host, port in hosts ]
    except ValueError:
        error("invalid or empty host line in '%s'" % options.hostfile)

    for host, port in hosts:
        if not (0 < port < 65536):
            error("invalid port number '%s' in hosts file", port)

    if not hosts:
        error("no valid search hosts specified")

    if options.host_utilization < 1.0:
        u_count = max(1, int(round(options.host_utilization * len(hosts))))
        random.shuffle(hosts)
        hosts = hosts[:u_count]
    hosts.sort()

    # check spectrum basename uniqueness, as corresponding sqt files will be
    # in a single directory
    base_spectrum_fns = [ os.path.basename(fn) for fn in spectrum_fns ]
    if len(base_spectrum_fns) != len(set(base_spectrum_fns)):
        error("base spectrum filenames must be unique")

    try:
        params_s = file_contents('sequest.params')
    except IOError, e:
        error("cannot read 'sequest.params' [%s]", e)
    dbs = re.findall(r'^database_name\s*=\s*(.+)', params_s, re.MULTILINE)
    if len(dbs) != 1:
        error("invalid 'database_name' entry in 'sequest.params'"
              " (none or more than one)")
    dbfn = dbs[0].strip()
    info("database is '%s'", dbfn)
    DBFN_MAXLEN = 80
    if len(dbfn) > DBFN_MAXLEN:
        error("SEQUEST cannot handle database_name longer than %s characters",
              DBFN_MAXLEN)

    db_checksum = file_sha1(dbfn)

    dbfn_contents = file_contents(dbfn)
    locus_count = sum(1 for m in re.finditer(r'^>', dbfn_contents,
                                             re.MULTILINE))
    del dbfn_contents

    # the '2' is intended to smooth things--does this matter?
    db_clump_index, effective_db_size \
        = create_sequence_clump_index(dbfn, 2*len(hosts))
    assert db_clump_index[0] == 0
    debug("database clump index: %s", db_clump_index)
    db_clump_count = len(db_clump_index)
    info("database has %s loci, %s clumps, %s effective residues",
         locus_count, db_clump_count, effective_db_size)

    # FIX: is it worth checking validity of ms2 files?
    spectrum_count = count_ms2_spectra(spectrum_fns)

    if not spectrum_count:
        error("no input spectra")
    info("%s spectra in %s files", spectrum_count, len(spectrum_fns))

    info(summarize_search_parameters(params_s))

    # abort if any of these absolute paths change during search
    watch_files = [ (fn if fn.startswith('/') else os.path.join(cwd, fn))
                    for fn in ['sequest.params', dbfn] + spectrum_fns ]
    watch_files_fingerprints = [ get_file_fingerprint(fn)
                                 for fn in watch_files ]

    # a sliver is a search of one spectrum against one clump
    sliver_count = db_clump_count * spectrum_count
    info("%s slivers to be searched", sliver_count)

    # now connect to chasers and loop until done (or timeout)
    command = "unify_sequest-1"
    if options.old_sequest:
        command = "unify_sequest"
    parameters = { 'pwd' : os.getcwd(),
                   'command' : command,
                   'sequence_db' : dbfn,
                   'sequence_db_checksum' : db_checksum,
                   'sequence_db_clump_index' : db_clump_index,
                   }

    chase_client.set_parameters(parameters)
    chase_client.set_spectra(generate_descriptors_from_ms2_files(spectrum_fns,
                                                          options.block_size),
                             len(db_clump_index), sliver_count)
    chase_client.set_result_files(temp_result_files, temp_result_file_indices)
    chase_client.client_count = len(hosts)
    if options.no_fragment:
        chase_client.do_not_fragment_spectra = True

    # adjust for performance (620000?)
    chase_client.ac_in_buffer_size = 62000
    chase_client.ac_out_buffer_size = 62000

    reserve_fds()               # reserve before opening chaser connections

    for host, port in hosts:
        chase_client(host, port)

    # retry time for dead clients (in seconds)
    retryafter = 60

    start_time = time.time()
    last_status_time = 0
    last_fingerprint_time = 0
    while True:
        asyncore.loop(count=1, timeout=retryafter, use_poll=True)
        if chase_client.slivers_to_go <= 0:
            break

        now = time.time()
        if now - last_status_time >= 10:
            message = ("%s to go, %s/%s/%s chasers"
                       % (chase_client.slivers_to_go,
                          len(chase_client.connected_clients),
                          len(asyncore.socket_map), len(hosts)))
            slivers_finished = sliver_count - chase_client.slivers_to_go
            if float(slivers_finished) / sliver_count >= 0.0001:
                eta_minutes = ((now - start_time)
                               / slivers_finished
                               * chase_client.slivers_to_go / 60.0)
                # guessed empirical correction (doesn't work very well)
                # eta_minutes = 5 * math.log(max(0.0001, eta_minutes))
                message += ", ETA %dm" % max(0, int(round(eta_minutes)))
            info(message)
            last_status_time = now

        if now - last_fingerprint_time >= 60:
            check_fingerprint(watch_files_fingerprints)
            last_fingerprint_time = now

        while (chase_client.dead_clients
               and chase_client.dead_clients[0][0] < now - retryafter):
            debug("died %s restart %s", chase_client.dead_clients[0][0], now)
            dt, host, port = heapq.heappop(chase_client.dead_clients)
            chase_client(host, port)
        if not asyncore.socket_map and chase_client.dead_clients:
            debug("sleeping")
            time.sleep(max(0, retryafter - (now
                                            - chase_client.dead_clients[0][0])))
        # failsafe, better performance?
        time.sleep(0.05)

    debug('searched pieces cache size was %s',
          len(chase_client.searched_pieces))

    # This excludes time spent by chasers still running at this point.  (FIX:
    # any good way to include this?)
    info("nonfinal nonproductive time: %.1fs", chase_client.wasted_time)

    # close all client sockets to immediately release chasers for other work
    for fd in asyncore.socket_map:
        try:
            os.close(fd)
        except OSError, e:
            info("close on fd %s failed", fd)

    info("preparing result files")
    write_sqt_files(spectrum_fns, temp_result_files, temp_result_file_indices)
    # FIX: need to explicitly close the index files here?

    check_fingerprint(watch_files_fingerprints) # be paranoid

    send_mail(True)
    info("finished (elapsed time %.1f minutes)"
         % ((time.time() - START_TIME) / 60))


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
            import cProfile
            import pstats
            report_fn = "sequest-rally.prof.%s" % os.getpid()
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
    except KeyboardInterrupt:
        error("received keyboard interrupt")
        send_mail(False, "received keyboard interrupt")
        sys.exit(1)
    except:
        logging.exception("unhandled exception")
        send_mail(False, "unhandled exception")
        sys.exit(1)
    finally:
        if lockfile_created:
            try:
                os.remove(LOCKFILENAME)
            except OSError, e:
                pass
        logging.shutdown()
