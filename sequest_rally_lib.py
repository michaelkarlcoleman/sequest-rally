# common sequest-rally functions and constants

##  sequest-rally, a cluster driver for SEQUEST, derived from greylag
##  Copyright (C) 2006-2009  Stowers Institute for Medical Research


from __future__ import with_statement


import logging; from logging import debug, info, warning
import stat
import os

try:
    from fadvise import posix_fadvise, POSIX_FADV_DONTNEED
except ImportError:
    def posix_fadvise(*args): pass
    POSIX_FADV_DONTNEED = None


VERSION = "0.2"


ELAPSED_TIME_PRECISION = 3

SYSTEM_LOG = '/var/log/sequest-rally.log'
SYSTEM_LOG_LEVEL = logging.INFO

def set_logging(options, also_to_system_log=False):
    """Set up logging at the specified level, optionally logging to a common
    system logfile as well.
    """
    log_level = logging.WARNING
    if options.quiet:
        log_level = logging.ERROR
    if options.verbose:
        log_level = logging.INFO
    if options.debug:
        log_level = logging.DEBUG

    username = os.getenv('USER', 'unknown')
    fmt = ('%(asctime)s ' + username
           + ' [%(process)d] %(levelname)s: %(message)s')
    datefmt = '%b %e %H:%M:%S'
    formatter = logging.Formatter(fmt, datefmt)

    root = logging.getLogger('')
    root.setLevel(min(log_level, SYSTEM_LOG_LEVEL))
    if options.logfile:
        ulog_handler = logging.FileHandler(options.logfile)
    else:
        ulog_handler = logging.StreamHandler()
    ulog_handler.setLevel(log_level)
    ulog_handler.setFormatter(formatter)
    root.addHandler(ulog_handler)

    if not also_to_system_log:
        return
    try:
        slog_handler = logging.FileHandler(SYSTEM_LOG)
    except IOError, e:
        warning("could not write to system log '%s' [%s]", SYSTEM_LOG, e)
        return
    slog_handler.setLevel(SYSTEM_LOG_LEVEL)
    slog_handler.setFormatter(formatter)
    root.addHandler(slog_handler)


def file_block_at_offset(f, offset, length=None):
    """Read a block at the specified position in file f.  If length is None,
    read rest of file."""
    f.seek(offset)
    if length is None:
        return f.read()
    r = f.read(length)
    assert len(r) == length
    return r


def file_contents(fn, offset=None, length=None, binary=False):
    """Return the entire contents of a file, or of a specified segment.  Open
    the file in binary mode if specified.

    >>> file_contents('/dev/null')
    ''
    """
    mode = 'r' if not binary else 'rb'
    with open(fn, mode) as f:
        if offset is not None:
            return file_block_at_offset(f, offset, length)
        if length is not None:
            s = f.read(length)
            assert len(s) == length
            return s
        s = f.read()

        # tell the OS it can discard these pages, because the aggregate size
        # of the spectrum files may be huge, and we want to avoid pushing more
        # important stuff out of memory
        if fn.endswith('.ms2'):
            posix_fadvise(f.fileno(), 0, 0, POSIX_FADV_DONTNEED)

        return s


def file_sha1(filename):
    """Return the (hex) SHA1 digest of the given file.

    >>> file_sha1('/dev/null')
    'da39a3ee5e6b4b0d3255bfef95601890afd80709'
    """
    try:
        import hashlib
    except:
        return "no checksum--libs missing"
    h = hashlib.sha1()
    h.update(file_contents(filename))
    return h.hexdigest()


def get_file_fingerprint(filename):
    """Return a 'fingerprint' of this filename.  The fingerprint is a
    (filename, size, mod_time) triple.  (The filename should probably be an
    absolute path.)

    >>> get_file_fingerprint('/dev/null')[:2]
    ('/dev/null', 0)
    """
    st = os.stat(filename)
    return (filename, st[stat.ST_SIZE], st[stat.ST_MTIME])
