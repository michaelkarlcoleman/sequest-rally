README for sequest-rally

This is an adaption of greylag's scheduler to run SEQUEST
(specifically, the Yates' labs version).

Because we also run SEQUEST at our site, we adapted greylag
to run this binary (as opposed to its native spectrum
matcher).  A number of additional features were added to
increase reliability of searching with a computational
cluster, and to reduce the running time of PTM searches.
The result seems to work well, and is currently in
production at SIMR.

============================================================

ORIGINAL REQUIREMENTS

* Must handle arbitrary errors on nodes, due to (for
  example) warewulf filesystem update problems, filesystems
  full or missing, etc.

  Jan 23 17:01:50 [840] WARNING: node0005:20081 gave 'error got error 'sh: error while loading shared libraries: libtermcap.so.2: cannot open shared object file: Permission denied'

  Jan 23 17:01:55 [840] WARNING: node0007:20078 gave 'error [<type 'exceptions.OSError'> "[Errno 13] Permission denied: '/tmp/sequest-chase558_Yf'"]'

* Must handle nodes repeatedly disappearing and reappearing
  during runs.

* Must deal reasonably with skew--substantial variation in
  search times between spectra.  In particular, need to
  dynamically split spectra during long PTM searches.

* Must handle 300MB sequence database.

* There needs to be a way for users to rapidly abort their
  jobs.  In particular, it would be very nice if Control-C
  did this.  Similarly, if their PC/etc crashes, the job
  should either properly keep running or abort, but not go
  into a "zombie" state where cluster time is being wasted.

* Should do something reasonable if user stops their job
  (with Control-Z or kill -STOP).  Ideally they should be
  able to continue it later without trouble.  Whatever
  happens, it must not block the cluster from use by other
  jobs.


============================================================

ANNOUNCEMENT

The program 'sequest-rally', which is a replacement for PDQ
('sequest-process') is ready for testing.  It has a number
of improvements that will make life easier for SEQUEST users
and administrators both.

I'm still adding features, but it seems good enough to use.
It may very well have a few lingering bugs, though I've
squashed all of the ones I can find.  Please try it out,
compare the results to those of PDQ, and tell me if you find
a way to break it.

============================================================

Here is the least that you need to know to get started:

1.  Run it like so:

    $ sequest-rally -v -l sequest.log *.ms2 &
    $ tail -F sequest.log

If you need to, you can kill it using the 'kill' command (or
Control-C if you didn't background it with the '&').
There's only one process per job, so nothing more to worry
about.

2.  It uses the new 'unify_sequest' (the one for the
Orbitrap) by default.  There doesn't seem to be any good
reason to use the old one, but if you really want to, you
can add the '--old-sequest' flag.

3.  It does spectrum splitting.  For example, spectrum #1
might be searched against the first half of the sequence
database on one node, and against the second half of the
database on another node, with the results being combined in
a manner similar to what 'sqt-merge' does.

This means that the output will be slightly different than
the normal output, though it should never be of lower
quality.  The exception is that SpRank information is
unfortunately lost; ignore it for split spectra.  The output
will also vary a little from run to run, because split
points will vary according to the dynamics of the run.
(More on this below.)

If you need to avoid this, you can add the '--no-fragment'
flag.  The downside of doing so is that your search may take
longer to complete, particularly if it is a PTM search.
With this flag, the results should be essentially identical
to those from PDQ.

============================================================


DETAILS
-------

PDQ, though an improvement over what came before, is causing
some inefficient use of our cluster resources, not to
mention user frustration.  To address these issues, the
greylag scheduler has been adapted as a SEQUEST scheduler,
sequest-rally.  Like greylag, it is structured as a "master"
program (sequest-rally) that talks to multiple "slave"
programs (sequest-chase, one per CPU).  If there are
multiple masters, they contend for resources in a natural
way.

Features:

- Immediate abort.  If a sequest-rally is killed (or
  crashes), it immediately releases the cluster resources it
  was using, making them available rapidly for other jobs to
  use.  (PDQ jobs sometimes continue on for many hours,
  wasting cluster resources.)

- Resilience.  Like greylag, sequest-rally is not bothered
  by cluster nodes that crash, hang, act strangely, etc.  As
  long as the 'sequest-rally' process continues, it will
  attempt to complete the search.  (PDQ itself sometimes
  gets strung up indefinitely under these conditions.)

- Much better cluster saturation.  This is accomplished by
  general optimizations, spectrum ordering, and spectrum
  splitting.  (PDQ's saturation is sometimes good and
  sometimes not very good at all.)

- Completion email.  Under most circumstances, sequest-rally
  will send an email upon termination, letting you know that
  your job is done and whether or not it succeeded.  (This
  won't happen if the master node crashes, or if
  sequest-rally is killed with "kill -9", etc.)

- CPU and RAM usage seems to be very low (lower than PDQ).
  This should help prevent the cluster from getting swamped.

- Job pausing.  You can stop a 'sequest-rally' process with
  'kill -STOP <pid>' (or Control-Z if it is in the
  foreground).  After a minute or two, its chasers will
  disconnect and be available for other sequest-rally
  processes to grab.  You can continue (un-pause) with 'kill
  -CONT <pid>' (or the shell's 'fg' command to put it in the
  foreground).  Among other things, this could be used to
  let others "pass" you.


Sharing between multiple users
------------------------------

When sequest-rally starts, it will attempt to grab all idle
chasers.  If there are none, it will keep looking,
approximately once per minute.  As chasers become available,
they will be grabbed by one of the running sequest-rally
processes.  The result of this is that if there is just one
job, it gets all of the resources; if there are many, they
share the resources in a somewhat arbitrary way.

To keep long-running jobs from monopolizing the cluster, the
chasers will disconnect after two hours, at which point they
will be split arbitrarily by all running sequest-rally
processes.  This should keep long PTM jobs from completely
starving shorter jobs.

In addition, 'sequest-rally' has a '-u' (or
'--host-utilization') flag.  If given, this sets a ceiling
on the proportion of the cluster that will be used.  For
example, '-u 0.5' would specify that at most half of the
cluster CPUs will be used (even if there are more CPUs
idle).  This might be appropriate for PTM searches expected
to run for days.

We can look at other adjustments if these sharing methods
are not sufficient.


Spectrum splitting
------------------

The specific motivation for spectrum splitting is to
eliminate the situation common during PTM searches where you
find yourself waiting for hours for that "one last spectrum"
to finish.  This problem, sometimes referred to as "skew",
is common in parallel processing situations.

Our solution to the problem is to split the search of
individual spectra into pieces.  These pieces are then
searched in parallel and the results combined at the end (as
with sqt-merge).  Specifically, the sequence database is
divided into parts and each search is performed against one
of the parts.

Because of the SEQUEST scoring algorithm and because we
cannot modify SEQUEST, there are some limitations to this
approach.  The most significant problem is that SpRank
information is lost.  Unfortunately there is simply not
enough information in the SEQUEST output to allow correct
combined SpRank rank information to be generated from the
parts.  (This is also true for sqt-merge'd files).  Since we
(apparently) don't make much use of SpRank information, this
is hopefully not too much of an issue.

Beyond that, the peptides that a split search finds will
sometimes be *better* (e.g., higher XCorr) than the ones
found in a traditional search.  The cause of this is
SEQUEST's two-stage scoring strategy.  Specifically, it
scores all candidate peptides using a cheap initial score
function ("Sp"), then takes the top 500 and scores those
using a more expensive score function ("XCorr").  When a
spectrum is searched in parts, in addition to the 500
candidates that would have been scored in a tradition
search, many other candidates may end up being evaluated
with the expensive function.  The resulting top-ranked
peptides may be better than what would have been found, or
they may be the same, but they should never be worse.
(Spectrum splitting with greylag, when implemented, will not
suffer from this problem, because of its single-stage
scoring algorithm.)

Finally, some differences are due to ties that are split
arbitrarily, and which may be decided differently in split
searches.  In addition, one difference is due to a SEQUEST
bug: the comparison count overflows, but this happens less
often in a split search.  These differences are unimportant.

As currently implemented, the sequence database is only
split down to (at most) single loci.  For databases having
many loci this is not a problem, but if there are only a few
loci, this will potentially limit the granularity of
splitting.  Splitting a locus is simple in principle, but
would be difficult to do efficiently with unmodified
SEQUEST.  The splits would have to overlap by the length of
the longest peptide that might be found, which would lead to
quite a bit of duplicated search.  (This will not be a
problem with greylag.)  For now, it appears that this is
good enough.

Spectrum splitting, though generally beneficial, is not
completely free.  There is some overhead involved in setting
up the splits, and some extra computation is performed
because of the additional peptides evaluated with the
expensive scoring function, as mentioned above.
Nonetheless, although a bit more CPU is used, elapsed time
is usually shorter.

The specific time advantage of spectrum splitting depends on
the length of the longest spectrum search time with respect
to other search times, and with respect to the length of the
entire search.  If the longest spectrum takes three hours to
search (unsplit), then the elapsed time for the whole search
must necessarily be at least three hours if no splitting is
performed.  Even if the whole search takes much longer, if
one unluckily schedules this spectrum as the last one
started, the search may run about three hours longer than it
otherwise might have.

Taken together, these insights imply that PTM searches
(which have spectrum search times with high variation) will
benefit from spectrum splitting.  (I've already seen 10x
speedups.)  Non-PTM searches may or may not benefit much,
but probably will not usually be harmed.  More testing is
needed on that point.

As mentioned above, if desired, spectrum splitting can be
skipped with the '--no-fragment' flag.


Spectrum ordering
-----------------

Our studies of spectrum search times suggest that log(search
time) is highly correlated with spectrum mass.  Combining
this insight with a greedy scheduling heuristic leads to a
strategy of searching spectra in order of decreasing mass.
This also helps reduce elapsed time.  (PDQ searches in scan
order, which is often worse than random order.)


Other details
-------------

- The algorithm generally tries to send work to the chasers
  in two-minute chunks.  This limits communication overhead
  while allowing for reasonably up-to-date progress reports.
  During transitions (a job completes, a node crashes), it
  may take a minute or two until everything catches up.

- Just like 'sequest-process', the 'sequest-rally' program
  uses a lockfile named '.lock' to prevent multiple searches
  from being run simultaneously in the same directory (which
  would lead to confusing results).  Usually when you get an
  error mentioning the lock file, it means that there is
  already a search in progress (possibly on a different
  cluster).  If not, you can remove the lockfile manually,
  like so:

    $ rm -f .lock

- The search times that SEQUEST reports in 'S' lines are
  sometimes incorrect (negative) because of a SEQUEST bug.
  These incorrect values are corrected by sequest-rally.

- sequest-rally tries *really* hard to keep going with the
  search, even when things go wrong.  Sometimes it won't
  stop even in circumstances under which there is no hope of
  finishing.  Generally in cases like this you will see a
  stream of error messages in the log and no progress being
  made for long periods.  If you encounter a situation like
  this, please report it.

  That notwithstanding, the results should be correct unless
  a fatal error message ("ERROR...") was given.  The SQT
  files are not written until the very end, so if they were
  written, they are probably correct.

- See 'sequest-rally --help' to see what all of the flags
  do.
