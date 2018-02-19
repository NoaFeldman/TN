# ==================================================================== #
# MPSPACK v3.0
# ==================================================================== #
  
  MPSPACK_v3.0pre (prelease) synced to $WBX : Feb 06, 2015

  For more recent changes,
  ==> see end of this readme file!

# -------------------------------------------------------------------- #
# new environmental variables with MPSPACK v3.0

  RC_STORE

    path to somewhere in /data where you would like to
    store all Clebsch Gordan related data
    (if your path ends with .../RCStore/ this last directory
    will be created automatically if it does not exist).

    Update 09/2016: by now, this path may contain multiple paths,
    e.g: export RC_STORE=PATH1:PATH2:...:PATHN
    with the following procedure: symmetry related data is looked
    up in PATH1 first; if not found, it proceeds to PATH2, etc.
    if not found all the way to PATHN, the required data is 
    computed and store in the LAST path PATHN. Therefore the paths
    PATH1 to PATH(N-1) are read-only.

    Background: RC_STORE gets large for larger SU(N) symmetries;
    therefore when working on a cluster, one does not want to
    copy the entire RC_STORE for every single job. Therefore
    when using large symmetries such as SU(N>3) [for SU(N<=3)
    one may as well just synchronize or regenerate symmetry
    related data], the recommended setting is:

       RC_STORE=NETWORK_PATH:LOCAL_PATH

    where NETWORK_PATH e.g. is some globally visible
    /project/theorie/.. network directory, whereas LOCAL_PATH
    e.g. is some job specific directory on /data/$USER/myjob/..
    This realizes differential storage, in that only new or
    updated data is stored locally. In cases where the RC_STORE
    on the NETWORK_PATH was already complete, the LOCAL_PATH
    will then remain empty.  Conversely, if significant new
    data was generated locally, one can update i.e. rsync the
    RCStore in the NETWORK_PATH from some LOCAL_PATH.
    It is recommended to do this after all jobs are finished
    to avoid interference with other running jobs.

    Note that instead of using /project/*, the NETWORK_PATH
    can also point to some LOCAL machine's /data/$USER/ ...
    This directory can then be mounted in the job directory using

        cd /data/$USER/myjob/...
        mkdir ./RCStore
        sshfs $USER@some_machine:/data/$USER/etc. ./RCStore
        run_my_program

    which should be unmounted once the job is finished

        fusermount -u ./RCStore

  CG_VERBOSE

    level of log output;
    general guideline on the value of CG_VERBOSE:

     = 0  basically log nothing concerning CG activity
    <= 4  log mainly information regarding newly generated RCStore data
     > 5  also log information on RCStore data that is (newly)
          read from RCStore.

# -------------------------------------------------------------------- #
# comments on itags (index tags)

  itag specficiation through cell array of strings.

  "conjugate" flags (= trailing '*') must be preserved.

  maximum string length for single specific itag is 8 characters
  (the reason being that internally, the itag is simply converted
  to unsigned long which has 8 chars)

# -------------------------------------------------------------------- #
# Notes regarding Clebsch Gordan tensors

# In general: do not even think about tinkering with the new
# info.cgr references to Clebsch Gordan spaces!

  e.g. must be careful about signs for degenerate q-labels
  which are left unchanged through permutation! e.g. for SU2:
  permute((1,1;0),[2 1 3]) acquires a minus sign!

# only exception for (small) CGC data during setup:
# set basic CData info fields (type,qset,qdir,cgw),
# while keeping the rest empty:

  c=emptystruct(A_someQSpace.info.cgr);
    c.type=qtype;
    c.qset=repmat(q,1,numel(Ak.Q));
    c.qdir=Ik.cgr(1).qdir;
    c.cgw=['1' 0];

  A_someother.info.cgr(:,j)=c;

# -------------------------------------------------------------------- #
# further comments regarding specific (new) routines

  @QSpace/getvac.m 

    allows to obtain vacuum state
    (through identity QSpace) for given symmetry setting,
    e.g. see setupSIAM_SU2x2.m: q=getvac(Z);

# -------------------------------------------------------------------- #
# removed environmental variable WB_CONTRACT_QFLAG
  and with it the contractQS( '-[lL]' options! // -l, -q
  => use automated lenient flag! // Wb,Nov22,14

# renamed routines for better naming consistency
  see $MEX/rename_qs.sh and archived versions therefrom

     for historical reasons
        MPSPACK_v1  => got QS label
        MPSPACK_v2  => got CG label

     mpsMakeUnique   => makeUniqueQS   # Wb,Nov22,14
     mpsNormQS       => normQS
     mpsPermuteQS    => permuteQS
     mpsPlusQS       => plusQS
     mpsSkipZerosQS  => skipZerosQS
     mpsGetDim       => getDimQS
     mpsGetQDim      => getQDimQS
     mpsIsIdentityCG => isIdentityCG
     mpsIsIdentityQS => isIdentityQS
     mpsTimesEl      => timesElQS
     mpsIsDiagQS     => isDiagQS

   * mpsMaxDiffQS    => maxDiffQS      # (mostly) unused
   * mpsGetDRange    => getDRangeQS    # (mostly) unused

     NRGWilsonCG     => NRGWilsonQS    # Wb,Nov24,14
     fdmNRG_CG       => fdmNRG_QS

     mpsOrthoQS      => orthoQS        # Wb,Aug25,15
     mpsGetSVD       => svdQS          # Wb,Aug25,15
     mpsEigenSymmQS  => eigQS          # Wb,Aug25,15

# -------------------------------------------------------------------- #
# Wb,May01,14
  outsourcing CGC's in QSpace into CGC library (single file)
  => not very flexible, thus replaced by the following:

# Wb,Dec15,14
  outsourcing CGC's in QSpace into CGC library (file system)
  => requires new environmental variable RC_STORE
  => since MANY files are generated, DO NOT store
     in backed up file systems such as HOME! => use /data

# -------------------------------------------------------------------- #
# Wb,Aug25,15 :: update / renamed routine

  mpsEigenSymmQS  => eigQS

  the behavior of this routine changed, in that in the presence
  of non-abelian symmetries the first return argument now contains
  a 2nd column which contains degeneracy (i.e. combined multiplet
  dimension)

# -------------------------------------------------------------------- #
# Wb,Aug24,15 :: reactivated / renamed routine

  mpsOrthoQS => orthoQS
  mpsGetSVD  => svdQS

  Here svdQS() is a specialization of orthoQS
  in that instead of returning [U,S,V] it returns [U,(SV')'=VS];
  Both routines accept truncation arguments (Nkeep, stol)

  Note that orthonormalization only makes sense in the presence
  of orthonormal input spaces! Hence the above routines now
  require in most cases that all indices are of the same type
  (i.e. all-in).

  The input index (idx) indicates which indices shall be 
  combined. Since non-abelian symmetries require proper multiplet
  fusion, here the routines orthoQS and svdQS have been 
  constrained such that idx eiter represents a single index,
  or (r-1) indices with r the rank (=dimensionality) of the 
  QSpace tensor.

# -------------------------------------------------------------------- #
# switched from matlab/2013a to matlab/2013b
# since it supports gcc-4.7.4 (-std=c++11)
# which now supports 'explicit' type conversion
# Wb,Jan04,16

  WRN! starting with matlab/2013b, flip() is also a matlab routine now!
  it has pretty much the same in syntax as my previous lib/flip.m which
  led to infinite recursive calls [e.g. my flip() calls flipud -> calls
  flip(x,1), which however again reroutes to my flip.m routine]
  Solution: simply removed my lib/flip.m routine.
  Wb,Jan04,16

  Drawback: 'clear all' appears to unload shared libraries
  => error message: shared object not open
  => requires to restart Matlab to return to regular operation

# -------------------------------------------------------------------- #
# altered treatment of decomposition of tensor-product of ireps
# Wb,Jan04,16

  further parallelized sparse tensor-product decomposition
  e.g. to speedup SU(6) calculations; this is mostly irrelevant
  in terms of speed for SU(N<=4) calculations.

  orthogonalization of inner multiplicity spaces is now 
  based on QR decomposition (the two-fold Gram-Schmidt used so far
  for numerical stability is numerically more expensive, yet still
  not as accurate anyway)

  => NB! this can lead to rotation in inner multiplicity spaces
     for SU(N>=3); therefore old RC_STORES need to be recomputed!

# -------------------------------------------------------------------- #
# minor adjustments for qlabel representation of non-abelian symmetries:
# Wb,Jan04,16

  For a more readable / compact representation of multiplet labels
  (q-labels) of non-abelian irreps, the numbers 10-35 are now replaced
  by the letters A-Z. E.g. the SU(4) multiplet (8 10 12) is 
  written in the compact form (8AC); similarly, this also affects
  SU(2) in that for example the multiplet q=2S=12 may be simply
  written as q=(C).

  This only affects log-outputs and file names in the RC_STORE.

# -------------------------------------------------------------------- #
# change in log behavior # Wb,Jan09,16

  log output referring to computation of generalized CGCs has
  been removed from STDOUT and are now written centralized
  with the RC_STORE; for larger CG_VERBOSE, the output to
  the log-file eventually is also replicated to STDOUT.

# -------------------------------------------------------------------- #
# I changed my default Matlab version from 2013b to to R2016a
# Major reason: severe SVD bug in Matlab 2013b .. 2015b
# Wb,Aug30,16

  Since in the meantime matlab also changed its graphics handles
  to actual class objects, many minor changes in scripts that
  produce graphics were incurred. Several plot scripts that
  I have used in the meantime, are updated. More to follow.

  I kept the old version in MPSPACK_v3.0pre/
  and moved the new version to MPSPACK_v3.0/

  Mostly affected directories

    bin/           please update all mex files
    util/          please update all mex files, etc.
    lib/           please change as far as necessary.
    Class/@QSpace
    Class/@SymOp

# -------------------------------------------------------------------- #
# Added Z_n symmetry
# this essentially acts like U(1) with charge q taken module n.
# Wb,Aug31,16

  New abelian symmetry labels that can be used to replace 'Acharge'
  (already also with getLocalSpace)

   * 'ZNcharge' with N=2,3,...
      total U(1) charge symmery reduced to a Z_N symmetry
      with q-labels = [U(1) charge modulo N] = {0,1,...,N-1}

   * 'Pcharge'  charge parity
      this is equivalent to 'Z2charge', except that symmetry
      labels are {-1,1} instead of {0,1}).

# -------------------------------------------------------------------- #
# Parallelization in QSpace // Wb,Feb07,17

  All of the parallelization in QSpace concerns threading
  in shared memory environments and is therefore based
  on the OpenMP standard. It is controlled via a set of
  environmental variables.

  Relevant environmental variables:

  OMP_NUM_THREADS
   * generic variable that specifies the default maximum number of
     threads when parallelizing
   * affects shared public libraries with internal parallelization
     such as lapack/blas
   * also used with QSpace via my wb::sparse() class to deal
     with generic sparse tensors (e.g. all Clebsch-Gordan data)

  MKL_NUM_THREADS
   * affecs the shared intel MKL libraries (similar to lapack/blas)
     note that these are delivered with the matlab API, and hence
     are the ones that are included in mex-files.

  QSP_NUM_THREADS

   * introduces a higher level parallelization within the QSpace
     tensor library. For example, consider the following piece
     of matlab-like pseudo-code

     | ee=cell(1,numel(A.data));
     |
     | parallel_for i=1:numel(A.data) <-- distribute this to QSP_NUM_THREADS threads
     |    ee{i}=parallel_eig(A.data{i}); <-- calculate using MKL_NUM_THREADS threads
     | end

     > each of the QSP_NUM_THREADS threads that gets assigned a
       unique value i in the loop therefore has a task [here eig()]
       which itself parallelizes, say using MKL_NUM_THREADS threads.
       This is referred to as `nested parallelization'.
  
     > hence the actual maximal number of threads is given by
       QSP_NUM_THREADS * max(OMP_NUM_THREADS, MKL_NUM_THREADS) !!

     > for cluster jobs: make sure that the number or requested
       cores matches, i.e. never exceeds (NTASK * QSP_NUM_THREADS).

   * the QSpace parallelization based on QSP_NUM_THREADS affects both,
     A.data{} as well as A.info.cgr(), e.g. parallelizes when working
     through reduced sets of reduced matrix elements as well as their
     corresponding generalized Clebsch Gordan coefficient data.

   * so far, QSP_NUM_THREADS mostly affects the contractQS() routine
     where by matching the symmtery sectors of two higher-rank tensors,
     say A and B, can quickly result in a huge list to work through.

   * a change in QSP_NUM_THREADS will be seen by mex routines
     if they are called after the changed (i.e. a 'clear functions'.
     may be required within matlab).

   * threads get locked (i.e. need to wait) when other threads
     e.g. already happen to work on the same Clebsch-Gordan tensor;
     this is based on standard OpenMP locking. Therefore large
     QSP_NUM_THREADS will results in overhead with less overall
     speedup.

   * I have tested QSP_NUM_THREADS with (NRG + fdmNRG) both, with
     abelian and non-abelian symmetries. If you see strange kinds of
     crashes try turning off QSP_NUM_THREADS (the reason being that
     different threads may get into each others way, despite that
     I tried to eliminate as much of possible of this scenario).

  The default value for the above variables is 1, i.e. do not parallelize
  (if you experience something different for mex-compiled routines,
  please let me know).

  Suggested usage:

   * set OMP_NUM_THREADS and MKL_NUM_THREADS the same, e.g.
     with NTASK the number of threads assigned to a given task,

        export OMP_NUM_THREADS $NTASK
        export MKL_NUM_THREADS $NTASK

        and within matlab: >> num_threads(NTASK);  % see MPS_PACK/lib

   * if one expects large A.data{} sets [which carry (reduced)
     matrix elements as in typical NRG calculations], you may simply
     ignore QSP_NUM_THREADS (default value for QSP_NUM_THREADS is 1).

   * if you are dealing with larger-rank tensors such as in PEPS
     that can have QSpace with lists of lengths hundreds to thousands,
     you man want to turn on QSP_NUM_THREADS while reducing NTASK
     accordingly.

# -------------------------------------------------------------------- #
# fixed memory leak in QSpace (data sector)

  thanks to Benedikt for pointing this out (e.g. within svdQS)
  Wb,Jun27,17

# -------------------------------------------------------------------- #
# altered behavior of Nkeep/stol in svdQS and orthoQS
# in response to comments by Wei Li // Wb,Sep17,17

  Nkeep is interpreted as Nkeep_max
   - if Nkeep is specified, but stol is not,
     the default value of stol (1E-8) changes to 0
   - added option 'Nkmin' to specify minimum number
     of states/multiplets to keep

# -------------------------------------------------------------------- #

