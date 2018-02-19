#!/usr/bin/perl
# Usage: eval "$(~/.matlab_setup.pl [opts])"
# 
#    this script generates the relevant bash settings in text form
#    for running matlab.
# 
# Options
# 
#   -v   verbose flag (verbose in that eval(...) prints extra information)
#   -q   quiet flag (no extra information written to STDERR)
#   -f   enforces (re)setup of variables
#   -F   enforces buildup of full LD_LIBRARY_PATH
# 
# Wb,Feb08,11

  sub getenv {
     if (@_<1) { die "\n  ERR getenv() invalid usage\n  ERR"; }
     my $v=shift; if (exists $ENV{$v}) { return $ENV{$v}; }
     if (@_) {
        if (@_>1) { die "\n  ERR getenv() invalid usage\n  ERR"; }
        return shift;
     }
     return '';
  }

# $,=' ';
  $\="\n  "; print;

  my $MLR='R2016a'; # current matlab version

  my $NR='NR_CALL_MATLAB_RC';
  my $nrcall=getenv($NR,0); ++$nrcall;

  my $fflag=0; my $vflag=0; my $qflag=0; # my $tflag=0;
  my $ldflag=0; my $useproj=0;

  $V='MVERSION'; if (exists $ENV{$V}) { $MLR=$ENV{$V}; $fflag=1; }
  
  while (@ARGV) { $_=$ARGV[0];
     if (/^-[h\?]$/) { die usage($0); }
     elsif ($_ eq '-v') { $vflag=1; shift; }
     elsif ($_ eq '-q') { $qflag=1; shift; }
   # elsif ($_ eq '-t') { $tflag=1; shift; }
     elsif ($_ eq '-f') { $fflag=1; shift; }
     elsif ($_ eq '-F') { $fflag=2; shift; }
     elsif ($_ eq '-ld') { $ldflag=1; shift; }
     elsif ($_ eq '-proj') { $useproj=1; shift; }
     elsif (/^-(\d+\w)/) { $MLR=$1; $fflag=1; shift; }
     else { last; }
  }
  if (@ARGV) { 
     die usage($0)."\n  Got extra arguments: ".join2args(@ARGV)."\n";
  }

  if ($useproj) {
     print("\n. ~/.switch2proj -ml\n");
   # WRN calls rebash => .matlab_setup
  }

  if ($MLR=~/^-*(\d+)(\w)/) { $MLR=sprintf('R2%03d%s',$1,$2); }
  my $ML=''; $_=$MLR; s/[A-Za-z]//;

  if ($_<2013) {
     $ML='/usr/local/dist/DIR/matlab-'.$MLR;
  }
  else {
   # $ML='/software/opt/precise/x86_64/matlab/2013a';
     $MLR=~s/R(\d)/matlab\/$1/; # remove leading R

   # my $cmd="bash -c 'module show matlab/$_'";
     my $cmd="/software/opt/noarch/Modules/1.147/modulecmd.tcl bash show $MLR";

     my @ll=`$cmd 2>&1`;
     if ($? || grep(/ERROR/,@ll)) {
        $_=$cmd; s/.*\'(.*)\'.*/$1/g;
        print "\n  WBERR=matlab_setup.pl:module\n";
        print STDERR "\n  $_ ($?) returned ...\n";
        die "  ".join("  ",@ll)."\n";
     }
     foreach (@ll) {
        if (/PATH\s*(.*)$/) {
           $ML=$1; $ML=~s/\/bin$//;
        }
     }
  }

  if ($ML eq '') {
     $_=$cmd; s/.*\'(.*)\'.*/$1/g;
     print STDERR "\n",@ll,"\n";
     die "\n  failed to determine MatLab path ($_)\n";
  }

  if ($vflag>1) { print STDERR $ML; exit 0; }

  print "export $NR=$nrcall\n";
  print "export MATLAB_VERSION=$MLR";  # to be used with 'module load ...'
  print "export ML_GCC_VERSION=gcc/4.7.4\n";  # to be used with 'module load ...'

  if ($ldflag && $fflag) {
     print "module load \$MATLAB_VERSION";
     print "module load \$ML_GCC_VERSION\n";
  }

# MYMATLAB is my Matlab home directory (used by Matlab/startup.m)
# MEXROOT is my MEX root directory (used by MEX/make)

  print 'export MYMATLAB=$HOME/Matlab';
  print 'export MEXROOT=$MYMATLAB/MEX'."\n";

# echo -e "\n  ${LD_LIBRARY_PATH//:/\\n  }" # ldpath

  if (exists $ENV{'PROJECT_ENV'}) {
         print 'export MLBLOG=${MYMATLAB/$HOME/$HOME_OLD}/mbatch.log'; }
  else { print 'export MLBLOG=$MYMATLAB/mbatch.log'; }

# NB! use LD_LIBRARY_PATH path only if really necessary
# eg. call to libfftw3 crashed since LD_LIBRARY_PATH was set
# incidentially for MatLab paths!! tags: ldflag LD_LIBRARY_PATH
# Wb,Jun05,07

  my $march='(unknown)'; my $mcc; my $t;
  foreach (`uname -m`) { chomp; $t=$_;
        if (/i686/  ) { $march='glnx86';  $mcc='binglx'; }
     elsif (/x86_64/) { $march='glnxa64'; $mcc='bina64'; }
     else { die "\n  ERR invalid uname $_ ??\n ERR"; }
  }

  print 'export MCC_BIN=$MCC/'.$mcc."\n";

# -------------------------------------------------------------------- #
# Ralph, Jun25,07
#    /usr/local/dist/DIR/...    should be used instead of
#    /amnt/dist/DIR/... 
# the later is an automatic mount point which may be system
# dependent (even though the same for the local system)
# -------------------------------------------------------------------- #

  print 'export MATLAB_ROOT='.$ML,"\n";

# -------------------------------------------------------------------- #
# needed by MatLab runtime library

  print 'export ARCH='.$march."\n";

# NB! mex requires environmental variable ARCH
# to avoid warnings of type
# I18N Runtime Warning: Missing ICU data file 
#    detected while processing $(MATLAB_ROOT)/bin/$(ARCH).
#    Hint: Check for a misconfigured environment or installation.
# -------------------------------------------------------------------- #

  my $LD=$ENV{'LD_LIBRARY_PATH'};
  my $ld='(null)'; if (($ldflag || $vflag) && $fflag<2) {
     my @ll=grep(/LD_LIBRARY_PATH/,
       `unset LD_LIBRARY_PATH ; $ML/bin/matlab -nodesktop -nodisplay -n`);
     if ($? || @ll!=1) { die
     "\n  ERR failed to determine matlab default LD_LIBRARY_PATH\n "; }
     $ld=shift(@ll);
  }

if ($nrcall==1 || $fflag) {

  print 'export SGE_HIST=~/.sge/sge_hist';

# export MCRROOT=$MATLAB_ROOT
# export MCROOT=$MATLAB_ROOT
# export MROOT=$MATLAB_ROOT

  if ($ldflag) {
     if (!($ld=~/extern\/lib/) || !($ld=~/runtime\/$ARCH/) || !($ld=~/sys\/os/)) {
     if ($fflag || !($LD=~/$march/) && !($LD=~/matlab/)) {
        print 'export LD_LIBRARY_PATH=".:'.
          '$MATLAB_ROOT/runtime/$ARCH:'.
          '$MATLAB_ROOT/bin/$ARCH:$MATLAB_ROOT/sys/os/$ARCH:'.
          '$MATLAB_ROOT/extern/lib/$ARCH"'."\n";
     }}
     else { ++$ldflag;
      # print 'export LD_LIBRARY_PATH=.:/usr/local/dist/lib:/usr/local/common/lib',"\n";
        print 'export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH',"\n";
        if (!$qflag) { print STDERR
        "\r".'# using default LD_LIBRARY_PATH set my matlab (runtime)'; }
     }
  }

# -------------------------------------------------------------------- #
# setup for mcc (formerly ~/.mcc_setup)

# removes a few warning regarding java
# present even though -R -nojvm was specified below

# if [ ]; then
# LD_LIBRARY_PATH actually appears to be adjusted correctly automatically

  if ($ldflag) {
     if (!($ld=~/extern\/lib/) || !($ld=~/runtime\/$ARCH/) || !($ld=~/sys\/os/)) {
     if ($fflag || !($LD=~/java/)) { # see mcc docu (6.4, p.102)

      # see mcc output -> start script .sh
      # MatLab 2008a switched from jre1.5.0 to jre 1.6.0
        print 'JAV="$MATLAB_ROOT/sys/java/jre/$ARCH"';

        my $JAV=$ML.'/sys/java/jre/'.$march;
        my $f=$JAV.'/jre.cfg'; # contains version number
        if (-f $f) {
         # for older matlab version (no longer there for Matlab > 2009)
         # Wb,Oct06,09
           my @ver=`cat $f`; $_=shift(@ver); chomp;
           $JAV.="/jre$_/lib";
        }
        else { $JAV.="/jre/lib"; }

        $_=$t;
           if (/i686/  ) { $JAV.='/i386'; }
        elsif (/x86_64/) { $JAV.='/amd64'; }
        else { die; }

        print 'export XAPPLRESDIR="$MATLAB_ROOT/X11/app-defaults"',"\n";
        print 'export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:'.
        $JAV.'/native_threads:'.$JAV.'/server/:'.$JAV.'/client/:'.$JAV,"\"\n";
     }}
     else { ++$ldflag; if (!$qflag) {
        print STDERR "\r".'# using default LD_LIBRARY_PATH set my matlab (java)';
        if ($ldflag==3) { $_=$ld; s/:/\n# | /g; s/^.*=\s+//;
           print STDERR "\r# | ",$_;
        }
     }}
  }

# OMP_NUM_THREADS
  $V='OMP_NUM_THREADS';
  if (!exists $ENV{$V}) {
     my $n=grep(/^processor/,`cat /proc/cpuinfo`);
     print 'export OMP_NUM_THREADS='.($n-4)."\n";
  }

} # of ($nrcall==1 || $fflag)

  if ($vflag) { print 'cat << EOT__

# -------------------------------------------------------------------------- #
# MatLab Environment (~/.matlab_setup.pl, np=$OMP_NUM_THREADS) on $(date)

  machine: $(hostname)/$(pwd)
  MATLAB : $MATLAB_ROOT
  ARCH   : $ARCH

# XAPPLRESDIR: $XAPPLRESDIR
# LD_LIBRARY_PATH:
';
   if ($ldflag==3 && $qflag) {
      $_=$ld; s/:/\n  /g; s/^.*=\s+//;
      print "\r# using matlab default LD_LIBRARY_PATH\n  ",$_;
   }

   print 
   '$(echo $LD_LIBRARY_PATH | sed \'s/:/\n  /g\' | \
    sed "s/${MATLAB_ROOT//\//\\/}/MATLAB\//g")'."\n\n".
'# -------------------------------------------------------------------------- #'.
   "\nEOT__\n";

 # print "\n  TST ldflag=$ldflag, qflag=$qflag, vflag=$vflag, ld=$ld\n";
}

# echo -e "\n  ${LD_LIBRARY_PATH//:/\\n  }" # ldpath

