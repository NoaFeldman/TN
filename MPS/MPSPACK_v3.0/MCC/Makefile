# -------------------------------------------------------------------- #
# Makefile for MatLabs standalone compiler mcc
# See also i make (or doku make.pdf)
#
# NB! MatLab 2006b appears to be the latest version of MatLab that
# does run on the cluster through mcc compilation or as batch
# all later MatLab version only start as zombie processes !?@*%@!
# Wb,Nov16,09
#
# Wb,May23,07
# -------------------------------------------------------------------- #

  SHELL = /bin/bash

  ifndef MEXROOT
      $(error missing environmental variable MEXROOT)
  endif

  ifndef MYMATLAB
      $(error missing environmental variable MYMATLAB)
  endif

  MCC_DIR = $(MEXROOT)/../MCC

# MAKEFLAGS += -B    # compile always
# export MVERSION=6b # use matlab version 2006b - tags: CTF_FLAG

  MCC := bin/mcc -C -m -w enable -R -nodisplay
  # MCC options:
  #
  #  -C   create CTF archive (crucial for MatLab > 2006b to work!
  #       tags: CTF_FLAG ; Wb,Nov18,09
  #
  #   -g  include debug options
  #   -v  verbose
  #
  #   -R -nojvm
  #   -R -nodisplay
  #
  # See also mcc -help.

  #-I $(MYMATLAB)/       \

  MLI := \
   -I $(MYMATLAB)/libx   \
   -I $(MYMATLAB)/lib    \
   -I $(MYMATLAB)/QCI    \
   -I $(MYMATLAB)/MPS    \
   -I $(MYMATLAB)/NRG    \
   -I $(MYMATLAB)/MEX    \
   -I $(MYMATLAB)/Class  \
   -I $(MYMATLAB)/Phys   \
   -I $(MYMATLAB)/Verstraete/PEPS \
   -I $(MEXROOT)/bin     \
   -I $(MEXROOT)/lib     \
   -I $(MEXROOT)/GetPsi  \

# NB! cto / cd only applies to the shell for given command line i.e. no
# need to return (where cd - would result in another pwd output line !?)

# @( cd binglx; pwd ; echo $(shell pwd) ; pwd )
# shell pwd shows the *old* path !?!

# See i make (or doku make.pdf)

  MINC = run save2 global basename inputname 

# -------------------------------------------------------------------- #
# implicit rule (assumes compilation of %.m file) // default

%:
	@( cd $(MCC_DIR) ; \
       if [ $(ARCH) = glnxa64 ]; then cd bina64; else cd binglx; fi; \
	   module matlab -v -ld -f ; \
	   printf '  cd '; pwd ;   \
	   printf "\n  $(MCC) $(@:.m=)\n  -a $(MINC)...\n\n"; \
	   $(WHICH_MATLAB)/$(MCC) $(@:.m=) -a $(MINC) $(MLI) \
	);

#	   printf "\n  cmd: $(WHICH_MATLAB)/$(MCC) $(@:.m=) $(MLI)\n\n" ; \

#	  -R -startmsg mcr_start -R -completemsg mcr_complete \
# with mcr_start.m and mcr_compete.m are functions defined
# on the matlab path ($ML) does not seem to have any effect

# -------------------------------------------------------------------- #

# will always be executed when listed as a dependence

safemsg:
	@( printf "\n" )

safemsg_old:
	@( printf "\n"; \
	   printf "  NB! make sure no MCC processes are left running!\n"; \
	   i=3; while [ 1 ]; do \
	   if [ $$i -le 0 ]; then break; fi; \
	   printf "\r  Starting compilation in $$i sec ...  "; \
	   i=$$[i-1]; sleep 1; done; printf "\n\n"; \
	);

hworld_old: safemsg
	@( cd $(MCC_DIR) ; \
	   if [ $(ARCH) = glnxa64 ]; then cd bina64; else cd binglx; fi; \
	   module matlab -v -ld -f ; \
	   printf '  cd '; pwd ;   \
	   printf "\n  $(MCC) $@ -o $@ -a run ...\n\n"; \
	    $(WHICH_MATLAB)/$(MCC) $@ -o $@ -a run $(MLI) \
	);


rnrg_SUN_kotliar: safemsg
	@( cd $(MCC_DIR) ; if [ $(ARCH) = glnxa64 ]; then \
	       cd bina64; else cd binglx; \
	   fi; \
	   module matlab -v -ld -f ; \
	   printf '  cd '; pwd ;   \
	   printf "\n  $(MCC) $@ -o $@ -a run getGC_JJH setglobal ...\n\n"; \
	    $(WHICH_MATLAB)/$(MCC) $@ -o $@ -a run getGC_JJH setglobal $(MLI) \
	);
	@rsync2proj.pl --mcc

tst_Hamilton1D: safemsg
	@( cd $(MCC_DIR) ; if [ $(ARCH) = glnxa64 ]; then \
	       cd bina64; else cd binglx; \
	   fi; \
	   module matlab -v -ld -f ; \
	   printf '  cd '; pwd ;   \
	   printf "\n  $(MCC) $@ -o $@ -a run setglobal ...\n\n"; \
	    $(WHICH_MATLAB)/$(MCC) $@ -o $@ -a run setglobal $(MLI) \
	);
	@rsync2proj.pl --mcc

# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #

