ifndef MAKEFILE_VARS
MAKEFILE_VARS := 1
#...OS name (no underscores)
L_EMPTY :=
SPACE = $(L_EMPTY) $(L_EMPTY)
# uname vars
UNAME_N         := $(shell uname -n)
UNAME_S         := $(shell uname -s)
UNAME_M         := $(shell uname -m)
# machine vars
L_OS_BASE := $(subst -,,$(subst $(SPACE),,$(subst _,,$(UNAME_S)$(UNAME_M))))
# strip off non-word and trailing digits and upper case for default names
L_NAME_SHORT := $(shell perl -e '($$a="$(UNAME_N)")=~s/\W.*$$//;$$a=~s/\d+$$//;$$a=~tr/a-z/A-Z/;print $$a;')
# L_CLASS                    = what class of machine you are on
# L_MACHINE                  = longer name of the machine you are on
# L_OS                       = short  name of the machine you are on
# L_END                      = FRONT or BACK
# L_BATCH_ARGS_D         = default flags to use in places (tests, myllogini, ... )
#                              (like account, queue, ...)
# rules for these names:
# L_OS: 3 characters in length max
L_MACHINE = M_$(L_NAME_SHORT)
L_CLASS   = C_$(L_NAME_SHORT)
L_OS      = $(L_OS_BASE)
L_GROUP   = flag
L_GROUP_SOURCE = shavano
L_GROUP_OTHER = $(L_GROUP_SOURCE)
STAT = stat
#...Linux
ifeq "$(L_OS_BASE)" "Linuxx8664"
  #...Cielo
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^(ci-)/' $(UNAME_N))" ""
    L_CLASS   = CIELO
    L_MACHINE = CIELO
    L_OS      = CI
    L_PPN     = 16
  endif
  #...Cielito
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^(ct-)/' $(UNAME_N))" ""
    L_CLASS   = CIELO
    L_MACHINE = CIELITO
    L_OS      = CI
    L_PPN     = 16
  endif
  #...Smog
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^(sm-)/' $(UNAME_N))" ""
    L_CLASS   = CIELO
    L_MACHINE = SMOG
    L_OS      = CI
    L_PPN     = 16
  endif
  #...ROADRUNNER3
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^rr/' $(UNAME_N))" ""
    L_CLASS   = RR3
    L_MACHINE = ROADRUNNER3
    L_OS      = RR3
  endif
  #...luna
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^lu/' $(UNAME_N))" ""
    L_CLASS   = TLCC
    L_MACHINE = LUNA
    L_OS      = LU
  endif
  #...Typhoon
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^ty/' $(UNAME_N))" ""
    L_CLASS   = TLCC
    L_MACHINE = TYPHOON
    L_OS      = TY
  endif
  #...Moonlight
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^ml/' $(UNAME_N))" ""
    L_CLASS   = TLCC
    L_MACHINE = MOONLIGHT
    L_OS      = ML
  endif
  #...mapache
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^mp/' $(UNAME_N))" ""
    L_CLASS   = TLCC
    L_MACHINE = MAPACHE
    L_OS      = MP
  endif
  #...mustang
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^mu/' $(UNAME_N))" ""
    L_CLASS   = TLCC
    L_MACHINE = MUSTANG
    L_OS      = MU
  endif
  #...Pinto
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^pi/' $(UNAME_N))" ""
    L_CLASS   = TLCC
    L_MACHINE = PINTO
    L_OS      = PI
  endif
  #...DARWIN_GPU: front
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /darwin\.lanl\.gov/' $(UNAME_N))" ""
    L_CLASS   = DARWIN_GPU
    L_MACHINE = DARWIN_GPU
    L_OS      = DGPU
  endif
  #...DARWIN_GPU
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^cn\d+/' $(UNAME_N))" ""
    L_CLASS   = DARWIN_GPU
    L_MACHINE = DARWIN_GPU
    L_OS      = DGPU
  endif
  #...Wolf
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^wf/' $(UNAME_N))" ""
    L_CLASS   = TLCC
    L_MACHINE = WOLF
    L_OS      = WF
  endif
endif
#...Mac
ifeq "$(UNAME_S)" "Darwin"
    L_CLASS   = DARWIN
    # if need be, hardwire this to 1
    #L_PPN     = 1
endif
#...DAWN: front end (dawn\d+)
ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^(dawn\d+)/' $(UNAME_N))" ""
  L_CLASS   = DAWN
  L_MACHINE = DAWN
  L_OS      = DAWN
  L_PPN     = 1
endif
#...SEQUOIA: front end (seq\d+)
ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^(seqlac\d+)/' $(UNAME_N))" ""
  L_CLASS   = SEQUOIA
  L_MACHINE = SEQUOIA
  L_OS      = SQ
endif
#...Vulcan: front end (vul\d+)
ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^(vulcanlac\d+)/' $(UNAME_N))" ""
  L_CLASS   = SEQUOIA
  L_MACHINE = VULCAN
  L_OS      = SQ
endif

#...L_PPN
#...L_PPN from sinfo command
ifeq "$(L_PPN)" ""
  TMP_CMD := $(shell which sinfo 2> /dev/null )
  ifneq "$(TMP_CMD)" ""
    # L_PPN = SOCKETS * CORES = CPUS
    #  Correct for all machines except sequoia (and cielo, which does not have sinfo)
    #  If CPUS has non-digit (seq returns 8K), just use CORES
    # %c: CPUS    (SOCKETS * CORES)
    # %X: SOCKETS 
    # %Y: CORES   (correct for sequoia)
    L_PPN := $(shell $(TMP_CMD) -o "%c %X %Y" 2> /dev/null | perl -ne 'if(/CPUS/){$$found="yes"} elsif( defined($$found) ){ @fields = split(/\s+/); if( $$fields[0] =~ /\D/ ) { print $$fields[2] } else{ print $$fields[0] }; exit;}')
  endif
endif
#...L_PPN
#...L_PPN from moab config file
ifeq "$(L_PPN)" ""
  MOAB_FILE = $(shell ls /opt/MOAB/moab.cfg 2> /dev/null )
  ifneq "$(MOAB_FILE)" ""
    L_PPN := $(shell perl -ne 'if( /^[^\#].*MAXPE=(\d+)/ ){print $$1;exit;}' $(MOAB_FILE))
  endif
endif
#...try also run lcstat
ifeq "$(L_PPN)" ""
  TMP_CMD := $(shell which lcstat 2> /dev/null )
  ifneq "$(TMP_CMD)" ""
    L_PPN := $(shell $(TMP_CMD) 2> /dev/null | perl -ne 'if(/\d+:(\d+)/){print $$1; exit;}')
  endif
endif
# L_PPN
ifeq "$(L_PPN)" ""
  L_PPN := $(shell cat /proc/cpuinfo 2> /dev/null | egrep '^processor' | wc -l | perl -pe 's/\s+//g')
  ifeq "$(L_PPN)" "0"
     L_PPN =
  endif
endif
# default L_PPN
ifeq "$(L_PPN)" ""
  L_PPN = 16
endif

# The following are created for determining if libraries and executables
# are compatible to load/run on different machines.
# Used by prebuilt libraries.
#
# L_TYPE      : T+<shortened name for compatible execs>
# L_TYPE_FULL : Full name used for generating L_TYPE
#                     L_CLASS + chip info    or
#                     L_MACHINE
# To allow for mappings, create the following files:
#    L_TYPE_FULL_$(L_TYPE)
#  Containing:
#     $(L_TYPE)=$(L_TYPE_FULL)
#  and soft link L_OS_$(L_OS) -> L_TYPE_FULL_$(L_TYPE)
L_TYPE_FULL:=
L_TYPE:=
# L_TYPE_FULL: from /proc/cpuinfo
ifeq "$(L_TYPE_FULL)" ""
  L_TYPE_FULL := $(shell cat /proc/cpuinfo 2> /dev/null | perl -ne 'if(/^(model name|cpu)\s*:\s*(\S.*)/){($$chip=$$2)=~s/\s+//g; print "cpuinfo=$$chip"; exit}')
endif
# L_TYPE_FULL: if still not defined, use L_MACHINE
ifeq "$(L_TYPE_FULL)" ""
  L_TYPE_FULL := L_MACHINE=$(L_MACHINE)
endif
# L_TYPE_FULL: prepend with L_CLASS
L_TYPE_FULL := L_CLASS=$(L_CLASS),$(L_TYPE_FULL)
ifeq "$(L_TYPE)" ""
  # full | cksum (prints checksum and number of bytes) |
  #        perl mod (digits) (ignores second part and digits)
  L_TYPE := $(shell echo '$(L_TYPE_FULL)' | cksum 2> /dev/null | perl -ne '$$digits=4;printf "%0$${digits}x",$$_%(16**$$digits)' )
endif
ifeq "$(L_TYPE)" ""
  L_TYPE = L_OS=$(L_OS)
endif
# ensure starts with letter (so "T")
L_TYPE := T$(L_TYPE)

#... L_END
# back if SLURM_JOB_ID set
ifeq "$(L_END)" ""
  ifneq "$(SLURM_JOB_ID)" ""
    L_END = BACK
  endif
endif
# back if PBS_JOBID set
ifeq "$(L_END)" ""
  ifneq "$(PBS_JOBID)" ""
    L_END = BACK
  endif
endif
# if above environment vars not set, but msub still exists, then front end
ifeq "$(L_END)" ""
  TMP_CMD := $(shell which msub 2> /dev/null )
  ifneq "$(TMP_CMD)" ""
    L_END = FRONT
  endif
endif
# front if machine named -fe
ifeq "$(L_END)" ""
  ifneq "$(shell perl -e '$$_=$$ARGV[0]; print if /^(\w+-fe)/' $(UNAME_N))" ""
    L_END = FRONT
  endif
endif
# default is to be on the "BACK" end - where you can run interactively
ifeq "$(L_END)" ""
  L_END = BACK
endif
# L_INSTALL_DIR can be gotten from the location of a special location script
ifeq "$(L_INSTALL_DIR)" ""
  L_INSTALL_DIR := $(dir $(shell eap_tools_dir 2> /dev/null))
endif
ifeq "$(L_INSTALL_DIR)" ""
  L_INSTALL_DIR := /ERROR/L_INSTALL_DIR/not/defined/a/b/c
endif
L_INSTALL_DIR := $(patsubst %/,%,$(L_INSTALL_DIR))
L_INSTALL_DIR_UP := $(dir $(L_INSTALL_DIR))
#...L_BATCH_ARGS_D: account for quick turnaround
# mdiag can take forever, so if already set, do not redo
ifndef L_BATCH_ARGS_D
  TMP_CMD := $(shell which mdiag 2> /dev/null )
  ifneq "$(TMP_CMD)" ""
    L_BATCH_ARGS_D_try := $(shell $(TMP_CMD) -u $(LOGNAME) 2> /dev/null | perl -ne 'if( /\b(access|drlanl|lanl-exec|lanlexec)\b/ ){print "-A $$1"; exit;}')
    L_BATCH_ARGS_DI_try := $(L_BATCH_ARGS_D_try)
  endif
  #...L_BATCH_ARGS_D: queue (seq has it)
  #...see if some automated way to get this in general
  ifeq "$(L_OS)" "SQ"
    # msub
    L_BATCH_ARGS_D_try  += -q psmall
    # myllogini
    L_BATCH_ARGS_DI_try += -p psmall
  endif
  ifeq "$(L_OS)" "DAWN"
    L_BATCH_ARGS_D_try += -q pdebug
  endif
  L_BATCH_ARGS_D = $(L_BATCH_ARGS_D_try)

  # set to a space so that we do not go through this again
  # for gmake, "" is the same as ifndef
  ifeq "$(L_BATCH_ARGS_D)" ""
    L_BATCH_ARGS_D = $(SPACE)
  endif

  # myllogini
  ifndef L_BATCH_ARGS_DI
    L_BATCH_ARGS_DI = $(L_BATCH_ARGS_DI_try)
  endif
endif

# find flavors of various execs (specifically mac flavors)
# do not know if this is the best place for it
#    $(shell date $(L_DATE_FLAG_SPECIFY) "yesterday" +%u.%a)
# 
ifeq "$(L_EXEC_DATE)" ""
  TMP_CMD := $(shell date --v 2> /dev/null)
  ifneq "$(TMP_CMD)" ""
    L_EXEC_DATE = gnu
  else
    L_EXEC_DATE = mac
  endif
endif
ifeq "$(L_EXEC_FIND)" ""
  TMP_CMD := $(shell find -version 2> /dev/null)
  ifneq "$(TMP_CMD)" ""
    L_EXEC_FIND = gnu
  else
    L_EXEC_FIND = mac
  endif
endif
ifeq "$(L_EXEC_TAIL)" ""
  TMP_CMD := $(shell tail --version 2> /dev/null)
  ifneq "$(TMP_CMD)" ""
    L_EXEC_TAIL = gnu
  else
    L_EXEC_TAIL = mac
  endif
endif

# flag for date
ifeq "$(L_EXEC_DATE)" "gnu"
  L_DATE_FLAG_SPECIFY = "-d"
else
  L_DATE_FLAG_SPECIFY = "-j -f"
endif

# eap specific 
L_LS := $(shell test -e Makefile_vars.mk && echo exists)
ifneq "$(L_LS)" ""
  L_ROOT_SOURCE_DIR := $(dir $(shell readlink Makefile_vars.mk))
  # will be just "." if not softlinked
  ifeq "$(L_ROOT_SOURCE_DIR)" ""
    L_ROOT_SOURCE_DIR = "."
  # otherwise, strip off 2 dirs from Tools.rh/Environment
  else
    L_ROOT_SOURCE_DIR := $(patsubst %/,%,$(dir $(patsubst %/,%,$(dir $(patsubst %/,%,${L_ROOT_SOURCE_DIR})))))
  endif
endif
ifeq "$(L_LS)" ""
  L_ROOT_SOURCE_DIR = /unknown/need/to/fix/Makefile_vars.mk
endif

# ============
# version info
# ============
L_V_HEADURL = $$HeadURL: file:///usr/projects/eap/svnroot/eap.rh/eap/tags/Test/NIGHTLY.2.Tue/eap.rh/Tools.rh/Environment/Makefile_vars.mk $$
L_V_REVISION = $$Revision: 1.12 $$
L_V_REVNUM = $(firstword $(patsubst $$Revision: 1.12 $(L_V_REVISION)))
L_V_MODVER_DEF = 1405
L_V_MODEXT_DEF = 12
# convuluted way to pull out <MODVER>.<MODEXT> from L_V_HEADURL
#   L_V_HEADURL = <stuff>$(_release_path)<MODVER>.<MODEXT>/<stuff>
_release_path = /tags/Releases/v_
L_V_MODEXT_VER := $(firstword $(subst /, ,$(patsubst $(_release_path)%,%,$(filter $(_release_path)%,$(subst $(_release_path), $(_release_path),$(L_V_HEADURL))))))
ifneq "$(L_V_MODEXT_VER)" ""
  L_V_MODVER = $(basename $(L_V_MODEXT_VER))
  L_V_MODEXT = $(patsubst .%,%,$(suffix $(L_V_MODEXT_VER)))
  L_V_REL    = release
else
  L_V_MODVER = $(L_V_MODVER_DEF)
  L_V_MODEXT = $(L_V_MODEXT_DEF)
  L_V_REL    = none
endif
DATE_YMD := $(shell date +%Y%m%d)
ifeq "$(L_V_REL)" "release"
  L_V_REL_DIR = $(L_INSTALL_DIR)/releases/$(L_V_MODVER).$(L_V_MODEXT)
else
  L_V_REL_DIR = $(L_INSTALL_DIR)/releases/special/$(DATE_YMD)
endif

# printing rules
# example:
#   gmake -f Makefile_vars.mk L_END
PACK_HELP_MSG = \
	@ echo "" ; \
	  echo "NOTE: To build a local developer executable:" ; \
	  echo "          ./pack.pl" ; \
	  echo "          cd <build directory>" ; \
	  echo "          $(MAKE) -j" ;

# first target is all to make "gmake" default to "gmake all"
all:

L_VARS = \
        L_BATCH_ARGS_D \
        L_BATCH_ARGS_DI \
        L_EXEC_DATE \
        L_EXEC_FIND \
        L_EXEC_TAIL \
	L_TYPE \
	L_TYPE_FULL \
	L_CLASS \
	L_INSTALL_DIR \
	L_INSTALL_DIR_UP \
        L_V_MODEXT \
        L_V_MODVER \
        L_V_REL \
        L_V_REVNUM \
	L_END \
	L_GROUP \
	L_GROUP_SOURCE \
	L_GROUP_OTHER \
	L_MACHINE \
	L_OS \
	L_OS_BASE \
	L_PPN \
	L_ROOT_SOURCE_DIR
.PHONY: $(L_VARS)

# parsed by various scripts
# or used in .cshrc in eval
# if DO_ENV is set, then
#  print environment variable setenv
#  All vars of the form L_<VAR>
#   easiest to remove all L then put all back (2 patsubst)
$(L_VARS):
ifeq "$(DO_ENV)" ""
  ifeq "$(val_only)" ""
	@ echo "$@=$($@)"
  else
	@ echo "$($@)"
  endif
else
	@ echo "setenv $(patsubst L_%,L_%,$(patsubst L_%,L_%,$@)) '$($@)' ; "
endif

# for performance, explicitly echo out
print_makefile_vars pmv:
	@ echo "L_BATCH_ARGS_D=$(L_BATCH_ARGS_D)"
	@ echo "L_BATCH_ARGS_DI=$(L_BATCH_ARGS_DI)"
	@ echo 'L_TYPE=$(L_TYPE)'
	@ echo 'L_TYPE_FULL=$(L_TYPE_FULL)'
	@ echo L_CLASS=$(L_CLASS)
	@ echo L_EXEC_DATE=$(L_EXEC_DATE)
	@ echo L_EXEC_FIND=$(L_EXEC_FIND)
	@ echo L_EXEC_TAIL=$(L_EXEC_TAIL)
	@ echo L_INSTALL_DIR=$(L_INSTALL_DIR)
	@ echo L_INSTALL_DIR_UP=$(L_INSTALL_DIR_UP)
	@ echo L_END=$(L_END)
	@ echo L_GROUP=$(L_GROUP)
	@ echo L_GROUP_SOURCE=$(L_GROUP_SOURCE)
	@ echo L_GROUP_OTHER=$(L_GROUP_OTHER)
	@ echo L_MACHINE=$(L_MACHINE)
	@ echo L_OS=$(L_OS)
	@ echo L_OS_BASE=$(L_OS_BASE)
	@ echo L_PPN=$(L_PPN)
	@ echo L_ROOT_SOURCE_DIR=$(L_ROOT_SOURCE_DIR)
endif # MAKEFILE_VARS
