# KPIC: Build scripts for calibration
include ../../Mk.instrument

override SYSNAM = kss/$(INSTRUMENT)/operations/Calibration
override VERNUM = 1.0

override ENABLE_PYTHON2 = False

ifdef TESTBED
	RELLIB = spec_scans_lab.py
else
	RELLIB = spec_scans.py spec_scans_combined.py kpic_background.py kpic_extraction.py kpic_trace.py VFN_Toolkit_summit.py extract_spec_flux.py
endif
LIBSUB = python
RELBIN = TT_Scan

################################################################################
# KROOT boilerplate:
# Include general make rules, using default values for the key environment
# variables if they are not already set.

ifndef KROOT
	KROOT = /kroot
endif

ifndef RELNAM
	RELNAM = default
endif

ifndef RELDIR
	RELDIR = $(KROOT)/rel/$(RELNAM)
endif

include $(RELDIR)/etc/defs.mk
################################################################################

ifneq "$(PYTHON3)" ""
	ENABLE_PYTHON3 = True
endif

################################################################################
include $(RELDIR)/etc/rules.mk
################################################################################
