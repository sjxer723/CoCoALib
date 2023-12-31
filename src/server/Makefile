# Makefile for CoCoALib/src/server directory

COCOA_ROOT=../..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=src/server/

DEPEND_FILE=Makefile_dependencies

SRCS=CoCoA5io.C \
     GlobalIO.C  SocketStream.C  \
     RegisterServerOps.C  RegisterServerOpsUser.C  RegisterServerOpsFrobby.C  \
     ServerOp.C  CoCoAServer.C

OBJS=$(SRCS:.C=.o)

EXECS=CoCoAServer

.PHONY: default
default: CoCoAServer


# This target assumes that the CoCoA library has already been compiled.
# Normally this target is used only as part of "make" from the
# CoCoALib root directory.
all: CoCoAServer
	touch all

CoCoAServer: $(OBJS) $(COCOA_LIB)
	$(COMPILE) $(OBJS) -o CoCoAServer $(LDLIBS)  # just LINKS the object files


.PHONY: benchmarks
benchmarks:
	@(cd benchmarks; $(MAKE) benchmarks)


.PHONY: clean clean-local
clean: clean-local
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -f core a.out all  ./*~  ./.*~  ./.\#*
	@/bin/rm -f  $(EXECS)  $(OBJS)  "$(DEPEND_FILE).old"  "$(DEPEND_FILE).new"
	@/bin/rm -rf  ./*.dSYM

.PHONY: veryclean veryclean-local
veryclean: clean-local
	@true >| "$(DEPEND_FILE)" # Empties $(DEPEND_FILE)
	@echo "Verycleaned CoCoALib/$(CWD)"


# This target checks whether we should really build DEPEND_FILE
.PHONY: dependencies
dependencies:
	@if [ \! -s "$(DEPEND_FILE)" -o \! "$(DEPEND_FILE)" -nt "$(COCOA_HDR)" ] ; \
	then \
	  $(MAKE) create-dependencies; \
	fi

# This is the target which really does rebuild the DEPEND_FILE
.PHONY: create-dependencies
create-dependencies:
	@/bin/rm -rf $(DEPEND_FILE).new
	@echo "Rebuilding dependencies in <CoCoALib>/$(CWD)"
	@echo "##################################################################"   >> $(DEPEND_FILE).new
	@echo "#       ---  CoCoAServer with CoCoALib-$(COCOALIB_VERSION)  ---"      >> $(DEPEND_FILE).new
	@echo "# Dependencies for directory $(CWD)"                                  >> $(DEPEND_FILE).new
	@echo "# DO NOT EDIT: this file is automatically generated; to regenerate"   >> $(DEPEND_FILE).new
	@echo "# run \`make dependencies' in the CoCoALib root directory."           >> $(DEPEND_FILE).new
	@echo "##################################################################"   >> $(DEPEND_FILE).new
	@echo                                                                        >> $(DEPEND_FILE).new
	@for srcfile in $(SRCS); \
	do \
	  objfile=`basename $$srcfile .C`.o; \
	  $(CXX) $(CXXFLAGS) $(CXXFLAGS_DEFINES) $(INCLUDE_PATHS)  -MM "$$srcfile" -MQ "$$objfile" >> $(DEPEND_FILE).new; \
	  echo >> $(DEPEND_FILE).new; \
	done
	@if [ -s "$(DEPEND_FILE)" ] ; \
	then \
	  /bin/mv -f "$(DEPEND_FILE)" "$(DEPEND_FILE).old" ; \
	fi
	@/bin/mv $(DEPEND_FILE).new  $(DEPEND_FILE)


$(OBJS): $(COCOA_ROOT)/configuration/autoconf.mk

$(DEPEND_FILE):
	@touch "$(DEPEND_FILE)"

include $(DEPEND_FILE)

#############################################################################

# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/server/Makefile,v 1.24 2020/10/27 09:44:35 abbott Exp $
# $Log: Makefile,v $
# Revision 1.24  2020/10/27 09:44:35  abbott
# Summary: Changed name of make var INCLUDE to INCLUDE_PATHS; removed very old var GMP_INCLUDE
#
# Revision 1.23  2020/06/17 18:43:44  abbott
# Summary: Removed cruft
#
# Revision 1.22  2020/03/06 20:35:15  abbott
# Summary: Improved behaviour of building Makefile_dependencies (even when that file does not exist)
#
# Revision 1.21  2020/01/27 19:59:08  abbott
# Summary: Removed some cruft
#
# Revision 1.20  2020/01/14 21:03:24  abbott
# Summary: Commented out dependency of Makefile_depend on configuration/version
#
# Revision 1.19  2019/12/18 08:16:49  abbott
# Summary: Revised configure; now must explicitly request to build Qt GUI; changed VERSION into COCOALIB_VERSION; conseq changes to Makefiles
#
# Revision 1.18  2019/03/20 16:20:55  bigatti
# -- updates for the new server
#
# Revision 1.17  2017/07/21 13:27:46  abbott
# Summary: Building of Makefile_dependencies should now be "clean" even when interrupted
#
# Revision 1.16  2017/02/14 17:06:29  abbott
# Summary: Updated clean/veryclean targets in all Makefiles
#
# Revision 1.15  2014/07/28 15:00:19  abbott
# Summary: Improved targets (very)clean-local; new scheme for building dependencies
# Author: JAA
#
# Revision 1.14  2014/06/28 15:04:20  abbott
# Summary: Revised in line with new configuration scripts; corrected generation of dependencies
# Author: JAA
#
# Revision 1.13  2014/06/17 10:20:05  abbott
# Summary: depend_file now depends on library.H (instead of version)
# Author: JAA
#
# Revision 1.12  2014/05/15 12:31:34  abbott
# Summary: Now using new files server/GlobalIO.HC (previously in CoCoA/io.H)
# Author: JAA
#
# Revision 1.11  2014/04/17 12:22:26  abbott
# Summary: Changed clean-local, veryclean-local; new approach for rebuilding dependencies
# Author: JAA
#
# Revision 1.10  2014/03/14 11:03:58  abbott
# Summary: clean and veryclean targets now delete .*~ files too
# Author: JAA
#
# Revision 1.9  2013/05/27 12:58:18  abbott
# Consequential change from moving all server-related code into src/server/
#
# Revision 1.8  2012/10/17 13:51:03  abbott
# Added helpful(!) comment.
#
# Revision 1.7  2011/05/03 09:39:24  abbott
# All Makefiles now offer both "clean" and "veryclean" as targets.
#
# Revision 1.6  2010/10/08 22:05:35  abbott
# Removed pointless CoCoALib prefix in definitions of CWD in Makefiles.
#
# Revision 1.5  2010/10/07 15:41:22  abbott
# Replaced explicit recursive calls to "make" bu calls to "$(MAKE)".
#
# Revision 1.4  2009/01/26 16:00:23  bigatti
# -- typo
#
# Revision 1.3  2008/12/16 10:14:13  bigatti
# -- changed makefiles for compatibility with Solaris make (no "-C" option)
#
# Revision 1.2  2007/05/14 16:37:07  bigatti
# -- removed dependencies from numerical libs
#
# Revision 1.1.1.1  2007/03/09 15:16:12  abbott
# Imported files
#
