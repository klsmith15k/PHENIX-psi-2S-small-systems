AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(MY_INSTALL)/include -I$(OPT_PHENIX)/include -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

AM_LDFLAGS = -L$(libdir) -L$(MY_INSTALL)/lib -L$(OPT_PHENIX)/lib -lodbc++ -L$(OFFLINE_MAIN)/lib -L$(ROOTSYS)/lib

lib_LTLIBRARIES = \
  libfill_ana_tree.la

noinst_HEADERS = \
  fill_ana_treeLinkDef.h

include_HEADERS = \
  mFillAnaTree.h 

libfill_ana_tree_la_SOURCES = \
  mFillAnaTree.cc 

nodist_libfill_ana_tree_la_SOURCES = fill_ana_tree_Dict.C

libfill_ana_tree_la_LIBADD = \
  -lmutgeom \
  -lmutoo_core \
  -lmutoo_base \
  -lmutoo_display \
  -lmutoo_interface \
  -lfvtxoo_interface \
  -lrpcoo_interface \
  -lmuioo_interface \
  -lMWG_interface \
  -lMWGOO \
  -lMWG_tools \
  -lSubsysReco \
  -lMinuit2 \
  -lphool \
	-lpicodst_object

# executables needed to make sure libraries are properly linked
noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libfill_ana_tree.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating CINT dictionaries from class headers
fill_ana_tree_Dict.C: \
  mFillAnaTree.h \
  fill_ana_treeLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
