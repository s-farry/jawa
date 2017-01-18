lib_LTLIBRARIES = libJawa.la
libJawa_la_LDFLAGS = -version-info 0:0:0
libJawa_la_SOURCES = JawaObj.cxx Expr.cxx Var.cxx Var2D.cxx Var3D.cxx EffVar.cxx EffVar2D.cxx ReweightVar.cxx Template.cxx MWTemplate.cxx Fitter.cxx Tree.cxx EntryList.cxx Eff.cxx EfficiencyClass.cxx AnalysisClass.cxx Tune.cxx Tune2D.cxx Utils.cxx
libJawa_la_CXXFLAGS = -std=c++11 @ROOTCFLAGS@
libJawa_la_LIBADD = @LIBS@ -L@ROOTLIBDIR@ @ROOTLIBS@
if WANT_PYTHON
libJawa_la_SOURCES += pyboost.cxx
libJawa_la_CXXFLAGS += @BOOST_PYTHON_CPPFLAGS@ -DWITHPYTHON
libJawa_la_LIBADD += @BOOST_PYTHON_LIBS@
endif

