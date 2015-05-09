CC = g++ # c compiler
CXXFLAGS = -O2 -ansi -W -Wall -fPIC -DWITHPYTHON -c # c flags
LDFLAGS = -shared
ROOTLIBS = -lCore -lRIO -lNet -lHist -lTree -lMatrix -lMathCore -lThread -lPhysics -lPyROOT -lMinuit -lMinuit2 -lGpad -lRooFit # root libraries

#python
PYTHON_VERSION = 2.7
PYTHON_INCLUDE = $(PYTHONHOME)/include/python$(PYTHON_VERSION)
PYTHON_LIB = $(PYTHONHOME)/lib/

#boost
BOOST_INC = /usr/local/include
BOOST_LIB = /usr/local/lib
BOOSTLIBS = -lboost_python-mt
PYTHONLIBS = -lpython$(PYTHON_VERSION)

RM = rm -f # rm command
TARGET_LIB = libanalysis.so # target lib

SRCS = Expr.cxx Var.cxx Var2D.cxx Var3D.cxx EffVar.cxx EffVar2D.cxx ReweightVar.cxx Template.cxx MWTemplate.cxx Fitter.cxx Tree.cxx EntryList.cxx Eff.cxx EfficiencyClass.cxx AnalysisClass.cxx Tune.cxx Tune2D.cxx Utils.cxx pyboost.cxx

vpath %.o lib/

OBJS = $(SRCS:.cxx=.o)

all	: ${TARGET_LIB}

$(TARGET_LIB) : $(OBJS)
	$(CC) $(LDFLAGS) $^ -L$(BOOST_LIB) $(BOOSTLIBS) -L$(ROOTSYS)/lib $(ROOTLIBS) -o lib/$@ 
$(OBJS): %.o:src/%.cxx src/%.h
	$(CC) $(CXXFLAGS) -Isrc/ -I/usr/include -I$(ROOTSYS)/include -I$(PYTHON_INCLUDE) -I$(BOOST_INC) $< -o lib/$@

clean	: 
	-${RM} lib/${TARGET_LIB} ${OBJS}
