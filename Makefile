CC = g++ # c compiler
CXXFLAGS = -O2 -ansi -W -Wall -fPIC -c # c flags
LDFLAGS = -shared
ROOTLIBS = -lRIO -lNet -lHist -lTree -lMatrix -lMathCore -lThread -lPhysics -lPyROOT -lMinuit -lMinuit2 -lGpad -lRooFit # root libraries

#boost
BOOST_INC = /usr/local/include
BOOST_LIB = /usr/local/lib
BOOSTLIBS = -lboost_python

RM = rm -f # rm command
TARGET_LIB = libJawa.so # target lib

ROOTLIBDIR = /batchsoft/root/root53407-x86_64-slc6-46/lib
ROOTINCDIR = /batchsoft/root/root53407-x86_64-slc6-46/include
PYTHON_INCLUDE_DIR = /usr/include/python2.6

SRCS = Expr.cxx Var.cxx Var2D.cxx Var3D.cxx EffVar.cxx EffVar2D.cxx ReweightVar.cxx Template.cxx MWTemplate.cxx Fitter.cxx Tree.cxx EntryList.cxx Eff.cxx EfficiencyClass.cxx AnalysisClass.cxx Tune.cxx Tune2D.cxx Utils.cxx

SRCS += pyboost.cxx
CXXFLAGS += -DWITHPYTHON

vpath %.o lib/

OBJS = $(SRCS:.cxx=.o)

all	: ${TARGET_LIB}

$(TARGET_LIB) : $(OBJS)
	$(CC) $(LDFLAGS) $^ $(BOOSTLIBS) -L$(ROOTLIBDIR) $(ROOTLIBS) -o lib/$@ 
$(OBJS): %.o:src/%.cxx src/%.h
	$(CC) $(CXXFLAGS) -Isrc/ -I/usr/include -I$(ROOTINCDIR) -I$(PYTHON_INCLUDE_DIR) -I$(BOOST_INC) $< -o lib/$@

clean	: 
	-${RM} lib/${TARGET_LIB} ${OBJS}
