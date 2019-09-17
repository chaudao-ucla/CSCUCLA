FLAGS=-std=c++11 -O5 -Wall -fPIC 
#flags necessary to compile an executable which contains access to root libraries
ROOTFLAGS=`root-config --cflags --libs` -L$(ROOTSYS)/lib

LIBDIR=lib
SRCDIR=src
INCDIR=include
PROJLIBS=$(LIBDIR)/CSCClasses_cpp.so $(LIBDIR)/CSCHelperFunctions_cpp.so $(LIBDIR)/ALCTHelperFunctions_cpp.so $(LIBDIR)/LUTClasses_cpp.so $(LIBDIR)/Processor_cpp.so $(LIBDIR)/StlCollectionProxy_cpp.so

#TODO: Wildcards here!!
# Assume it contains a main() function from https://gist.github.com/ghl3/3975167
#all: $(PROJLIBS) $(SRCDIR)/PatternFinder $(SRCDIR)/printPatternCC $(SRCDIR)/CLCTLayerAnalyzer $(SRCDIR)/BayesPatternAnalysis $(SRCDIR)/testTMBEmulation $(SRCDIR)/MultiplicityStudy $(SRCDIR)/ThreeLayerCLCTEmulationAnalyzer $(SRCDIR)/LUTBuilderTEMPLATE
all: $(PROJLIBS) $(SRCDIR)/LUTBuilder $(patsubst %.cpp,%,$(wildcard $(SRCDIR)/*Tester.cpp)) $(patsubst %.cpp,%,$(wildcard $(SRCDIR)/*Analyzer.cpp)) $(patsubst %.cpp,%,$(wildcard $(SRCDIR)/*TEMPLATE.cpp)) $(SRCDIR)/PatternPrinter $(SRCDIR)/ALCTChamberPrinter $(SRCDIR)/ALCTEmulationTreeCreator


# Make shared libraries to minimize code compilation, but primarily to
# allow access of C++ classes / functions through python using
#
#  >> ROOT.gSystem.Load(<lib>)  //without .so
#
# then access with
#
#  >> ROOT.<function>(<args>)
#
$(LIBDIR)/%_cpp.so : $(SRCDIR)/%.cpp $(INCDIR)/%.h $(INCDIR)/CSCConstants.h
	##### Making shared library #####
	echo ".L $<++" | root -b -l #make the shared library using ACLiC
	mv $(SRCDIR)/*_cpp* $(LIBDIR) #move it to the libs dir
	$(shell export PROJLIBS="$(PROJLIBS)")
	$(shell export CSCPROJLIBS="$(PROJLIBS)") #for other projects
	
	
clean:
	rm $(LIBDIR)/*.so $(LIBDIR)/*.pcm $(LIBDIR)/*.d

#Makes whatever script you name, add $(ROOTLIBS) after -L$(ROOTSYS)/lib if you need it
$(SRCDIR)/%: $(SRCDIR)/%.cpp $(PROJLIBS) $(INCDIR)/CSCConstants.h
	##### Making executable #####
	$(shell export LD_LIBRARY_PATH=$(LIBDIR):${LD_LIBRARY_PATH})
	g++ $(FLAGS) $(ROOTFLAGS) -L$(LIBDIR) $(PROJLIBS) $< -o $@ 


	