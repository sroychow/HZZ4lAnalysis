UNAME  = $(shell uname)
EXE    = crslctn
 
VPATH  = .:./interface
vpath %.h ./interface

CSUF   = cc
HSUF   = h
DICTC  = Dict.$(CSUF)
DICTH  = $(patsubst %.$(CSUF),%.h,$(DICTC))

SRCS   = src/PhysicsObjects.cc src/AnaUtil.cc src/HZZ4lUtil.cc src/AnaBase.cc src/PhysicsObjSelector.cc src/CRSelection.cc src/crSelection.cc
OBJS   = $(patsubst %.$(CSUF), %.o, $(SRCS))

LDFLAGS  = -g
SOFLAGS  = -shared 
CXXFLAGS = -I./interface -I$(CMSSW_BASE)/src/ZZMatrixElement/MEMCalculators/interface/ -I./

CXX       = g++
CXXFLAGS += -g -Wall -Wno-deprecated
LIBS =  -L$(CMSSW_BASE)/lib/$(SCRAM_ARCH) -lZZMatrixElementMEMCalculators
HDRS_DICT = interface/PhysicsObjects.h interface/LinkDef.h

bin: $(EXE) 
all: 
	make cint 
	make bin 
cling: $(DICTC) 
 
$(EXE): $(OBJS) src/Dict.o
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS) `root-config --libs`  -lTMVA
 
$(DICTC): $(HDRS_DICT)
	@echo "Generating dictionary $(DICTC) ..."
	rootcling -f $@ -rmf interface/AnalysisSpaceTreeMaker_xr.rootmap -c $(CXXFLAGS) $^
	@mv $(DICTC) src/ 
	@mv Dict_rdict.pcm src/	 
 
# Create object files
%.o : %.$(CSUF)
	$(CXX) $(CXXFLAGS) -fPIC -std=c++11 `root-config --cflags` -o $@ -c $<

# makedepend
depend: $(SRCS:.$(CSUF)=.$(CSUF).dep)
	@cat $(notdir $^) > Makefile.dep
	@-rm -f $(notdir $^) $(patsubst %,%.bak,$(notdir $^))
%.dep:
	@touch $(notdir $@)
	rmkdepend -f$(notdir $@) -- $(CXXFLAGS) `root-config --cflags` -- $*
include Makefile.dep

# Clean 
.PHONY   : clean 
clean : 
	@-rm $(OBJS) $(EXE) src/$(DICTC) src/Dict.o
