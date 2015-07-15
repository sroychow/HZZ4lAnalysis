UNAME  = $(shell uname)
EXE    = crsel
 
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
CXXFLAGS = -I./interface -I./

CXX       = g++
CXXFLAGS += -g -Wall -Wno-deprecated
LIBS = 
HDRS_DICT = interface/PhysicsObjects.h interface/LinkDef.h

bin: $(EXE) 
all: 
	make cint 
	make bin 
cint: $(DICTC) 
 
$(EXE): $(OBJS) src/Dict.o
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS) `root-config --libs`  -lTMVA
 
$(DICTC): $(HDRS_DICT)
	echo "Generating dictionary $(DICTC) and $(DICTH) ..."
	rootcint -f $@ $(CXXFLAGS) $^
	perl -pi -e 's#interface/##' $(DICTH)
	mv $(DICTC) src/ 
	mv $(DICTH) interface/
 
# Create object files
%.o : %.$(CSUF)
	$(CXX) $(CXXFLAGS) -std=c++11 `root-config --cflags` -o $@ -c $<

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
	@-rm $(OBJS) $(EXE) interface/$(DICTH) src/$(DICTC) src/Dict.o
