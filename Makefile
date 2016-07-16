#this is much more Makefile then test3, have the root libs.
#
EXE = O
TOP = $(shell pwd)
OBJ = $(TOP)/obj
BIN = $(TOP)/bin
SRC = $(TOP)/src
INCLUDE = $(TOP)/include

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CPPLIBS = $(ROOTLIBS) $(ROOTGLIBS) -lfftw3 -lm

CPP = g++
CFLAGS = -g -Wall -fPIC -I$(INCLUDE)/ -I$(ROOTCFLAGS)


all:$(EXE)
$(EXE):$(patsubst $(SRC)/%.cpp,$(OBJ)/%.o,$(wildcard $(SRC)/*.cpp))
	$(CPP) $^ $(CPPLIBS) -o $(BIN)/$(notdir $@)			
	@echo			

$(OBJ)/%.o:$(SRC)/%.cpp
	$(CPP) $(CFLAGS) -c $(SRC)/$(notdir $<) -o $(OBJ)/$(notdir $@)
	@echo	

.PHONY:clean
clean:
	rm -f $(OBJ)/*.o $(BIN)/$(EXE) $(INCLUDE)/*h~ $(SRC)/*c~ *e~
