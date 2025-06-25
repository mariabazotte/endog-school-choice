CXX	:= c++
CXXFLAGS := -Wall -Wextra -MMD -O2 -m64 -g 
INCLUDES := /home/gurobi/12.0.0/linux64/include/  
LIBSINCLUDES := /home/gurobi/12.0.0/linux64/lib/ 
LIBS = gurobi_c++ gurobi120 m

CPPFLAGS += $(foreach includedir,$(INCLUDES),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(LIBSINCLUDES),-L$(librarydir))
LDFLAGS += $(foreach library,$(LIBS),-l$(library))

OBJECTS := main.o instance/instance.o input/input.o algorithm/algorithm.o algorithm/saa.o transform/nonstrategic/nonstrategictransform.o transform/strategic/strategictransform.o transform/strategic/enum/enumtransform.o transform/strategic/compact/compacttransform.o transform/strategic/simple/simpletransform.o solver/abstractsolver.o solver/cutoffscore/cutoffscoresolver.o solver/lconstraints/lconstraintssolver.o solver/asgheuristic/asgheuristicsolver.o solver/status/status.o solver/test1st/test1stsolver.o solver/localsearch/localsearch.o solver/simulannealing/simulannealing.o solver/lagrangian/lagrangian.o solver/integerlshaped/integerlshaped.o

DIRS = $(OBJECTS:.o=.d)

EXE:=exe

$(EXE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(EXE) $(OBJECTS) $(LDFLAGS) 

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	$(RM) *.o *~ $(EXE)
	rm -rf $(OBJECTS)
	rm -rf $(DIRS)