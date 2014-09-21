# The PROGRAM macro defines the name of the program or project.

PROGRAM = ParallelUnionFind

# The INCLUDEDIRS macro contains a list of include directories to pass to the compiler so it can find necessary header files.
# The LIBDIRS macro contains a list of library directories to pass to the linker so it can find necessary libraries.
# The LIBS macro contains a list of libraries that the the executable must be linked against.

INCLUDEDIRS = -I$(HOME)/usr/include

LIBDIRS = -L$(HOME)/usr/lib

#LIBS = -lfftw3_mpi -lfftw3

# The CXXSOURCES macro contains a list of source files.
# The CXXOBJECTS macro converts the CXXSOURCES macro into a list of object files.
# The CXXFLAGS macro contains a list of options to be passed to the compiler.  
# Adding "-g" to this line will cause the compiler to add debugging information to the executable.
# The CXX macro defines the C++ compiler.
# The LDFLAGS macro contains all of the library and library
# directory information to be passed to the linker.

CXXSOURCES = main.cpp ParallelUnionFindImpl.cpp TestParallelUnionFind.cpp ParallelUnionFind2DStripes.cpp SendColumnStrategy.cpp WeightedUnionFind.cpp ParallelUnionFind.cpp SendLeftColumnStrategy.cpp ParallelUnionFindFactory.cpp SendRightColumnStrategy.cpp
 
CXXOBJECTS = $(CXXSOURCES:.cpp=.o)  # expands to list of object files
CXXFLAGS = $(INCLUDEDIRS) #-D _DEBUG
CXX = mpic++

LDFLAGS = $(LIBDIRS) $(LIBS)

# Link target: automatically builds its object dependencies before executing its link command.

$(PROGRAM): $(CXXOBJECTS)
	$(CXX) -o $@ $(CXXOBJECTS) $(LDFLAGS)
	mv $(PROGRAM) ../Release
	rm *.o

clean:
	rm -f $(CXXOBJECTS) $(PROGRAM)
