CXX = g++
CXXFLAGS = -g -O3 -std=c++11
LIBS =

#CXX = clang++ #-v
#CXXFLAGS = -g -O2 -std=c++11 -stdlib=libc++
#LIBS =

#CXX = icpc
#CXXFLAGS = -g -O3
#LIBS =

LD = $(CXX)
LDFLAGS = $(CXXFLAGS)

TESTS = memory_hooks

all: $(TESTS)

%: %.cxx
	$(CXX) $(CXXFLAGS) $< $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TESTS)