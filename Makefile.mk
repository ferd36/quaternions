# Makefile outside of cmake

BOOST_INCLUDE_DIR = /usr/local/include
BOOST_LIB_DIR     = /usr/local/lib
CXX               = g++
CXX_FLAGS         = -std=c++11 -Wall -DNDEBUG -O3
DEPS              = src/quaternion.h
OBJS              = unit_tests.o
#LIBS              = -lstdc++

all: quaternions

quaternions: $(OBJS)
	$(CXX) -o $@ $^ $(CXX_FLAGS) -L$(BOOST_LIB_DIR) $(LIBS)

%.o: src/%.cpp $(DEPS)
	$(CXX) $(CXX_FLAGS) -c -o $@ $<

clean:
	rm *.o quaternions
