# Makefile outside of cmake

#CXX               = g++
CXX_FLAGS         = -std=c++11 -Wall -DNDEBUG -O3 --coverage
DEPS              = src/quaternion.h
OBJS              = unit_tests.o

all: unit_tests

unit_tests: $(OBJS)
	$(CXX) -o $@ $^ $(CXX_FLAGS)
	lcov --capture --directory . --output-file coverage.info

%.o: src/%.cpp $(DEPS)
	$(CXX) $(CXX_FLAGS) -c -o $@ $<

clean:
	rm *.o
