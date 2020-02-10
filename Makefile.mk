# Makefile outside of cmake

CXX_FLAGS         = -std=c++11 -Wall --coverage
DEPS              = include/quaternion.h include/quaternion_io.h include/quaternion_utils.h
OBJS              = unit_tests.o

all: unit_tests

unit_tests: $(OBJS)
	$(CXX) -o $@ $^ $(CXX_FLAGS)

%.o: test/%.cpp $(DEPS)
	$(CXX) $(CXX_FLAGS) -c -o $@ $<

clean:
	-rm unit_tests *.o *.gcda *.gcno *.gcov
