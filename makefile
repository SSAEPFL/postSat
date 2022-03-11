CXX = g++
CXXFLAGS = -g
CXXFLAGS += -pedantic
CXXFLAGS += -std=c++17
CXXFLAGS += -Wall
CXXFLAGS += -O2

simulation: simulation.o Ephemeris.o ELP2000.hpp Calendar.o configuration.in
	$(CXX) $(CXXFLAGS) $^

simulation.o: simulation.cpp Ephemeris.o ELP2000.hpp Calendar.o configuration.in
		$(CXX) $(CXXFLAGS) -c simulation.cpp -o $@

test: test.o Ephemeris.o ELP2000.hpp Calendar.o
	$(CXX) $(CXXFLAGS) $^

test.o: test.cc Ephemeris.o ELP2000.hpp Calendar.o
	$(CXX) $(CXXFLAGS) -c test.cc -o $@

Ephemeris.o: Ephemeris.cpp Ephemeris.hpp
	$(CXX) $(CXXFLAGS) -c Ephemeris.cpp -o $@

Calendar.o: Calendar.cpp Calendar.hpp
	$(CXX) $(CXXFLAGS) -c Calendar.cpp -o $@
