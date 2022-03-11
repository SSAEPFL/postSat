CXX = g++
CXXFLAGS = -g
CXXFLAGS += -pedantic
CXXFLAGS += -std=c++17
CXXFLAGS += -Wall
CXXFLAGS += -O2

simulation: simulation.o Ephemeris.o ephemeris/ELP2000.hpp Calendar.o
	$(CXX) $(CXXFLAGS) $^

simulation.o: simulation.cpp Ephemeris.o ephemeris/ELP2000.hpp Calendar.o 
		$(CXX) $(CXXFLAGS) -c simulation.cpp -o $@

test: test.o Ephemeris.o ELP2000.hpp Calendar.o
	$(CXX) $(CXXFLAGS) $^

test.o: test.cc Ephemeris.o ephemeris/ELP2000.hpp Calendar.o
	$(CXX) $(CXXFLAGS) -c test.cc -o $@

Ephemeris.o: ephemeris/Ephemeris.cpp ephemeris/Ephemeris.hpp
	$(CXX) $(CXXFLAGS) -c ephemeris/Ephemeris.cpp -o $@

Calendar.o: ephemeris/Calendar.cpp ephemeris/Calendar.hpp
	$(CXX) $(CXXFLAGS) -c ephemeris/Calendar.cpp -o $@
