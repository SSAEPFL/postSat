#include <iostream>
#include "ephemeris/Ephemeris.hpp"
using namespace std;

int main()
{
	Ephemeris::setLocationOnEarth(48,50,11, -2,20,14);
	int day=4, month=3,year=2022,hour=14,minute=42,second=45;
	SolarSystemObject planet = Ephemeris::solarSystemObjectAtDateAndTime(EarthsMoon, day, month, year, hour, minute, second);
    cout << "Les valeurs pour le système équatoriel vaut : " << planet.equaCoordinates.ra << " et " << planet.equaCoordinates.dec<<" et" << planet.distance<<endl;
	return 0;
}
