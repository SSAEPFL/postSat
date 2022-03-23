#include <iostream>
#include "ephemeris/Ephemeris.hpp"
#include <vector>
using namespace std;

void chgmtemps(int& day, int& month, int& year, int& heure,int& minute, int& second )
{
	while(second > 60){
		second -= 60;
		minute += 1;
		while(minute > 60){
			minute -= 60;
			heure += 1;
			while (heure > 24) {
				heure -= 24;
				day += 1;
				if ((month == 1 or month ==3 or month == 5 or month ==7 or month ==8 or month == 10) && day > 31) {

						day -=31;
						month += 1;
					}
				else if ((month == 4 or month == 6 or month ==9 or month == 11) && day > 30) {

						day -=30;
						month += 1;

				}
				else if (month == 2 && day > 28) {
						day -=28;
						month += 1;
				}
				else if (month == 12 && day > 31){
					day -= 31;
					month -= 11;
					year += 1;
				}
			}
		}
	}
}
vector<double> equatorialtocartesian(double ascension,\
double declinaison, double distance){ // Transformation du système equatorial à Cartesien par rapport au point VERNAL
vector<double> array3 = vector<double>(3);
  array3[0] = distance*cos(0.2618*ascension)*cos(0.2618*declinaison); // premier composante
  array3[1] = distance*sin(0.2618*ascension)*cos(0.2618*declinaison); // deuxieme composante
  array3[2] = distance*sin(0.2618*declinaison); // troisieme composante

  return array3;
}
int main()
{
Ephemeris::setLocationOnEarth(4124,12414);
	int day=121, month=3,year=2022,hour=14,minute=24,second=01;
	SolarSystemObject planet = Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
	vector coord = equatorialtocartesian(planet.equaCoordinates.ra ,planet.equaCoordinates.dec,planet.distance*1.496e11);
    cout << "Les valeurs pour le système équatoriel vaut : " << coord[0] << " et " << coord[1]<<" et" <<coord[2]<<endl;
		month += 6;
		planet = Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
		coord = equatorialtocartesian(planet.equaCoordinates.ra ,planet.equaCoordinates.dec,planet.distance*1.496e11);
		  cout << "Les valeurs pour le système équatoriel vaut : " << coord[0] << " et " << coord[1]<<" et" <<coord[2]<<endl;
month -=6 ;
year += 1;
planet = Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
coord = equatorialtocartesian(planet.equaCoordinates.ra ,planet.equaCoordinates.dec,planet.distance*1.496e11);
  cout << "Les valeurs pour le système équatoriel vaut : " << coord[0] << " et " << coord[1]<<" et" <<coord[2]<<endl;
	month += 6;
	planet = Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
	coord = equatorialtocartesian(planet.equaCoordinates.ra ,planet.equaCoordinates.dec,planet.distance*1.496e11);
		cout << "Les valeurs pour le système équatoriel vaut : " << coord[0] << " et " << coord[1]<<" et" <<coord[2]<<endl;
		return 0;
}
