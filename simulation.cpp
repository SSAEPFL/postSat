#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include <iterator>
#include <vector>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs
                          // Fichier .tpp car inclut fonctions template
#include <array>
#include<stdio.h>
#include "ephemeris/Ephemeris.hpp"
#include <chrono>

using namespace std; // ouvrir un namespace avec la librairie c++ de base

double scalarProduct(valarray<double> const& array1,\
valarray<double> const& array2){
  // compute and return the norm2 of a valarray
  return (array1*array2).sum();
}

// definir a fonction template pour calculer la norm2 d'un valarray
template<typename T> T norm2(valarray<T> const& array){
  // compute and return the norm2 of a valarray
  return sqrt((array*array).sum());
}
// Change les mois
void chgmtmonth(int& day, int& month, int& year, int& heure,int& minute, double& second)
{
  if ((month == 1 or month ==3 or month == 5 or month ==7 or month ==8 or month == 10) && day > 31) {
    day -=31;
    month += 1;

    chgmtmonth(day,month, year, heure,minute,second);

  }
  else if ((month == 4 or month == 6 or month ==9 or month == 11) && day > 30) {
	day -=30;
	month += 1;

  chgmtmonth(day,month, year, heure,minute,second);
  }
  else if (month == 2 && day > 28 && year % 4 != 0) {
	day -=28;
	month += 1;

  chgmtmonth(day,month, year, heure,minute,second);
  }
  // Année bissextile
  else if (month == 2 && day > 28 && year % 4 == 0){
	day -=29;
	month += 1;

  chgmtmonth(day,month, year, heure,minute,second);
  }
  else if (month == 12 && day > 31){
	day -= 31;
	month -= 11;
	year += 1;

  chgmtmonth(day,month, year, heure,minute,second);
  }
  else {
    return ;
  }
}
// Convertit le surplus de seconde (i.e. plus de 60 seconde) en minute, heure, jour, mois, année
void chgmtemps(int& day, int& month, int& year, int& heure,int& minute, double& second )
{
if (second < 0){
  second += 60;
  minute -= 1;
  if(minute < 0){
    minute += 60;
    heure -= 1;
    if(heure < 0){
      heure += 24;
      day -= 1;
      }
    }
  }

  while(second >= 60){
    second -= 60;
    minute += 1;
  }

  while(minute >= 60){
	minute -= 60;
	heure += 1;

  }

  while (heure >= 24) {
    heure -= 24;
    day += 1;

  }
  /*minute += int(second) / 60;
  second = int(second) % 60 + (second - int(second));
  heure += minute / 60;
  minute = minute % 60;
  day += heure / 24;
*/
  chgmtmonth(day,month, year, heure,minute,second);

}
//Fonction rendant continu RA et DEC

template<typename T> valarray<T> vectorProduct(valarray<T> const& array1,\
valarray<T> const& array2){
  // initialiser le nouveau valarray
  valarray<T> array3=valarray<T>(3);
  // calculer le produit vecteur
  array3[0] = array1[1]*array2[2] - array1[2]*array2[1]; // premiere composante
  array3[1] = array1[2]*array2[0] - array1[0]*array2[2]; // deuxieme composante
  array3[2] = array1[0]*array2[1] - array1[1]*array2[0]; // troisieme composante
  // retourner le array3
  return array3;
}

valarray<double> equatorialtocartesian(double ascension,\
double declinaison, double distance){ // Transformation du système equatorial à Cartesien (axe X passe par le point VERNAL, axe Z passe par les pôles)
valarray<double> array3 = valarray<double>(3);

double hr2rad = 2*3.1415926535897932384626433832795028841971/24.0; // Conversion d'angle en heure -> radians
double deg2rad = 3.1415926535897932384626433832795028841971/180.0; // Conversion d'angle en degré -> radians
  array3[0] = distance*cos(hr2rad*ascension)*cos(deg2rad*declinaison); // premier composante
  array3[1] = distance*sin(hr2rad*ascension)*cos(deg2rad*declinaison); // deuxieme composante
  array3[2] = distance*sin(deg2rad*declinaison); // troisieme composante

  return array3;
}

valarray<double> cartesiantoequatorial(double const& x, double const& y, double const& z){ // Transformation de coordonnées cartésiennes en coordonnées equatoriales (axe X passe par le point VERNAL, axe z passe par les pôles)
	// Données en radian!!!
	valarray<double> equa = valarray<double>(3);
	const double pi=3.1415926535897932384626433832795028841971;
	double r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
	equa[0] = r;
	equa[1] = asin(z/r); // Déclinaison
	equa[2] = atan(y/x); // Ascension droite
	//cout << x << ' ' << y << ' ' << atan(y/x) << endl;
	if ((x*y > 0 and x < 0) or (x*y < 0 and x < 0) ){
	equa[2] +=	pi;
	}
	return equa;
}
valarray<double> continuationRADEC(double dt_,int day, int month, int year, int heure,int minute, double second  ){
chgmtemps(day,month,year,heure,minute,second);
SolarSystemObject Soleil_avant= Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, heure, minute, second);
  minute += 1;

  chgmtemps(day,month,year,heure,minute,second);
  SolarSystemObject Soleil_apres=Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, heure, minute, second);
  valarray<double> continu = valarray<double>(3);
  double ra = Soleil_avant.equaCoordinates.ra + (-Soleil_avant.equaCoordinates.ra+Soleil_apres.equaCoordinates.ra)/50.0 * dt_;
  double dec =  Soleil_avant.equaCoordinates.dec + (-Soleil_avant.equaCoordinates.dec+Soleil_apres.equaCoordinates.dec)/50.0 * dt_;

  continu = equatorialtocartesian(ra, dec, 1.496e11*Soleil_avant.distance);

  return continu;
}
/* La class Engine est le moteur principale de ce code. Il contient
   les methodes de base pour lire / initialiser les inputs,
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{
/*ifstream is("Harris_Priester_tab.txt");
	istream_iterator<double> start(is), end;
	vector<double> numbers(start, end); // Liste des coefficients (0 mod 3 -> altitude, 1 mod 3 -> rho_m , 2 mod 3 -> rho_M)
	is.close();
	valarray<double> coeff(numbers.size()); // Conversion d'un vector en valarray
	for(size_t i(0); i < coeff.size();i++){
			coeff[i] = numbers[i];
		}
	size_t N = coeff.size()/3;
		// Réorganisation en 3 valarray
	valarray<double> h_ = coeff[slice(0,N,3)];
	valarray<double> rho_m = coeff[slice(1,N,3)];
	valarray<double> rho_M = coeff[slice(2,N,3)];*/
private:
  // definition des constantes
  const double pi=3.1415926535897932384626433832795028841971e0;
  const double G=6.674e-11;
  const double masse_soleil =1.988435e30 ;
  const double masse_lune = 7.3459e22;
  const double masse_terre = 5.97e24;
  const double rayon_terre = 6378.13630e3;
  const double celeritas =  299792458;
  const double Gm_t = 398600.4415e9;

  // definition des variables
  double tfin=0.e0;      // Temps final
  unsigned int nsteps=1; // Nombre de pas de temps

  valarray<double> x0=valarray<double>(0.e0,6); // vecteur contenant la position et vitesse initiale du Satellite en
			 		        // utilisant la sequence index-valeur: 0-2 : position 3-5 vitesse

  unsigned int sampling=1; // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie
  // Coefficient pour la densité atmosphérique (Harris-Priester)
  valarray<double> h_ = valarray<double>(0.,50);
  valarray<double> rho_m = valarray<double>(0.,50);
  valarray<double> rho_M = valarray<double>(0.,50);
  // Coefficient pour le géopotentiel
  vector<vector<double>> C_nm;
  vector<vector<double>> S_nm;
  unsigned int ordre;
  double dr; // Valeur de la quantité infinitésimal dr pour calculer le gradient
  double dphi; 
  void printOut(bool write,vector< vector<double> >& matrix)
  {
    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {

   //    *outputFile << t  << " "<<x[0]<< " "<< x[1] << " " << x[2]  << " "<<dt<< " "<<endl; // write output on file

    matrix[0].push_back(t);
    matrix[1].push_back(x[0]);
    matrix[2].push_back(x[1]);
    matrix[3].push_back(x[2]);
    matrix[4].push_back(x[3]);
    matrix[5].push_back(x[4]);
    matrix[6].push_back(x[5]);
    matrix[7].push_back(dt);


      last = 1;
       // fin +=1;
    }
    else
    {
      last++;
    }
    }
protected:
  // Iteration temporelle, a definir au niveau des classes filles
  virtual void step(valarray<double>&,valarray<double>&, double)=0;
  virtual void step2(double&)=0;
  double mass ;  // masses des corps
  SolarSystemObject Soleil;
  SolarSystemObject Lune ;
int day,month,year,hour,minute;
double second;
const double astronomical_unit= 1.496e11;
double Solar_area; // Aire du satellite recevant la lumière du soleil
double Drag_area; // Surface de frottement du satellite avec l'atmosphère
double C_d; // Coefficient de frottement (Drag coefficient)
valarray<double> x1=valarray<double>(0.e0,12); // Position des astres
// donnes internes
  double t,dt,tol;        // Temps courant pas de temps
  valarray<double> x=valarray<double>(0.e0,6); // Position des astres
  template<typename T> valarray<T> heliocentricshift(valarray<double> position)
  {
    valarray<double> shifting = valarray<double>(0.e0,3);
    //shifting = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance);
    shifting += position;
    return shifting;
  }

valarray<double> distance(int j, valarray<double> const& x_,valarray<double> x1_) const { //  On donne les deux vecteurs de postions et l'indice de l'astre (1 : Soleil, 2: Lune, 3: Terre, 4: Autre)
   // La distance est toujours par rapport au Satellite
  valarray<double> diff= valarray<double>(0.e0,3);
  valarray<double> astre= valarray<double>(0.e0,12);
  astre = x1_;
  int d = 3*(j-1);

  diff[0] = x_[0]-x1_[d];
  diff[1] = x_[1]-x1_[d+1];
  diff[2] = x_[2]-x1_[d+2];

  return diff;
}

double altitude() const{
	return norm2(distance(3,x,x1))-rayon_terre;
}
valarray<double> ForceGravitationSoleil(valarray<double> const& x_,valarray<double> const& x1_) const {
  valarray<double> force= valarray<double>(1.e0,3);

  force *= -G*masse_soleil;

  valarray<double> soleil= x1_[slice(0,3,1)];
  double norme3 = norm2(soleil)*norm2(soleil)*norm2(soleil);
  force[0] = force[0] * ((x_[0]-x1_[0])/(norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_)))+x1_[0]/(norme3));
  force[1] = force[1] * ((x_[1]-x1_[1])/(norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_)))+x1_[1]/(norme3));
  force[2] = force[2] * ((x_[2]-x1_[2])/(norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_)))+x1_[2]/(norme3));

/*
force[0] = force[0] * ((x_[0]-x1_[0])/(norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_))));
force[1] = force[1] * ((x_[1]-x1_[1])/(norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_))));
force[2] = force[2] * ((x_[2]-x1_[2])/(norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_))));
*/

  return force;
}

valarray<double> ForceGravitationLune(valarray<double> const& x_,valarray<double> const& x1_) const {

  valarray<double> force= valarray<double>(0.e0,3);

valarray<double> lune = x1_[slice(3,3,1)];
double norme3 = norm2(lune)*norm2(lune)*norm2(lune);
  force = -G*masse_lune/(norm2(distance(2,x_,x1_))*norm2(distance(2,x_,x1_))*norm2(distance(2,x_,x1_)));
  force[0] = force[0] * (x_[0]-x1_[3]);
  force[1] = force[1] * (x_[1]-x1_[4]);
  force[2] = force[2] * (x_[2]-x1_[5]);
  force[0] -= G*masse_lune*x1_[3]/norme3;
  force[1] -= G*masse_lune*x1_[4]/norme3;
  force[2] -= G*masse_lune*x1_[5]/norme3;

    return force;
}

valarray<double> ForceGravitationTerre(valarray<double> const& x_,valarray<double> const& x1_) const {

  valarray<double> force= valarray<double>(0.e0,3);
  force = G*masse_terre/(norm2(distance(3,x_,x1_))*norm2(distance(3,x_,x1_))*norm2(distance(3,x_,x1_)));
  force[0] = force[0] * (x1_[6] - x_[0]);
  force[1] = force[1] * (x1_[7] - x_[1]);
  force[2] = force[2] * (x1_[8] - x_[2]);

	return force;
}

double P_norm(size_t const& n, size_t const& m, long double const& x)const{

	// Polynome de Legendre associé normalisé


	double k = 0;

	if (m == 0){
		k = 1;
	}
	double C = sqrt((2-k)*(2*n+1)*tgamma(n-m+1)/tgamma(n+m+1));
	//cout << C <<' ' << m << ' ' << n <<endl;
	//cout << setprecision(30) << C*assoc_legendrel(n,m,x) << endl;
	return C*assoc_legendrel(n,m,x);
	//C*assoc_legendrel(n,m,x);
}

double Geopot(valarray<double> const& equatorial) const {

	// Calcul le potentiel terrestre U. Les coordonées sont dans le système equatorial
	// Formules tirées du livre Satellite Orbits (Montenbruck and Gill) p.57-58

	double r = equatorial[0];
	double phi = equatorial[1];
	double lambda = equatorial[2];
	//cout << r <<' '<< phi << ' ' << lambda << ' ' <<endl;
	double double_somme(0);

	for(size_t n(0); n <= ordre ; n++){
		for(size_t m(0); m <= n ; m++){
			double_somme += pow(rayon_terre/r,n)*P_norm(n,m,sin(phi))*(C_nm[n][m]*cos(m*lambda)+S_nm[n][m]*sin(m*lambda));
			//cout << "xd ";
			//cout << pow(rayon_terre/r,n) <<' ' <<P_norm(n,m,sin(phi)) << ' ' <<C_nm[n][m]*cos(m*lambda) <<' ' << S_nm[n][m]*sin(m*lambda) << endl;
		}}
		//cout << endl;
	//cout << pow(rayon_terre/r,n)*P_norm(n,m,sin(phi)) << ' ' << (C_nm[n][m]*cos(m*lambda)+S_nm[n][m]*sin(m*lambda)) <<endl;
	//cout << double_somme
	//cout << double_somme<< endl;
	double U = G*masse_terre/r*double_somme;

	return U;
}
double dUdr_geo(valarray<double> const& equatorial) const{
	
	double dUdr;
	
	double r = equatorial[0];
	double phi = equatorial[1];
	double lambda = equatorial[2];
	double double_somme(0);
	// 
	
	for(size_t n(0); n <= ordre ; n++){
		for(size_t m(0); m <= n ; m++){
			double_somme += (n+1)*pow(rayon_terre/r,n)*P_norm(n,m,sin(phi))*(C_nm[n][m]*cos(m*lambda)+S_nm[n][m]*sin(m*lambda));
		}}
		
	dUdr = -Gm_t/pow(r,2)*double_somme;
	
	return dUdr;
}

double dUdphi_geo(valarray<double> const& equatorial) const{
	
	double dUdphi;
	
	double r = equatorial[0];
	double phi = equatorial[1];
	double lambda = equatorial[2];
	double double_somme(0);
	
	for(size_t n(0); n <= ordre ; n++){
		for(size_t m(0); m <= n ; m++){
			double dP_norm;
			dP_norm = (P_norm(n,m,sin(phi + dphi)) - P_norm(n,m,sin(phi)))/dphi;
			double_somme += pow(rayon_terre/r,n)*dP_norm*cos(phi)*(C_nm[n][m]*cos(m*lambda)+S_nm[n][m]*sin(m*lambda));
		}}
	dUdphi = Gm_t/r*double_somme;
	return dUdphi;
}

double dUdlambda_geo(valarray<double> const& equatorial)const{
	
	double dUdlambda;
	
	double r = equatorial[0];
	double phi = equatorial[1];
	double lambda = equatorial[2];
	double double_somme(0);
	
	for(size_t n(0); n <= ordre ; n++){
		for(size_t m(0); m <= n ; m++){
			double_somme += pow(rayon_terre/r,n)*m*P_norm(n,m,sin(phi))*(S_nm[n][m]*cos(m*lambda)-C_nm[n][m]*sin(m*lambda));
		}}
		
	dUdlambda = Gm_t/r*double_somme;
	return dUdlambda;
}
valarray<double> Acceleration_Geopotentiel(valarray<double> const& x_, valarray<double> const& x1_) const{

	valarray<double> x_vect = x_[slice(0,3,1)];
	valarray<double> eq_vect(3);
	eq_vect = cartesiantoequatorial(x_vect[0],x_vect[1],x_vect[2]);

	// Quantités infinitésimales

	valarray<double> vec_dr(3); // Vecteur infinitésimal
	vec_dr[0] = dr;
	vec_dr[1] = 0;
	vec_dr[2] = 0;

	valarray<double> vec_dphi(3);
	vec_dphi[0] = 0;
	vec_dphi[1] = dphi;
	vec_dphi[2] = 0;

	valarray<double> dlambda(3);
	dlambda[0] = 0;
	dlambda[1] = 0;
	dlambda[2] = dphi;

	double U = Geopot(eq_vect);

	// Définition de grandeurs pratiques

	double x = x_vect[0];
	double y = x_vect[1];
	double z = x_vect[2];
	double r = eq_vect[0];

	// Dérivées de r

	double drdx = x/r;
	double drdy = y/r;
	double drdz = z/r;

	// Dérivées de phi

	double dphidx = -x/(pow(r,3)*sqrt(1-pow(z/r,2)));
	double dphidy = -y/(pow(r,3)*sqrt(1-pow(z/r,2)));
	double dphidz = 2*z/r/sqrt(1-pow(z/r,2))*(1/r - pow(z,2)/pow(r,3));

	// Dérivées de lambda

	double dlambdadx = -y/(pow(x,2) + pow(y,2));
	double dlambdady = x/(pow(x,2) + pow(y,2));
	double dlambdadz = 0;

	// Calcul du gradient

	//cout << Geopot(eq_vect + vec_dr) << ' ' << U << endl;
	//cout << eq_vect[0] << ' ' << eq_vect[1] << ' '  << eq_vect[2] << endl;
	//cout << vec_dr[0] << ' ' << vec_dr[1] << ' ' << vec_dr[2] << endl;

	/*
	double dUdr = (Geopot(eq_vect + vec_dr)-U)/dr;
	double dUdphi = (Geopot(eq_vect + vec_dphi)-U)/vec_dphi[1];
	double dUdlambda = (Geopot(eq_vect + dlambda)-U)/dlambda[2];
	*/
	
	double dUdr = dUdr_geo(eq_vect);
	double dUdphi = dUdphi_geo(eq_vect);
	double dUdlambda = dUdlambda_geo(eq_vect);

	//cout << dUdr << ' ' << dUdphi << ' ' << dUdlambda << ' ' << drdx << ' ' << drdy << ' ' << drdz << ' ' << dphidx << ' ' << dphidy << ' ' << dphidz << ' ' << dlambdadx << ' ' <<dlambdady << ' ' << dlambdadz <<endl;

	double dUdx = dUdr*drdx + dUdphi*dphidx + dUdlambda*dlambdadx;
	double dUdy = dUdr*drdy + dUdphi*dphidy	+ dUdlambda*dlambdady;
	double dUdz = dUdr*drdz + dUdphi*dphidz + dUdlambda*dlambdadz;

	valarray<double> r_p_p = valarray<double> (0.e0,3);
	r_p_p[0] = dUdx;
	r_p_p[1] = dUdy;
	r_p_p[2] = dUdz;

	valarray<double> f(3);
	f = ForceGravitationTerre(x_,x1_);

	//Pour débuguer

	//cout << f[0] << ' ' << f[1] << ' ' << f[2] << endl;
	//cout << r_p_p[0] << ' ' << r_p_p[1] << ' ' << r_p_p[2] << endl;

	return r_p_p;
}

valarray<double> Acceleration_Geopotentiel2(valarray<double> const& x_, valarray<double> const& x1_)const{
	
	// Juste pour ordre 2 en calculant explicitement les dérivées des polynômes de Legendre
	
	valarray<double> x_vect = x_[slice(0,3,1)];
	valarray<double> eq_vect(3);
	eq_vect = cartesiantoequatorial(x_vect[0],x_vect[1],x_vect[2]);
	
	

	// Quantités infinitésimales

	valarray<double> vec_dr(3); // Vecteur infinitésimal
	vec_dr[0] = dr;
	vec_dr[1] = 0;
	vec_dr[2] = 0;

	valarray<double> vec_dphi(3);
	vec_dphi[0] = 0;
	vec_dphi[1] = dphi;
	vec_dphi[2] = 0;

	valarray<double> dlambda(3);
	dlambda[0] = 0;
	dlambda[1] = 0;
	dlambda[2] = dphi;

	double U = Geopot(eq_vect);

	// Définition de grandeurs pratiques

	double x = x_vect[0];
	double y = x_vect[1];
	double z = x_vect[2];
	double r = eq_vect[0];
	double phi = eq_vect[1];
	double lambda = eq_vect[2];

	// Dérivées de r

	double drdx = x/r;
	double drdy = y/r;
	double drdz = z/r;

	// Dérivées de phi

	double dphidx = -x/(pow(r,3)*sqrt(1-pow(z/r,2)));
	double dphidy = -y/(pow(r,3)*sqrt(1-pow(z/r,2)));
	double dphidz = 2*z/r/sqrt(1-pow(z/r,2))*(1/r - pow(z,2)/pow(r,3));

	// Dérivées de lambda

	double dlambdadx = -y/(pow(x,2) + pow(y,2));
	double dlambdady = x/(pow(x,2) + pow(y,2));
	double dlambdadz = 0;

	// Calcul du gradient

	/*
	double dUdr = (Geopot(eq_vect + vec_dr)-U)/dr;
	double dUdphi = (Geopot(eq_vect + vec_dphi)-U)/vec_dphi[1];
	double dUdlambda = (Geopot(eq_vect + dlambda)-U)/dlambda[2];
	*/
	
	double dUdr = dUdr_geo(eq_vect);
	double dUdlambda = dUdlambda_geo(eq_vect);
	
	// Calcul de dUdphi "à la main"
	double dUdphi;
	double sum_n2(0);
	double terme_20_phi = sqrt(5)*3*sin(phi)*cos(phi)*(-484.165368);
	double terme_21_phi = sqrt(5/3)*3*cos(2*phi)*(-0.000187*cos(lambda) + 0.001195*sin(lambda));
	double terme_22_phi = sqrt(5/24)*(-6*sin(phi))*cos(phi)*(2.439261*cos(2*lambda) + -1.400266*sin(2*lambda));
	sum_n2 = pow(rayon_terre/r,2)*(terme_20_phi + terme_21_phi + terme_22_phi)*1.e-6;
	dUdphi = G*masse_terre/r*sum_n2;

	//cout << dUdr << ' ' << dUdphi << ' ' << dUdlambda << ' ' << drdx << ' ' << drdy << ' ' << drdz << ' ' << dphidx << ' ' << dphidy << ' ' << dphidz << ' ' << dlambdadx << ' ' <<dlambdady << ' ' << dlambdadz <<endl;

	double dUdx = dUdr*drdx + dUdphi*dphidx + dUdlambda*dlambdadx;
	double dUdy = dUdr*drdy + dUdphi*dphidy	+ dUdlambda*dlambdady;
	double dUdz = dUdr*drdz + dUdphi*dphidz + dUdlambda*dlambdadz;

	valarray<double> r_p_p = valarray<double> (0.e0,3);
	r_p_p[0] = dUdx;
	r_p_p[1] = dUdy;
	r_p_p[2] = dUdz;

	valarray<double> f(3);
	f = ForceGravitationTerre(x_,x1_);

	//Pour débuguer

	//cout << f[0] << ' ' << f[1] << ' ' << f[2] << endl;
	//cout << r_p_p[0] << ' ' << r_p_p[1] << ' ' << r_p_p[2] << endl;

	return r_p_p;
}

valarray<double> ForceFrottement(valarray<double> const& x_,valarray<double> const& x1_) const {
	// Force de frottement atmosphérique pour des satellites
	// x_ : Position du satellite (taille 6, 3 pos + 3 vit), x1_ : Positions des astres (taille 12, 3*4 pos)

  valarray<double> e_v(3);  // Vecteur vitesse normalisé du satellite
	valarray<double> x_vect = x_[slice(0,3,1)];
	valarray<double> v_vect = x_[slice(3,3,1)];
	double v = norm2(v_vect);
	e_v = v_vect/v;
	valarray<double> omega_T(3); // Vecteur vitesse angulaire terrestre. ATTENTION! Axe z le long de l'axe de rotation de la Terre (En fait c'est bien)
	omega_T[0] = 0;
	omega_T[1] = 0;
	omega_T[2] = 0.7292e-4; // rad/s
	valarray<double> v_r = valarray<double>(0.e0,3); // Vecteur vitesse relative du satellite par rapport à celle de l'atmosphère
	v_r = v_vect - vectorProduct(omega_T,x_vect);
	valarray<double> force =-0.5*C_d*Drag_area*rho(x_,x1_,6)*pow(norm2(v_r),2)*e_v;
	//cout << e_v[0] << ' '  << e_v[1] << ' ' << e_v[2] <<endl;
	//cout << force[0] << ' '  << force[1] << ' ' << force[2] <<endl;
	//cout << norm2(force) << endl;
	//cout << pow(norm2(v_r),2)<< endl;
	//cout << rho(x_,x1_,6) << endl;

  //valarray<double> force = valarray<double>(0.e0,3);
	return force;
}
double rho(valarray<double> const& x_, valarray<double> const& x1_,int const& n) const {
	// Implémentation du modèle Harris-Priester pour décrire la pression atmosphérique entre 100km et 1000km d'altitude.
	//cout << "lol" << endl;
	size_t i = 0;
	double h = altitude();
	if(h <= 1000.e3){
		while (h_[i+1]*1000. <= h){
			i++;
		}
	}else{return 0;}
	double H_m = (h_[i] - h_[i+1])/(log(rho_m[i+1]/rho_m[i]));
	double H_M = (h_[i] - h_[i+1])/(log(rho_M[i+1]/rho_M[i]));
	double rho_m_h = rho_m[i]*exp((h_[i]-h/1000.)/H_m);
	double rho_M_h = rho_M[i]*exp((h_[i]-h/1000.)/H_M);
	//cout << H_m << ' '  << H_M << ' ' <<  rho_m[i] <<  ' '  << rho_M[i] << ' ' << rho_m_h << ' ' << rho_M_h << endl;

	valarray<double> e_r = distance(3,x_,x1_)/norm2(distance(3,x_,x1_)); // Vecteur position Terre->Satellite normalisé
	valarray<double> r_T_S = x1_[slice(0,3,1)]; // Vecteur position Terre->Soleil
	valarray<double> e_r_sol = r_T_S/norm2(r_T_S);

	/*Valeurs des cosinus et sinus de la position du soleil
	 * cos_dec_cos_asc = r_T_S[0];
	cos_dec_sin_asc = r_T_S[1];
	sin_dec = r_T_S[2];*/

	double lambda = 30*3.1415926535897932384626433832795028841971/180.0;

	valarray<double> e_b(3);
	e_b[0] = e_r_sol[0] * cos(lambda) - e_r_sol[1] * sin(lambda);
	e_b[1] = e_r_sol[1] * cos(lambda) + e_r_sol[0] * sin(lambda);
	e_b[2] = e_r_sol[2];
	double cos_psi = sqrt(0.5*(1 + scalarProduct(e_r,e_b)));
	//cout << rho_m_h + (rho_M_h - rho_m_h)*pow(cos_psi,n) << endl;
	double rho = (rho_m_h + (rho_M_h - rho_m_h)*pow(cos_psi,n))*1.e-12;
	return rho;
 }
valarray<double> ForceSolaire(valarray<double> const& x_,valarray<double> const& x1_) const {

valarray<double> force= valarray<double>(0.e0,3);
// Force de radiation :
// Cst Solar formula : sigma T^4 (R_s/R)^2
double sigma = 5.670374e-8; double T_s = 5778; double R_s = 6.957e8;
valarray <double> position = x_[slice(0,3,1)];
double solarCst = sigma*pow(T_s,4)*R_s/norm2(position)*R_s/norm2(position);

force = solarCst/celeritas * Solar_area;
force[0] = x_[0]/norm2(position);
force[1] = x_[1]/norm2(position);
force[2] = x_[2]/norm2(position);

return force;
}
valarray<double> ForceCoriolis(valarray<double> const& x_,valarray<double> const& x1_,double dt_, double second, int minute, int heure, int jour, int mois, int annee) const {

  valarray<double> vitesse_terre = valarray<double>(0.e0,3);
  SolarSystemObject Soleil= Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
//  vitesse_terre = (continuationRADEC(dt_, jour, mois, annee, heure, minute, second) - x1_[slice(0, 3, 1)])/ dt_;
vitesse_terre = (equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)- x1_[slice(0, 3, 1)])/ dt_;
  valarray<double> vitesser = x_[slice(3,3,1)];

  // -2m \omega x vitessesatellite
  double distance = sqrt(x1_[0] * x1_[0] + x1_[1] * x1_[1] + x1_[2] * x1_[2]);
 double vitesseangulaire = norm2(vitesse_terre)/distance;
 //double vitesseangulaire = 1;
  valarray<double> ez = valarray<double>(0.e0,3);
  ez[2] = vitesseangulaire;
//ez[2] /= distance;
  valarray<double> resultante = valarray<double>(0.e0, 3);
 resultante = vectorProduct(ez, vitesser);

  return -2.0*resultante;
}

valarray<double> ForceCentrifuge(valarray<double> const& x_,valarray<double> const& x1_,double dt_, double second, int minute, int heure, int jour, int mois, int annee) const {
  // OK
  // \omega cross \omega cross position

  valarray<double> vect_der = valarray<double>(0.e0,3);
  valarray<double> position = valarray<double>(0.e0,3);
  valarray<double> resultante = valarray<double>(0.e0,3);

  position = x_[slice(0,3,1)];
  //vect_der = (continuationRADEC(dt_,jour, mois, annee, heure, minute,second+dt_) - x1_[slice(0,3,1)]) / dt_;
  vect_der = (equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)- x1_[slice(0, 3, 1)])/ dt_;

  // -2m \omega x vitessesatellite
  double distance = sqrt(x1_[0]*x1_[0] + x1_[1]*x1_[1]+x1_[2]*x1_[2]);
  double vitesseangulaire = norm2( vect_der)/distance;
//  vitesseangulaire = 29800/distance;
  valarray<double> ez = valarray<double>(0.e0,3);
  ez[2] = 2*-vitesseangulaire;
  resultante = vectorProduct(ez,position);
  resultante = vectorProduct(resultante,ez);

  return resultante;
}

valarray<double> ForceEuler(valarray<double> const& x_,valarray<double> const& x1_,double dt_, double second, int minute, int heure, int jour, int mois, int annee) const {

  valarray<double> vect_der = valarray<double>(0.e0,3);
  return vect_der;
  valarray<double> position = x_[slice(0,3,1)];

  // A optimiser
  if (t>2*dt_) {
    vect_der = continuationRADEC(dt_,jour, mois, annee, heure, minute,second+dt_);
    vect_der += continuationRADEC(dt_,jour, mois, annee, heure, minute,second-dt_)  ;
    vect_der -= 2.0*continuationRADEC(dt_,jour, mois, annee, heure, minute,t);
    vect_der /= dt_*dt_;
  }
//  cout << "2e der " << norm2(vect_der) << endl; // problème est là, l'accélération angulaire est beaucoup trop grande !!

  vect_der = -vectorProduct(vect_der,position);

  return vect_der;
}


valarray<double> acceleration(valarray<double> const& x_,valarray<double> const& x1_,double dt_, double second, int minute, int heure, int jour, int mois, int annee) const{
  valarray<double> accelere(0.e1,3);
  /*accelere[0] = ForceGravitationSoleil(x_,x1_)[0]+ForceGravitationTerre(x_,x1_)[0]+ForceGravitationLune(x_,x1_)[0]+ForceFrottement(x_,x1_)[0]/mass+ForceSolaire(x_,x1_)[0]/mass ;
  accelere[1] = ForceGravitationSoleil(x_,x1_)[1]+ForceGravitationTerre(x_,x1_)[1]+ForceGravitationLune(x_,x1_)[1]+ForceFrottement(x_,x1_)[1]/mass+ForceSolaire(x_,x1_)[1]/mass;
  accelere[2] = ForceGravitationSoleil(x_,x1_)[2]+ForceGravitationTerre(x_,x1_)[2]+ForceGravitationLune(x_,x1_)[2]+ForceFrottement(x_,x1_)[2]/mass+ForceSolaire(x_,x1_)[2]/mass;
*/

//accelere = ForceGravitationSoleil(x_,x1_)+ForceGravitationTerre(x_,x1_)+ForceGravitationLune(x_,x1_)+ForceFrottement(x_,x1_)/mass+ForceSolaire(x_,x1_)/mass;
accelere = ForceGravitationSoleil(x_,x1_)+Acceleration_Geopotentiel(x_,x1_)+ForceGravitationLune(x_,x1_)+ForceFrottement(x_,x1_)/mass+ForceSolaire(x_,x1_)/mass;
//accelere = ForceGravitationSoleil(x_,x1_)+ForceGravitationTerre(x_,x1_)+ForceFrottement(x_,x1_)/mass;
//accelere = ForceGravitationTerre(x_,x1_) + ForceFrottement(x_,x1_)/mass;
//accelere = ForceGravitationTerre(x_,x1_) + ForceFrottement(x_,x1_)/mass+ ForceGravitationSoleil(x_,x1_)+ForceGravitationLune(x_,x1_)+ForceSolaire(x_,x1_)/mass+ForceCoriolis(x_,x1_,dt_,second,minute, heure,jour, mois, annee)+ForceCentrifuge(x_,x1_,dt_,second,minute, heure,jour, mois, annee)+ForceEuler(x_,x1_,dt_,second,minute, heure,jour, mois, annee);
//accelere = ForceGravitationTerre(x_,x1_)+ForceGravitationSoleil(x_,x1_)+ForceGravitationLune(x_,x1_)+ForceFrottement(x_,x1_)/mass+ForceSolaire(x_,x1_)/mass;
//cout <<"Terre2 : " << Acceleration_Geopotentiel(x_,x1_)[0] <<"Terre : " << ForceGravitationTerre(x_,x1_)[0] <<" Lune : " << ForceGravitationLune(x_,x1_)[0]<< " Soleil : " << ForceGravitationSoleil(x_,x1_)[0]<< " Frottement " << ForceFrottement(x_,x1_)[0]/mass << " Inertie : " <<ForceCoriolis(x_,x1_,dt_,second,minute, heure,jour, mois, annee)[0]+ForceCentrifuge(x_,x1_,dt_,second,minute, heure,jour, mois, annee)[0]+ForceEuler(x_,x1_,dt_,second,minute, heure,jour, mois, annee)[0] << " Radiation : " <<ForceSolaire(x_,x1_)[0]/mass << endl;
  return accelere;
}



public:

  valarray<double> fonction(valarray<double>const& x_,valarray<double> const& x1_,double dt_, double second,int minute, int heure, int jour, int mois, int annee) const {
      valarray<double> f = valarray<double>(0.e0, 6);
      valarray<double> accelere = valarray<double>(0.e0, 3);
      accelere = acceleration(x_,x1_,dt_,second,minute, heure,jour, mois, annee);
      for (int i(0); i<3; i++){
        f[i] = x_[i+3];
        f[i+3] = accelere[i];
  }
      return f;
    }
  /* Constructeur de la classe Engine
     inputs:
       configFile: (ConfigFile) handler du fichier d'input
  */

  Engine(ConfigFile configFile)
  {
    // variable locale

    // Stockage des parametres de simulation dans les attributs de la classe
    tfin     = configFile.get<double>("tfin");		 // lire la temps totale de simulation
    nsteps   = configFile.get<unsigned int>("nsteps");   // lire la nombre de pas de temps
    mass     = configFile.get<double>("mass1");		 // lire la mass du corps 1
    x0[0]    = configFile.get<double>("x01");		 // lire composante x position initiale Satellite
    x0[1]    = configFile.get<double>("y01");		 // lire composante y position initiale Satellite
    x0[2]    = configFile.get<double>("z01");		 // lire composante z position initiale Satellite
    x0[3]    = configFile.get<double>("vx01");		 // lire composante x vitesse initiale Satellite
    x0[4]    = configFile.get<double>("vy01");		 // lire composante y vitesse initiale Satellite
    x0[5]    = configFile.get<double>("vz01");		 // lire composante z vitesse initiale Satellite
    //Solar_area   = configFile.get<double>("Solar_area");		 // lire composante de l'aire
   // Drag_area = configFile.get<double>("Drag_area"); // lire la surface de frottement

    C_d = configFile.get<double>("C_d");             // lire le drag coefficient
    day = configFile.get<double>("day"); // Lire l'heure de la détection
    month = configFile.get<double>("month"); // Lire l'heure de la détection
    year = configFile.get<double>("year"); // Lire l'heure de la détection
    hour = configFile.get<double>("hour"); // Lire l'heure de la détection
    minute = configFile.get<double>("minute"); // Lire l'heure de la détection
    second = configFile.get<double>("second"); // Lire l'heure de la détection
    sampling = configFile.get<unsigned int>("sampling"); // lire le parametre de sampling
    tol = configFile.get<double>("tol");
    ordre = configFile.get<unsigned int>("ordre");

    dt = tfin / nsteps;          // calculer le time step
    
    dr = configFile.get<double>("dr"); // Quantité infinitésimal dr
    dphi = configFile.get<double>("dphi"); // Quantité infinitésimal dphi


    /*Soleil= Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
    Lune = Ephemeris::solarSystemObjectAtDateAndTime(EarthsMoon, day, month, year, hour, minute, second); // initialisation du soleil et de la Lune
*/
	// Création de la table de coefficients pour rho (H-P)
	ifstream is("Harris_Priester_tab.txt");
	istream_iterator<double> start(is), end;
	vector<double> numbers(start, end); // Liste des coefficients (0 mod 3 -> altitude, 1 mod 3 -> rho_m , 2 mod 3 -> rho_M)
	is.close();
	valarray<double> coeff(numbers.size()); // Conversion d'un vector en valarray
	for(size_t i(0); i < coeff.size();i++){
			coeff[i] = numbers[i];
		}
	size_t N = coeff.size()/3;
		// Réorganisation en 3 valarray
	h_ = coeff[slice(0,N,3)];
	rho_m = coeff[slice(1,N,3)];
	rho_M = coeff[slice(2,N,3)];

	// Création des tables de coefficients pour le géopotentiel

	ifstream cs("C_nm.txt");
	istream_iterator<double> start_c(cs), end_c;
	vector<double> numbers_c(start_c, end_c);
	cs.close();
	size_t n(ordre);
	vector<vector<double>> tab_C(n+1);
	for (size_t i(0) ; i <= n; i++){
		for (size_t j(0) ; j <= i ; j++){
			tab_C[i].push_back(1.e-6*numbers_c[i*(i+1)/2 + j]);
			}
		}

	C_nm = tab_C;

	ifstream ss("S_nm.txt");
	istream_iterator<double> start_s(ss), end_s;
	vector<double> numbers_s(start_s, end_s);
	ss.close();
	vector<vector<double>> tab_S(n+1);
	for (size_t i(0) ; i <= n; i++){
		for (size_t j(0) ; j <= i ; j++){
			tab_S[i].push_back(1.e-6*numbers_s[i*(i+1)/2 + j]);
			}
		}
	S_nm = tab_S;

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(20); // Les nombres seront ecrits avec 15 decimales
  };

  // Destructeur virtuel
  virtual ~Engine()
  {
    cout <<hour<< " "<< minute << " " <<second << " " << day << " "<< month << " " << year <<endl; // write output on file
    outputFile->close();
    delete outputFile;

  };
  // Simulation complete
  void run()
  {
    t = 0.e0; // initialiser le temps

    Ephemeris::setLocationOnEarth(46.61910958572537, 6.220300715342668);
   Soleil= Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
  Lune = Ephemeris::solarSystemObjectAtDateAndTime(EarthsMoon, day, month, year, hour, minute, second); // initialisation du soleil et de la Lune
  valarray<double> positionLune = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance); // Géocentrique (Terre -> Lune)
  x = x0; // Initialisation du satellite héliocentrique
    // Initialisation des positions
    x1[6]    = 0.0;		 // lire composante x position initiale Terre
    x1[7]    = 0.0;		 // lire composante y position initiale Terre
    x1[8]    = 0.0;		 // lire composante z position initiale Terre
    x1[0]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[0];		 // lire composante x position initiale Corps Soleil
    x1[1]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[1];		 // lire composante y position initiale Corps Soleil
    x1[2]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[2];	 // lire composante z position initiale Corps Soleil
    x1[3]    = positionLune[0] ;		 //  Lune Héliocentrique
    x1[4]    = positionLune[1];	 // Lune Héliocentrique
    x1[5]    = positionLune[2];		 //Lune Héliocentrique

    last = 0; // initialise le parametre d'ecriture
    vector< vector<double> >  matrix = {{},{},{},{},{},{},{},{}};
    matrix.reserve(nsteps);
        int compteur = 0;
    printOut(true,matrix); // ne pas ecrire premier pas de temps
  if (tol == 0){
    while(t < tfin -dt ) // boucle sur tout pas de temps
    {
      step(x,x1,dt);  // faire la mise a jour de la simulation
      t+=dt;

      printOut(false,matrix); // ecrire pas de temps actuel

      cout << t << "\r";
    }
    if (tfin > t){
      dt = tfin-t;
      step(x,x1 ,dt);
      t+=dt;
      printOut(false,matrix);
    }
    printOut(true,matrix); // ecrire dernier pas de temps
  }
  else {
    while(t < tfin -dt ) // boucle sur tout pas de temps
    {

      step2(dt);  // faire la mise a jour de la simulation
      t+=dt;

      printOut(false,matrix); // ecrire pas de temps actuel
    }
    if ( tfin > t){
    dt = tfin-t;
    step(x,x1,dt);
    t+=dt;
    printOut(false,matrix);}

    printOut(true,matrix); // ecrire dernier pas de temps

  }

  for(int i(0); i<matrix[0].size(); i++){
 *outputFile << matrix[0][i]  << " "<<matrix[1][i] << " "<< matrix[2][i] << " " <<matrix[3][i]   << " "<<matrix[4][i] << " "<< matrix[5][i] << " " <<matrix[6][i]   << " "<<matrix[7][i] << " "<<endl; // write output on file
}
};

};


// Extension de la class Engine implementant l'integrateur Runge-Kutta 4
class EngineRungeKutta4: public Engine
{
public:

  // construire la class Engine
  EngineRungeKutta4(ConfigFile configFile): Engine(configFile) {}



  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le scheme: Runge-kutta 4
  */
  void step(valarray<double>& x_,valarray<double>& x1_, double dt_) override
  {

    valarray<double> k1 =valarray<double>(6);
    valarray<double> k2 =valarray<double>(6);
    valarray<double> k3 =valarray<double>(6);
    valarray<double> k4 =valarray<double>(6);
    double sec = second; int min = minute; int hr = hour; int jr = day; int mois = month; int yr = year;
    k1 = fonction(x_,x1_,dt_,sec,min, hr,jr, mois, yr);
    k1 *= dt_;
    sec += dt_/2.0;
    chgmtemps(jr,mois,yr,hr,min,sec);
    k2 = fonction(x_+0.5*k1,x1_,dt_,sec,min, hr,jr, mois, yr);
    k2 *= dt_;
    k3 =fonction(x_+0.5*k2,x1_,dt_,sec,min, hr,jr, mois, yr);
    k3 *= dt_;
    sec += dt_/2.0;
    chgmtemps(jr,mois,yr,hr,min,sec);
    k4 = fonction(x_+k3,x1_,dt_,sec,min, hr,jr, mois, yr);
    k4 *= dt_;
    x_ +=(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
    second += dt_;
    chgmtemps(day,month,year,hour,minute,second);
    Soleil= Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
    Lune = Ephemeris::solarSystemObjectAtDateAndTime(EarthsMoon, day, month, year, hour, minute, second);
    //Mise a jour des positions des astres
    x1_[6]    = 0.0;		 // lire composante x position Terre
    x1_[7]    = 0.0;		 // lire composante y position
    x1_[8]    = 0.0;		 // lire composante z position
    x1_[0]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[0];		 // lire composante x position  Corps Soleil
    x1_[1]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[1];		 // lire composante y position  Corps Soleil
    x1_[2]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[2];		 // lire composante z position  Corps Soleil
    x1_[3]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[0];		 // lire composante x position Lune
    x1_[4]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[1];	 // lire composante y position Lune
    x1_[5]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[2];		 // lire composante z position Lune

//
/*
    // Tentative avec Euler-Chromer
    valarray<double> positiontemporaire = x_[slice(0,3,1)];
    valarray<double> vitessetemporaire = x_[slice(3,3,1)];
    valarray<double> d_terresoleil =valarray<double>(0e0,3);
    positiontemporaire += vitessetemporaire*dt_;
    vitessetemporaire += acceleration(x_,x1_,dt_,second+dt_,minute, hour,day, month, year)*dt_;


    second += dt_;
    //Mise à jour des astres
    chgmtemps(day,month,year,hour,minute,second);
    Soleil= Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
    Lune = Ephemeris::solarSystemObjectAtDateAndTime(EarthsMoon, day, month, year, hour, minute, second);
    // La nouvelle position, celui du satellite
    x_[slice(0,3,1)] = positiontemporaire;
    x_[slice(3,3,1)] = vitessetemporaire;
    //Mise a jour des positions des astres
    x1_[6]    = 0.0;		 // lire composante x position Terre
    x1_[7]    = 0.0;		 // lire composante y position
    x1_[8]    = 0.0;		 // lire composante z position
    x1_[0]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[0];		 // lire composante x position  Corps Soleil
    x1_[1]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[1];		 // lire composante y position  Corps Soleil
    x1_[2]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[2];		 // lire composante z position  Corps Soleil
    x1_[3]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[0];		 // lire composante x position Lune
    x1_[4]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[1];	 // lire composante y position Lune
    x1_[5]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[2];		 // lire composante z position Lune
*/

  }

  void step2(double& dt_) override // Methode utile pour le pas de temps adaptatif
  {
    valarray<double> x_ =valarray<double>(6);
    valarray<double> x__ =valarray<double>(6);
    x_=x;
    x__=x;
    double d(0.0);

    step(x_,x1,dt_);

    step(x__,x1,dt_*0.5);

    step(x__,x1,dt_*0.5);
    valarray<double> diff =valarray<double>(30);
    diff = x_-x__;
    diff = diff[slice(0,3,1)];
    d = norm2(diff);
      if(d > tol){
        while (d > tol && d >= 1.e-29) { // EVITER LA DIVISON PAR 0
          x_=x;
          x__=x;
          dt_ *= 0.99*pow((tol/d),1.0/5.0);
          step(x_,x1,dt_);
          step(x__,x1,dt_*0.5);
          step(x__,x1,dt_*0.5);
          diff = x_-x__;
          diff = diff[slice(0,3,1)];
          d = norm2(diff);
        }
        step(x,x1,dt_);
      }
      else if(d < tol && d >= 1.e-29){
         while (d <tol && d >= 1.e-29) {
            x_=x;
            x__=x;
            dt_ *= (2-0.99)*pow((tol/d),1.0/5.0);
            step(x_,x1,dt_);
            step(x__,x1,dt_*0.5);
            step(x__,x1,dt_*0.5);
            diff = x_-x__;
            diff = diff[slice(0,3,1)];
            d = norm2(diff);
          }
          step(x,x1,dt_);

      }
 dt = dt_;


  }

};
class EngineStormerVerlet : public Engine // Méthode de Stormer Verlet
{
public:
  EngineStormerVerlet(ConfigFile configFile): Engine(configFile) {}

  void step(valarray<double>& x_,valarray<double>& x1_, double dt_) override
  {

  valarray<double> temp = valarray<double>(0.0,6);
valarray<double> q = x_[slice(0,3,1)];
valarray<double> p = x_[slice(3,3,1)];
valarray<double> avant = x1_;
    temp[slice(0,3,1)] = dt_*p+dt_*dt_*acceleration(x_,x1_,dt_,second,minute, hour,day, month, year)/(2.0);
    valarray<double> positionSoleil =  -equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance);

    second+=dt_;
    chgmtemps(day,  month,  year,  hour, minute, second );
    Soleil= Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
    Lune = Ephemeris::solarSystemObjectAtDateAndTime(EarthsMoon, day, month, year, hour, minute, second);
    x1_[6]    = 0.0;		 // lire composante x position
    x1_[7]    = 0.0;		 // lire composante y position
    x1_[8]    = 0.0;		 // lire composante z position
    x1_[0]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[0];		 // lire composante x position  Corps
    x1_[1]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[1];		 // lire composante y position  Corps
    x1_[2]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[2];		 // lire composante z position  Corps
    x1_[3]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[0] ;		 // lire composante x position Lune
    x1_[4]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[1] ;	 // lire composante y position Lune
    x1_[5]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[2] ;
temp[slice(3,3,1)] +=dt_/(2.0)*(acceleration(x_+temp,x1_,dt_,second,minute, hour,day, month, year)+acceleration(x_,avant,dt_,second-dt,minute,  hour,day, month, year)); // erreur possible
valarray<double> d_terresoleil = valarray<double>(0e0,3);


x_ += temp;

  }
void  step2(double& dt_) override{ // Même fonction
	// DON'T WORK !!
  valarray<double> x_ =valarray<double>(3);
  valarray<double> x__ =valarray<double>(3);
  x_=x;
  x__=x;
  double d(0.0);

  step(x_,x1,dt_);

  step(x__,x1,dt_*0.5);

  step(x__,x1,dt_*0.5);
  valarray<double> diff =valarray<double>(30);
  diff = x_-x__;
  diff = diff[slice(0,3,1)];
  d = norm2(diff);
    if(d > tol){
      while (d > tol && d >= 1.e-29) {
        x_=x;
        x__=x;
        dt_ *= 0.99*pow((tol/d),1.0/5.0);
        step(x_,x1,dt_);
        step(x__,x1,dt_*0.5);
        step(x__,x1,dt_*0.5);
        diff = x_-x__;
        diff = diff[slice(0,3,1)];
        d = norm2(diff);
      }
      step(x,x1,dt_);
    }
    else if(d < tol && d >= 1.e-29){
       while (d <tol && d >= 1.e-29) {
          x_=x;
          x__=x;
          dt_ *= (2-0.99)*pow((tol/d),1.0/5.0);
          step(x_,x1,dt_);
          step(x__,x1,dt_*0.5);
          step(x__,x1,dt_*0.5);
          diff = x_-x__;
          diff = diff[slice(0,3,1)];
          d = norm2(diff);
        }
        step(x,x1,dt_);

    }
dt = dt_;

  }
};
// programme
int main(int argc, char* argv[])
{
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  string inputPath("configuration.dat"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Schema numerique
  string schema(configFile.get<string>("schema"));

  Engine* engine; // definer la class pour la simulation
  // choisir quel schema utiliser

  if(schema == "RungeKutta4" || schema == "RK4")
  {
    // initialiser une simulation avec schema runge-kutta 4
    engine = new EngineRungeKutta4(configFile);
  } else if(schema == "StormerVerlet" || schema == "SV") {
	  // initialiser une simulation avec schema Stormer-Verlet
    engine = new EngineStormerVerlet(configFile);
  }
  else
  {
    cerr << "Schema inconnu" << endl;
    return -1;
  }

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation
  // cout << "Fin de la simulation." << endl;
end = std::chrono::system_clock::now();
std::chrono::duration<double> elapsed_seconds = end - start;
cout << elapsed_seconds.count() << endl;
  return 0;
}
