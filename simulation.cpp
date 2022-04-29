#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include <iterator>
#include <vector>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs
                          // Fichier .tpp car inclut fonctions template
#include<stdio.h>
#include "ephemeris/Ephemeris.hpp"

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

  /*while(second >= 60){
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

  }*/
  minute += int(second) / 60;
  second = int(second) % 60 + (second - int(second));
  heure += minute / 60;
  minute = minute % 60;
  day += heure / 24;

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
double declinaison, double distance){ // Transformation du système equatorial à Cartesien par rapport au point VERNAL
valarray<double> array3 = valarray<double>(3);
double hr2rad = 2*3.1415926535897932384626433832795028841971/24.0; // Conversion d'angle en heure -> radians
double deg2rad = 3.1415926535897932384626433832795028841971/180.0; // Conversion d'angle en degré -> radians
  array3[0] = distance*cos(hr2rad*ascension)*cos(deg2rad*declinaison); // premier composante
  array3[1] = distance*sin(hr2rad*ascension)*cos(deg2rad*declinaison); // deuxieme composante
  array3[2] = distance*sin(deg2rad*declinaison); // troisieme composante
  return array3;
}

valarray<double> continuationRADEC(double dt_,int day, int month, int year, int heure,int minute, double second  ){
chgmtemps(day,month,year,heure,minute,second);
SolarSystemObject Soleil_avant= Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, heure, minute, second);
  second += 50.0;
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

private:
  // definition des constantes
  const double pi=3.1415926535897932384626433832795028841971e0;
  const double G=6.674e-11;
  const double masse_soleil =1.988435e30 ;
  const double masse_lune = 7.3459e22;
  const double masse_terre = 5.97e24;
  const double rayon_terre = 6371.009;
  const double celeritas =  299792458;
  // definition des variables
  double tfin=0.e0;      // Temps final
  unsigned int nsteps=1; // Nombre de pas de temps

  valarray<double> x0=valarray<double>(0.e0,6); // vecteur contenant la position et vitesse initiale du Satellite en
			 		        // utilisant la sequence index-valeur: 0-2 : position 3-5 vitesse

  unsigned int sampling=1; // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  void printOut(bool write)
  {
    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
    *outputFile << t  << " "<<x[0]<< " "<< x[1] << " " << x[2]  << " "<<x1[6]<< " "<< x1[7] << " " << x1[8] << " " <<dt<< " "<<endl; // write output on file
      last = 1;
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
double area; // Aire du satellite recevant la lumière du soleil
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
  valarray<double> force= valarray<double>(0.e0,3);
  
  force = -G*masse_soleil/(norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_)));
  force[0] = force[0] * (x_[0]-x1_[0]);
  force[1] = force[1] * (x_[1]-x1_[1]);
  force[2] = force[2] * (x_[2]-x1_[2]);

  return force;
}

valarray<double> ForceGravitationLune(valarray<double> const& x_,valarray<double> const& x1_) const {

  valarray<double> force= valarray<double>(0.e0,3);
  force = -G*masse_lune/(norm2(distance(2,x_,x1_))*norm2(distance(2,x_,x1_))*norm2(distance(2,x_,x1_)));
  force[0] = force[0] * (x_[0]-x1_[3]);
  force[1] = force[1] * (x_[1]-x1_[4]);
  force[2] = force[2] * (x_[2]-x1_[5]);

    return force;
}

valarray<double> ForceGravitationTerre(valarray<double> const& x_,valarray<double> const& x1_) const {

  valarray<double> force= valarray<double>(0.e0,3);
  force = G*masse_terre/(norm2(distance(3,x_,x1_))*norm2(distance(3,x_,x1_))*norm2  (distance(3,x_,x1_)));
  force[0] = force[0] * (x1_[6] - x_[0]);
  force[1] = force[1] * (x1_[7] - x_[1]);
  force[2] = force[2] * (x1_[8] - x_[2]);

return force;
}

valarray<double> ForceFrottement(valarray<double> const& x_,valarray<double> const& x1_) const {
	// Force de frottement atmosphérique pour des satellites
	// x_ : Position du satellite (taille 6, 3 pos + 3 vit), x1_ : Positions des astres (taille 12, 3*4 pos)
	/*valarray<double> e_v = valarray<double>(0.e0,3); // Vecteur vitesse normalisé du satellite
	Valarray<double omega_T = valarray<double(0.e0,3); // Vecteur vitesse angulaire terrestre. ATTENTION! Axe z le long de l'axe de rotation de la Terre (En fait c'est bien)
	omega_T[0] = 0;
	omega_T[1] = 0;
	omega_T[2] = 0.7292*10^(-4); // rad/s
	valarray<double> v_r = valarray<double>(0.e0,3); // Vecteur vitesse relative du satellite par rapport à celle de l'atmosphère
	v_r = x_[slice(3,3,1)] - vectorProduct(omega_T,x_[slice(0,3,1)]);
	force = -0.5*C_d*A_drag*rho(x_,x1_)*v_r^2*e_v;*/
	valarray<double> force= valarray<double>(0.e0,3);
	return force;
	/*TO DO :
	 * Implémenter rho, C_d, A_drag, e_v*/
}
double rho(valarray<double> const& x_, valarray<double> const& x1_,int const& n) const {
	// Implémentation du modèle Harris-Priester pour décrire la pression atmosphérique entre 100km et 1000km d'altitude.
	// Arguments: Position du soleil, Vecteur position, type d'orbite: n=2 pour des orbites proches de l'équateur et n=6 pour des orbites polaires
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
	valarray<double> h_ = coeff[slice(0,N,3)]; 
	valarray<double> rho_m = coeff[slice(1,N,3)];
	valarray<double> rho_M = coeff[slice(2,N,3)];
	size_t i = 0;
	double h = altitude();
	if(h <= 1000){
		while (h_[i+1] <= h){ // Erreur si le satellite va plus haut que 1000 km
			i++;
		}
	}else{return 0;}
	double H_m = (h_[i] - h_[i+1])/(log(rho_m[i+1]/rho_m[i]));
	double H_M = (h_[i] - h_[i+1])/(log(rho_M[i+1]/rho_M[i]));
	double rho_m_h = rho_m[i]*exp((h_[i]-h)/H_m);
	double rho_M_h = rho_M[i]*exp((h_[i]-h)/H_M);
	double cos_psi = 1; // A MODIFIER
	return rho_m_h + (rho_M_h - rho_m_h)*pow(cos_psi,n);
 }
valarray<double> ForceSolaire(valarray<double> const& x_,valarray<double> const& x1_) const {

valarray<double> force= valarray<double>(0.e0,3);
// Force de radiation :
// Cst Solar formula : sigma T^4 (R_s/R)^2
/*sigma = 5.670374e-8;T_s = 5778;R_s = 6.957e8;
solarCst = sigma*pow(5778,4)*pow((R_s/norm2(x_),2));
force = solarCst/celeritas * area ;
force[0] = x_[0]/norm2(x_);
force[1] = x_[1]/norm2(x_);
force[2] = x_[2]/norm2(x_);*/
return force;
}
valarray<double> ForceCoriolis(valarray<double> const& x_,valarray<double> const& x1_,double dt_, double second, int minute, int heure, int jour, int mois, int annee) const {

  valarray<double> vitesse_terre = valarray<double>(0.e0,3);
  vitesse_terre = (continuationRADEC(dt_, jour, mois, annee, heure, minute, second) - x1_[slice(0, 3, 1)])/ dt_;

  valarray<double> vitesser = x_[slice(3,3,1)];
  
  // -2m \omega x vitessesatellite
  double distance = sqrt(x1_[0] * x1_[0] + x1_[1] * x1_[1] + x1_[2] * x1_[2]);
  double vitesseangulaire = norm2(vitesse_terre)/distance;

  valarray<double> ez = valarray<double>(0.e0,3);
  ez[2] = vitesseangulaire;

  valarray<double> resultante = valarray<double>(0.e0, 3);
  resultante = vectorProduct(ez, vitesser);

  return -2*resultante;
}

valarray<double> ForceCentrifuge(valarray<double> const& x_,valarray<double> const& x1_,double dt_, double second, int minute, int heure, int jour, int mois, int annee) const {
  // OK
  // \omega cross \omega cross position

  valarray<double> vect_der = valarray<double>(0.e0,3);
  valarray<double> position = valarray<double>(0.e0,3);
  valarray<double> resultante = valarray<double>(0.e0,3);

  position = x_[slice(0,3,1)];
  vect_der = (continuationRADEC(dt_,jour, mois, annee, heure, minute,second+dt_) - x1_[slice(0,3,1)]) / dt_;

  // -2m \omega x vitessesatellite
  double distance = sqrt(x1_[0]*x1_[0] + x1_[1]*x1_[1]+x1_[2]*x1_[2]);
  double vitesseangulaire = norm2( vect_der)/distance;
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
  cout << "2e der " << norm2(vect_der) << endl; // problème est là, l'accélération angulaire est beaucoup trop grande !!

  vect_der = -vectorProduct(vect_der,position);

  return vect_der;
}


valarray<double> acceleration(valarray<double> const& x_,valarray<double> const& x1_,double dt_, double second, int minute, int heure, int jour, int mois, int annee) const{
  valarray<double> accelere(0.e1,3);
  /*accelere[0] = ForceGravitationSoleil(x_,x1_)[0]+ForceGravitationTerre(x_,x1_)[0]+ForceGravitationLune(x_,x1_)[0]+ForceFrottement(x_,x1_)[0]/mass+ForceSolaire(x_,x1_)[0]/mass ;
  accelere[1] = ForceGravitationSoleil(x_,x1_)[1]+ForceGravitationTerre(x_,x1_)[1]+ForceGravitationLune(x_,x1_)[1]+ForceFrottement(x_,x1_)[1]/mass+ForceSolaire(x_,x1_)[1]/mass;
  accelere[2] = ForceGravitationSoleil(x_,x1_)[2]+ForceGravitationTerre(x_,x1_)[2]+ForceGravitationLune(x_,x1_)[2]+ForceFrottement(x_,x1_)[2]/mass+ForceSolaire(x_,x1_)[2]/mass;
*/

accelere = ForceGravitationSoleil(x_,x1_)+ForceGravitationTerre(x_,x1_)+ForceGravitationLune(x_,x1_)+ForceFrottement(x_,x1_)/mass+ForceSolaire(x_,x1_)/mass +ForceCoriolis(x_,x1_,dt_,second,minute, heure,jour, mois, annee)+ForceCentrifuge(x_,x1_,dt_,second,minute, heure,jour, mois, annee)+ForceEuler(x_,x1_,dt_,second,minute, heure,jour, mois, annee);
  return accelere;
}

public:

  /*valarray<double> fonction(valarray<double>const& x_,valarray<double> const& x1_) const {
      valarray<double> f = valarray<double>(0.e0, 6);
      valarray<double> accelere = valarray<double>(0.e0, 3);
      accelere = acceleration(x_,x1_);
      for (int i(0); i<3; i++){
        f[i] = x_[i+3];
        f[i+3] = accelere[i];
  }
      return f;
    }*/
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
    area   = configFile.get<double>("Area");		 // lire composante de l'aire
    day = configFile.get<double>("day"); // Lire l'heure de la détection
    month = configFile.get<double>("month"); // Lire l'heure de la détection
    year = configFile.get<double>("year"); // Lire l'heure de la détection
    hour = configFile.get<double>("hour"); // Lire l'heure de la détection
    minute = configFile.get<double>("minute"); // Lire l'heure de la détection
    second = configFile.get<double>("second"); // Lire l'heure de la détection
    sampling = configFile.get<unsigned int>("sampling"); // lire le parametre de sampling
    tol = configFile.get<double>("tol");
    dt = tfin / nsteps;          // calculer le time step
    /*Soleil= Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
    Lune = Ephemeris::solarSystemObjectAtDateAndTime(EarthsMoon, day, month, year, hour, minute, second); // initialisation du soleil et de la Lune
*/

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales
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
    printOut(true); // ne pas ecrire premier pas de temps
  if (tol == 0){
    while(t < tfin -dt ) // boucle sur tout pas de temps
    {
      step(x,x1,dt);  // faire la mise a jour de la simulation
      t+=dt;
      printOut(false); // ecrire pas de temps actuel
      cout << t << "\r";
    }
    if (tfin > t){
      dt = tfin-t;
      step(x,x1 ,dt);
      t+=dt;
      printOut(false);
    }
    printOut(true); // ecrire dernier pas de temps
  }
  else {
    while(t < tfin -dt ) // boucle sur tout pas de temps
    {

      step2(dt);  // faire la mise a jour de la simulation
      t+=dt;

      printOut(false); // ecrire pas de temps actuel
    }
    if ( tfin > t){
    dt = tfin-t;
    step(x,x1,dt);
    t+=dt;
    printOut(false);}
    printOut(true); // ecrire dernier pas de temps

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
    /*valarray<double> k1 =valarray<double>(6);
    valarray<double> k2 =valarray<double>(6);
    valarray<double> k3 =valarray<double>(6);
    valarray<double> k4 =valarray<double>(6);
    k1 = fonction(x_,x1_);
    k1 *= dt_;
    k2 = fonction(x_+0.5*k1,x1_);
    k2 *= dt_;
    k3 =fonction(x_+0.5*k2,x1_);
    k3 *= dt_;
    k4 = fonction(x_+k3,x1_);
    k4 *= dt_;
    x_ +=(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
    second += dt_;
    chgmtemps(day,month,year,hour,minute,second);
    Soleil= Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
    Lune = Ephemeris::solarSystemObjectAtDateAndTime(EarthsMoon, day, month, year, hour, minute, second);
    //Mise a jour des positions des astres
    x1_[0]    = 0.0;		 // lire composante x position Soleil
    x1_[1]    = 0.0;		 // lire composante y position
    x1_[2]    = 0.0;		 // lire composante z position
    x1_[6]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[0];		 // lire composante x position  Corps Terre
    x1_[7]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[1];		 // lire composante y position  Corps Terre
    x1_[8]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[2];		 // lire composante z position  Corps Terre
    x1_[3]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[0] + x1_[6];		 // lire composante x position Lune
    x1_[4]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[1] + x1_[7];	 // lire composante y position Lune
    x1_[5]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[2] + x1_[8];		 // lire composante z position Lune
    */
//
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
    x1_[6]    = 0.0;		 // lire composante x position Soleil
    x1_[7]    = 0.0;		 // lire composante y position
    x1_[8]    = 0.0;		 // lire composante z position
    x1_[0]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[0];		 // lire composante x position  Corps Terre
    x1_[1]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[1];		 // lire composante y position  Corps Terre
    x1_[2]    = equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[2];		 // lire composante z position  Corps Terre
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
  cout << "Fin de la simulation." << endl;


  return 0;
}
