#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs
                          // Fichier .tpp car inclut fonctions template
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

// Convertit le surplus de seconde (i.e. plus de 60 seconde) en minute, heure, jour, mois, année
void chgmtemps(int& day, int& month, int& year, int& heure,int& minute, int& second )
{
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
  if ((month == 1 or month ==3 or month == 5 or month ==7 or month ==8 or month == 10) && day > 31) {
    day -=31;
    month += 1;
  }
  else if ((month == 4 or month == 6 or month ==9 or month == 11) && day > 30) {
	day -=30;
	month += 1;
  }
  else if (month == 2 && day > 28 && year % 4 != 0) {
	day -=28;
	month += 1;
  }
  // Année bissextile
  else if (month == 2 && day > 28 && year % 4 == 0){
	day -=29;
	month += 1;
  }
  else if (month == 12 && day > 31){
	day -= 31;
	month -= 11;
	year += 1;
  }
}
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
double hr2rad = 2*3.1415926535897932384626433832795028841971/24; // Conversion d'angle en heure -> radians
double deg2rad = 3.1415926535897932384626433832795028841971/180; // Conversion d'angle en degré -> radians
  array3[0] = distance*cos(hr2rad*ascension)*cos(deg2rad*declinaison); // premier composante
  array3[1] = distance*sin(hr2rad*ascension)*cos(deg2rad*declinaison); // deuxieme composante
  array3[2] = distance*sin(deg2rad*declinaison); // troisieme composante
  return array3;
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
int day,month,year,hour,minute,second;
const double astronomical_unit= 1.496e11;
valarray<double> x1=valarray<double>(0.e0,12); // Position des astres

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
valarray<double> ForceGravitationSoleil(valarray<double> const& x_,valarray<double> const& x1_) const {
  valarray<double> force= valarray<double>(0.e0,3);

/*  force = G*masse_soleil/(norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_))*norm2(distance(1,x_,x1_)));
  force[0] = force[0] * (x_[0]-x1_[0]);
  force[1] = force[1] * (x_[1]-x1_[1]);
  force[2] = force[2] * (x_[2]-x1_[2]);*/
  return force;
}
valarray<double> ForceGravitationLune(valarray<double> const& x_,valarray<double> const& x1_) const {

  valarray<double> force= valarray<double>(0.e0,3);
  /*force = G*masse_lune/(norm2(distance(2,x_,x1_))*norm2(distance(2,x_,x1_))*norm2(distance(2,x_,x1_)));
  force[0] = force[0] * (x_[0]-x1_[3]);
  force[1] = force[1] * (x_[1]-x1_[4]);
  force[2] = force[2] * (x_[2]-x1_[5]);*/
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
valarray<double> ForceFrottement(valarray<double> const& x_,valarray<double> const& x1_) const {

valarray<double> force= valarray<double>(0.e0,3);
return force;
}
valarray<double> ForceSolaire(valarray<double> const& x_,valarray<double> const& x1_) const {

valarray<double> force= valarray<double>(0.e0,3);
return force;
}
valarray<double> acceleration(valarray<double> const& x_,valarray<double> const& x1_) const{
  valarray<double> accelere(0.e1,3);
  accelere[0] = ForceGravitationSoleil(x_,x1_)[0]+ForceGravitationTerre(x_,x1_)[0]+ForceGravitationLune(x_,x1_)[0]+ForceFrottement(x_,x1_)[0]+ForceSolaire(x_,x1_)[0];
  accelere[1] = ForceGravitationSoleil(x_,x1_)[1]+ForceGravitationTerre(x_,x1_)[1]+ForceGravitationLune(x_,x1_)[1]+ForceFrottement(x_,x1_)[1]+ForceSolaire(x_,x1_)[1];
  accelere[2] = ForceGravitationSoleil(x_,x1_)[2]+ForceGravitationTerre(x_,x1_)[2]+ForceGravitationLune(x_,x1_)[2]+ForceFrottement(x_,x1_)[2]+ForceSolaire(x_,x1_)[2];

  return accelere;
}
  // donnes internes
  double t,dt,tol;        // Temps courant pas de temps
  valarray<double> x=valarray<double>(0.e0,6); // Position des astres

public:

  valarray<double> fonction(valarray<double>const& x_,valarray<double> const& x1_) const {
      valarray<double> f = valarray<double>(0.e0, 6);
      valarray<double> accelere = valarray<double>(0.e0, 3);
      accelere = acceleration(x_,x1_);
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
   valarray<double> positionSoleil =  -equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance);
  Lune = Ephemeris::solarSystemObjectAtDateAndTime(EarthsMoon, day, month, year, hour, minute, second); // initialisation du soleil et de la Lune
  valarray<double> positionLune = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance); // Géocentrique
  x = x0; // Initialisation du satellite héliocentrique
  x[slice(0,3,1)] += positionSoleil;
    // Initialisation des positions
    x1[6]    = positionSoleil[0];		 // lire composante x position initiale Terre
    x1[7]    = positionSoleil[1];		 // lire composante y position initiale Terre
    x1[8]    = positionSoleil[2];		 // lire composante z position initiale Terre
    x1[0]    = 0.0;		 // lire composante x position initiale Corps Soleil
    x1[1]    = 0.0;		 // lire composante y position initiale Corps Soleil
    x1[2]    = 0.0;		 // lire composante z position initiale Corps Soleil
    x1[3]    = positionLune[0] + positionSoleil[0];		 //  Lune Héliocentrique
    x1[4]    = positionLune[1] + positionSoleil[1];	 // Lune Héliocentrique
    x1[5]    = positionLune[2] + positionSoleil[2];		 //Lune Héliocentrique

    last = 0; // initialise le parametre d'ecriture
    printOut(true); // ne pas ecrire premier pas de temps
  if (tol == 0){
    while(t < tfin -dt ) // boucle sur tout pas de temps
    {
      step(x,x1,dt);  // faire la mise a jour de la simulation
      t+=dt;
      printOut(false); // ecrire pas de temps actuel
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
    valarray<double> positionSoleil =  -equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance);
    valarray<double> positiontemporaire = x_[slice(0,3,1)];
    valarray<double> vitessetemporaire = x_[slice(3,3,1)];
    valarray<double> d_terresoleil =valarray<double>(0e0,3);

     positiontemporaire += vitessetemporaire*dt_;
    vitessetemporaire += acceleration(positiontemporaire,x1_)*dt_;



    positiontemporaire = positiontemporaire - positionSoleil; // Passage en géocentrique
    second += dt_;
    //Mise à jour des astres
    chgmtemps(day,month,year,hour,minute,second);
    Soleil= Ephemeris::solarSystemObjectAtDateAndTime(Sun, day, month, year, hour, minute, second);
    Lune = Ephemeris::solarSystemObjectAtDateAndTime(EarthsMoon, day, month, year, hour, minute, second);
    //Déplacement de la terre
     d_terresoleil = (-equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)-positionSoleil);
    // La nouvelle position est donc le déplacement de la terre + celui du satellite
    x_[slice(0,3,1)] = positiontemporaire+d_terresoleil;
    x_[slice(3,3,1)] = vitessetemporaire;
    //Mise a jour des positions des astres
    x1_[0]    = 0.0;		 // lire composante x position Soleil
    x1_[1]    = 0.0;		 // lire composante y position
    x1_[2]    = 0.0;		 // lire composante z position
    x1_[6]    = -equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[0];		 // lire composante x position  Corps Terre
    x1_[7]    = -equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[1];		 // lire composante y position  Corps Terre
    x1_[8]    = -equatorialtocartesian(Soleil.equaCoordinates.ra, Soleil.equaCoordinates.dec, astronomical_unit*Soleil.distance)[2];		 // lire composante z position  Corps Terre
    x1_[3]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[0] + x1_[6];		 // lire composante x position Lune
    x1_[4]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[1] + x1_[7];	 // lire composante y position Lune
    x1_[5]    = equatorialtocartesian(Lune.equaCoordinates.ra, Lune.equaCoordinates.dec, astronomical_unit*Lune.distance)[2] + x1_[8];		 // lire composante z position Lune
    x_[slice(0,3,1)] += positionSoleil; //Retour en héliocentrique

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

  valarray<double> temp = valarray<double>(5);
valarray<double> q = x_[slice(0,3,1)];
valarray<double> p = x_[slice(3,3,1)];
    temp[slice(0,3,1)] = dt*p+dt_*dt_*acceleration(x_,x1_)/(2.0);
temp[slice(3,3,1)] +=dt_/(2.0)*(acceleration(x_+temp,x1_)+acceleration(x_,x1_));

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
