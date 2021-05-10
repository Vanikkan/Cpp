#include <iostream>
#include <cmath>
#include <fstream>


using namespace std;

//On va commencer par déclarer la classe point, qui va nous permettre de générer notre pendule: 

class point{
  
public: //On commence par déclarer les variables qui materialiserons le point.

  //On peut également ici définir l'angle en degrès que fera le pendule avec la verticale:
  double Theta;
  
  //La dérivée de la fonction angle par rapport au temps
  double dTheta;

  //Les coordonnées cartesiennes de notre point

  double x,y; 

  //La distance du point par rapport à l'origine (la longueurs du pendule)
  double l;

  //La masse de la amsse asscoié à ce point;
  double m; 

  //On définis maintenant les constructeurs associé à cette classe point:
  
  point(); //Le constructeur: il va "créer" l'objet "point" quand celui ci sera appellé: il prépare la RAM a recevoir des instructions. 
  point( double, double, double, double ); // On déclare la fonction point associé au constructeur point.
  
  //Ici on aura deux pendule identiques. On peut donc utiliser un constructeur de copie, qui va nous permettre de copier le point déjà définis à la ligne d'au dessus.
  
  point ( point &);//Constructeur de copie
  
  void init(double, double, double, double); //On crée la fonction initialisation qui va initialiser nos variables

  //On déclare maintenant les variables au moyen de la commande "Get", ce sont nos accesseurs; ils vont nous permettre d'acceder aux données de notre classe privée. 
  //On leur attribut la variable qu'il devront retourner. La commande "const" permet de fixer ces variables, car on ne veut pas de variables dynamique ici. 

  double GetTheta()  {return Theta;}
  double GetdTheta() {return dTheta;}
  double Getl()  {return l;}
  double Getm()  {return m;}
  double Getx() {return l*sin(Theta*M_PI/180.0);}
  double Gety() {return -l*cos(Theta*M_PI/180.0);}

  //On passe maintenant aux mutateurs. Ce sont des fonctions membres qui permettent d'attribuer et de changer la valeur d'une variable protégée. Cela nous permet ici d'initialiser nos variables.  


  double Setx(double x0) {return x=x0;}
  double Sety(double y0) {return y=y0;}
  double SetTheta(double Theta0) {return Theta=Theta0;}
  double SetdTheta(double dTheta0) {return dTheta=dTheta0;}
  double Setl(double l0) {return l=l0;}
  double Setm(double m0) {return m=m0;}



  
};

//Maintenant, on va initialiser nos fonctions membres.

//On commence par initialiser la fonction de notre premier point.

  point::point(  double Theta0, double dTheta0, double l0, double m0) {
 
    Theta=Theta0;
    dTheta=dTheta0;
    l=l0;
    m=m0;
  };
 
//On définis maintenant notre constructeur de copie.

point::point( point &point2){

  Theta=point2.GetTheta();
  dTheta=point2.GetdTheta();
  l=point2.Getl();
  m=point2.Getl();
};

// Pour éviter tout problème, on va demander à notre programme de commencer par initialiser toutes les variables de la classe point à 0 par défaut.

point::point(){

  
  Theta=0;
  dTheta=0;
  l=0;
  m=0;
};

// On initialise notre classe point, maintenant que l'on part bien de 0.

void point::init( double Theta0, double dTheta0, double l0, double m0){

   
    Theta=Theta0;
    dTheta=dTheta0;
    l=l0;
    m=m0;
};

  //Equations differentiels double pendule 
double F( double g, double  M1, double M2, double L1, double L2, double theta, double dtheta, double alpha, double dalpha){
  
  return (dtheta*dtheta*M2*L1*cos((alpha-theta))*sin((alpha-theta))+dalpha*dalpha*M2*L2*sin((alpha-theta))-(M1+M2)*g*sin(theta)+M2*cos((alpha-theta))*g*sin(alpha))/((M1+M2)*L1-M2*L1*cos((alpha-theta)*M_PI/180)*cos((alpha-theta)*M_PI/180));
    
  };
  
double G( double g, double  M1, double M2, double L1, double L2, double theta, double dtheta, double alpha, double dalpha){
  
  return (-dalpha*dalpha*M2*L2*cos((alpha-theta))*sin((alpha-theta))+(M1+M2)*(g*sin(theta)*cos((alpha-theta))-L1*dtheta*dtheta*sin((alpha-theta))-g*sin(alpha)))/((M1+M2)*L2-M2*L2*cos((alpha-theta)*M_PI/180)*cos((alpha-theta)*M_PI/180));
    
  };

//On peut maintenant passer au main.

int main(){
  double g=9.81; 
  double t=0; //Le temps
  double h=0.001; //Le pas de RK4
  double T=10;  //Le temps maximum.  

  //On initialise les grandeurs utilisé dans RK4.

  double kx1=0, ky1=0, ku1=0, kv1=0;
  double kx2=0, ky2=0, ku2=0, kv2=0;
  double kx3=0, ky3=0, ku3=0, kv3=0;
  double kx4=0, ky4=0, ku4=0, kv4=0;
  
//On définis ici les grandeurs qui vont caractériser nos deux point;

  //Le premier point A, symbolisant la première masse.
  double Theta1, dTheta1, l1, m1;

  //Le deuxième point B, symbolisant la deuxième masse.

  double Theta2, dTheta2, l2, m2; 

  //On peut pour plus de modularité inscrire ces grandeurs dans un fichier de donnée.

  ifstream fichier("cond_init.txt"); 

  //On en profite pour définir les deux fichiers de sortie (Pour les angles, et les coordonnées).

  ofstream fichier1("pendule_double.out");
  ofstream fichier2("pendule_double_point.out");

  //Grace à notre fichier d'initialisation, on peut attribuer des valeurs à nos différentes variables caractérisant nos points. 
  fichier>>Theta1>>dTheta1>>l1>>m1;
  fichier>>Theta2>>dTheta2>>l2>>m2;
 
  //Passons maintenant à la création de notre double pendule. On commence par définir un point origine, qui sers de point d'ancrage au premier pendule.
  double xo=0;
  double yo=0;
  
  //On créer ensuite nos deux points grace à la classe point.
  
  point A(Theta1, dTheta1, l1, m1); 
  point B(Theta2, dTheta2, l2, m2); 

  //On initialise ces deux points.

  A.init(Theta1, dTheta1, l1, m1);
  B.init(Theta2, dTheta2, l2, m2);
  
  //On affiche ces grandeurs pour voir si l'initialisation à partir du fichier d'entrée s'est bien faite.

  cout<<"Theta1= "<<Theta1<<"dTheta1= "<<"l1= "<<l1<<"m1= "<<m1<<endl; 
  cout<<"Theta1= "<<Theta1<<"dTheta1= "<<"l1= "<<l1<<"m1= "<<m1<<endl; 

  //On passe maintenant à la méthode de résolution RK4 pour l'ordre 2 pour les équa diff couplé.
  A.Theta=A.GetTheta()*M_PI/180.0;
  B.Theta=B.GetTheta()*M_PI/180.0;
  
  while (t<=T){

  kx1=A.GetdTheta()*h;
  ky1=B.GetdTheta()*h;
  ku1=F( g, A.Getm(), B.Getm(), A.Getl() , B.Getl(), A.GetTheta(), A.GetdTheta(), B.GetTheta(), B.GetdTheta())*h;
  kv1=G( g, A.Getm(), B.Getm(), A.Getl() , B.Getl(), A.GetTheta(), A.GetdTheta(), B.GetTheta(), B.GetdTheta())*h;
  
  kx2=(A.GetdTheta()+(1/2.0)*ku1)*h;
  ky2=(B.GetdTheta()+(1/2.0)*kv1)*h;
  ku2=F( g, A.Getm(), B.Getm(), A.Getl() , B.Getl(), A.GetTheta()+kx1/2.0, A.GetdTheta()+ky1/2.0, B.GetTheta()+ku1/2.0, B.GetdTheta()+kv1/2.0)*h;
  kv2=G( g, A.Getm(), B.Getm(), A.Getl() , B.Getl(), A.GetTheta()+kx1/2.0, A.GetdTheta()+ky1/2.0, B.GetTheta()+ku1/2.0, B.GetdTheta()+kv1/2.0)*h;
 
  kx3=(A.GetdTheta()+(1/2.0)*ku2)*h;
  ky3=(B.GetdTheta()+(1/2.0)*kv2)*h;
  ku3=F( g, A.Getm(), B.Getm(), A.Getl() , B.Getl(), A.GetTheta()+kx2/2.0, A.GetdTheta()+ky2/2.0, B.GetTheta()+ku2/2.0, B.GetdTheta()+kv2/2.0)*h;
  kv3=G( g, A.Getm(), B.Getm(), A.Getl() , B.Getl(), A.GetTheta()+kx2/2.0, A.GetdTheta()+ky2/2.0, B.GetTheta()+ku2/2.0, B.GetdTheta()+kv2/2.0)*h;
 
  
  kx4=(A.GetdTheta()+ku3)*h;
  ky4=(B.GetdTheta()+kv3)*h;
  ku4=F( g, A.Getm(), B.Getm(), A.Getl() , B.Getl(), A.GetTheta()+kx3, A.GetdTheta()+ky3, B.GetTheta()+ku3, B.GetdTheta()+kv3)*h;
  kv4=G( g, A.Getm(), B.Getm(), A.Getl() , B.Getl(), A.GetTheta()+kx3, A.GetdTheta()+ky3, B.GetTheta()+ku3, B.GetdTheta()+kv3)*h;

  A.Theta+=(kx1+2*kx2+2*kx3+kx4)/6.0;
  B.Theta+=(ky1+2*ky2+2*ky3+ky4)/6.0;
  A.dTheta+=(ku1+2*ku2+2*ku3+ku4)/6.0;
  B.dTheta+=(kv1+2*kv2+2*kv3+kv4)/6.0;

  fichier1<<" "<<t<<" "<<A.Theta<<" "<<A.dTheta<<endl; 
  fichier1<<" "<<t<<" "<<B.Theta<<" "<<B.dTheta<<endl; 
  fichier1<<" "<<endl;
  fichier1<<" "<<endl; 
  t+=h; 

  A.x=A.Getx();
  A.y=A.Gety();
  B.x=A.Getx()+B.Getx();
  B.y=A.Gety()+B.Gety();

  fichier2<<" "<<t<<" "<<xo<<" "<<yo<<endl;
  fichier2<<" "<<t<<" "<<A.x<<" "<<A.y<<endl; 
  fichier2<<" "<<t<<" "<<B.x<<" "<<B.y<<endl; 
  fichier2<<" "<<endl;
  fichier2<<" "<<endl; 

  };
  // Création du gif
   ofstream animation("animationpendule.gnu"); // création du fichier de sortie pour le gif

    animation<<"set title \"Animation du pendule au cours du temps\""<<endl;
    animation<<"set terminal gif animate delay "<<h/100<<endl;
    animation<<"set output \"penduleanime.gif\"" <<endl;
    animation<<"set xrange ["<<-(A.Getl()+B.Getl())<<":"<<(A.Getl()+B.Getl())<<"]"<<endl; //Permet d'adapter la fenetre à la longeur des fils en x
   
    animation<<"set yrange ["<<-(A.Getl()+B.Getl())<<":"<<(A.Getl()+B.Getl())<<"]"<<endl; //Permet d'adapter la fenetre à la longeur des fils en y
   
 
   
    animation<<"set tics nomirror"<<endl; //Impossible d'avoir une image qui apparait 2 fois
    animation<<"set pointintervalbox 3"<<endl; //Définis un style de ligne
    animation<<"do for [i=1:999] {plot \"pendule_double_point.out\" index i using 2:3  title 'double-pendule' w linespoints ls 1}"<<endl;

    system("gnuplot animationpendule.gnu");
    
    




  

}
