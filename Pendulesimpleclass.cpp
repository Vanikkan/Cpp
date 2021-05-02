#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

//On va passer par une classe

class eqd {
public:
  float theta, dtheta, t, h, T, k1, k2, k3, k4, g, l;

  float f(float t, float theta, float dtheta) {
    return -(g/l)*sin(theta*M_PI/180);
  }
};
//On declare la classe point

class point {
  
 public:
  float x, y;

};

int main () {
  //On déclare nos points
  point A;
  point B; 
  //On déclare la classe utilisé
  eqd eqd1;

  //On déclare nos variables
  eqd1.g=9.81;
  eqd1.l=1;
  eqd1.theta=45;
  eqd1.dtheta=0;
  eqd1.T=20;
  eqd1.h=0.01;
  eqd1.k1;
  eqd1.k2;
  eqd1.k3;
  eqd1.k4;
  A.x=0;
  A.y=0;
  B.x=sin(eqd1.theta*M_PI/180);
  B.y=-cos(eqd1.theta*M_PI/180);
  ofstream fichier {"pendule_simple.out"};
  ofstream fichier1{"pendule_simple_point.out"};
  //Solveur 
while (eqd1.t<=eqd1.T)
  {
 eqd1.k1=eqd1.f(eqd1.t,eqd1.theta,eqd1.dtheta);
 eqd1.k2=eqd1.f(eqd1.t+(eqd1.h/2),eqd1.theta+(eqd1.h/2)*eqd1.dtheta,eqd1.dtheta+(eqd1.h*eqd1.k1/2));
 eqd1.k3=eqd1.f(eqd1.t+(eqd1.h/2),eqd1.theta+(eqd1.h/2)*eqd1.dtheta+(eqd1.h*eqd1.h/4)*eqd1.k1,eqd1.dtheta+(eqd1.h*eqd1.k2)/2);
 eqd1.k4=eqd1.f(eqd1.t+eqd1.h,eqd1.theta+eqd1.h*eqd1.dtheta+(eqd1.h*eqd1.h/2)*eqd1.k2,eqd1.dtheta+eqd1.h*eqd1.k3);
 eqd1.theta=eqd1.theta+eqd1.h*eqd1.dtheta+(eqd1.h*eqd1.h)*(eqd1.k1+eqd1.k2+eqd1.k3)/6;
 eqd1.dtheta=eqd1.dtheta+eqd1.h*(eqd1.k1+2*eqd1.k2+2*eqd1.k3+eqd1.k4)/6;
 eqd1.t=eqd1.t+eqd1.h;
 B.x=sin(eqd1.theta*M_PI/180);
 B.y=-cos(eqd1.theta*M_PI/180);
 fichier<<"  "<<eqd1.t<<"  "<<eqd1.theta<<"  "<<eqd1.dtheta<<endl;		  	  	   
 fichier1<<" "<<A.x<<" "<<A.y<<endl;
 fichier1<<" "<<B.x<<" "<<B.y<<endl;
 fichier1<<" "<<endl;
 fichier1<<" "<<endl; 
    
 }
ofstream animation("animationpendule.gnu"); // création du fichier de sortie pour le gif

    animation<<"set title \"Animation du pendule au cours du temps\""<<endl;
    animation<<"set terminal gif animate delay "<<eqd1.h*100<<endl;
    animation<<"set output \"penduleanime.gif\"" <<endl;
    animation<<"set xrange ["<<-(eqd1.l)<<":"<<(eqd1.l)<<"]"<<endl; //Permet d'adapter la fenetre à la longeur des fils en x
   
    animation<<"set yrange ["<<-(eqd1.l)<<":"<<(eqd1.l)<<"]"<<endl; //Permet d'adapter la fenetre à la longeur des fils en y
   
 
   
    animation<<"set tics nomirror"<<endl; //Impossible d'avoir une image qui apparait 2 fois
    animation<<"set pointintervalbox 3"<<endl; //Définis un style de ligne
    animation<<"do for [i=1:999] {plot \"pendule_simple_point.out\" index i using 1:2  title 'double-pendule' w linespoints ls 1}"<<endl;

    system("gnuplot animationpendule.gnu");

 
 return 0;
}
