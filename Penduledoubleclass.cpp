#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

//On va passer par une classe

class eqd {
public:
  double theta, dtheta, alpha, dalpha, M1, M2,  t, h, T, g, l1, l2, ktheta1, kalpha1, kdtheta1, kdalpha1, ktheta2, kalpha2, kdtheta2, kdalpha2, ktheta3, kalpha3, kdtheta3, kdalpha3, ktheta4, kalpha4, kdtheta4, kdalpha4 ;

  double F(double t, double dtheta, double dalpha, double theta, double alpha) {
    return (pow(dtheta,2)*M2*l1*cos((alpha-theta)*M_PI/180)*sin((alpha-theta)*M_PI/180)+pow(dalpha,2)*M2*l2*sin((alpha-theta)*M_PI/180)-(M1+M2)*g*sin(theta*M_PI/180)+M2*cos((alpha-theta)*M_PI/180)*g*sin(alpha*M_PI/180))/((M1+M2)*l1-M2*l1*pow(cos((alpha-theta)*M_PI/180),2));
  }
  
  double G(double t, double dtheta, double dalpha, double theta, double alpha) {
    return (-pow(dalpha,2)*M2*l2*cos((alpha-theta)*M_PI/180)*sin((alpha-theta)*M_PI/180)+(M1+M2)*(g*sin(theta*M_PI/180)*cos((alpha-theta)*M_PI/180)-l1*pow(dtheta,2)*sin((alpha-theta)*M_PI/180)-g*sin(alpha*M_PI/180)))/((M1+M2)*l2-M2*l2*pow(cos((alpha-theta)*M_PI/180),2));
  }
};
//On declare la classe point

class point {
  
 public:
  double x, y;

};

int main () {
  //On déclare nos points
  point A;
  point B;
  point C;
  
  //On déclare la classe utilisé
  eqd eqd1;
  
  //On déclare nos variables:

  //PAS
  eqd1.h=0.1;
  
  eqd1.g=9.81;
  //Longeurs des deux fils
  eqd1.l1=1;
  eqd1.l2=2;
  
  //Masse des deux masses 
  eqd1.M1=25;
  eqd1.M2=25;
  
  //Conditions initiales pour le premier pendule
  eqd1.theta=45;
  eqd1.alpha=45;
  
  //Conditions initiale pour le deuxième pendule
  eqd1.dtheta=0;
  eqd1.dalpha=0;
  
  //Conditions de temps
  eqd1.t=0;
  eqd1.T=100;
  
  //Déclaration variable RK4
  eqd1.ktheta1;
  eqd1.kalpha1;
  eqd1.kdtheta1;
  eqd1.kdalpha1;

  eqd1.ktheta2;
  eqd1.kalpha2;
  eqd1.kdtheta2;
  eqd1.kdalpha2;

  eqd1.ktheta3;
  eqd1.kalpha3;
  eqd1.kdtheta3;
  eqd1.kdalpha3;

  eqd1.ktheta4;
  eqd1.kalpha4;
  eqd1.kdtheta4;
  eqd1.kdalpha4;
  
  A.x=0;
  A.y=0;
  
  B.x=eqd1.l1*sin(eqd1.theta*M_PI/180);
  B.y=-eqd1.l1*cos(eqd1.theta*M_PI/180);

  C.x=eqd1.l1*sin(eqd1.theta*M_PI/180)+eqd1.l2*sin(eqd1.alpha*M_PI/180);
  C.y=-eqd1.l1*cos(eqd1.theta*M_PI/180)-eqd1.l2*cos(eqd1.alpha*M_PI/180);
  
  ofstream fichier {"pendule_double.out"};
  ofstream fichier1{"pendule_double_point.out"};
  //Solveur
  
  while (eqd1.t<=eqd1.T)
    {
      /*
  eqd1.ktheta1=eqd1.dtheta*eqd1.h;
  eqd1.kalpha1=eqd1.dalpha*eqd1.h;
  eqd1.kdtheta1=eqd1.F(eqd1.t, eqd1.theta, eqd1.alpha, eqd1.dtheta, eqd1.dalpha)*eqd1.h;
  eqd1.kdalpha1=eqd1.G(eqd1.t, eqd1.theta, eqd1.alpha, eqd1.dtheta, eqd1.dalpha)*eqd1.h;

  
  eqd1.ktheta2=(eqd1.dtheta+eqd1.kdtheta1/2.)*eqd1.h;
  eqd1.kalpha2=(eqd1.dalpha+eqd1.kdalpha1/2.)*eqd1.h;
  eqd1.kdtheta2=eqd1.F(eqd1.t+ eqd1.h*0.5, eqd1.theta+eqd1.ktheta1*0.5, eqd1.alpha+eqd1.kalpha1*0.5, eqd1.dtheta+eqd1.kdtheta1*0.5, eqd1.dalpha+eqd1.kdalpha1*0.5)*eqd1.h;
  eqd1.kdtheta2=eqd1.G(eqd1.t+ eqd1.h*0.5, eqd1.theta+eqd1.ktheta1*0.5, eqd1.alpha+eqd1.kalpha1*0.5, eqd1.dtheta+eqd1.kdtheta1*0.5, eqd1.dalpha+eqd1.kdalpha1*0.5)*eqd1.h;
 
  
  eqd1.ktheta3=(eqd1.dtheta+eqd1.kdtheta2*0.5)*eqd1.h;
  eqd1.kalpha3=(eqd1.dalpha+eqd1.kdalpha2*0.5)*eqd1.h;
  eqd1.kdtheta3=eqd1.F(eqd1.t+ eqd1.h*0.5, eqd1.theta+eqd1.ktheta2*0.5, eqd1.alpha+eqd1.kalpha2*0.5, eqd1.dtheta+eqd1.kdtheta2*0.5, eqd1.dalpha+eqd1.kdalpha2*0.5)*eqd1.h;
  eqd1.kdtheta3=eqd1.G(eqd1.t+ eqd1.h*0.5, eqd1.theta+eqd1.ktheta2*0.5, eqd1.alpha+eqd1.kalpha2*0.5, eqd1.dtheta+eqd1.kdtheta2*0.5, eqd1.dalpha+eqd1.kdalpha2*0.5)*eqd1.h;
 
  
  eqd1.ktheta4=(eqd1.dtheta+eqd1.kdtheta3)*eqd1.h;
  eqd1.kalpha4=(eqd1.dalpha+eqd1.kdalpha3)*eqd1.h;
  eqd1.kdtheta4=eqd1.F(eqd1.t+ eqd1.h, eqd1.theta+eqd1.ktheta3, eqd1.alpha+eqd1.kalpha3, eqd1.dtheta+eqd1.kdtheta3, eqd1.dalpha+eqd1.kdalpha3)*eqd1.h;
  eqd1.kdtheta4=eqd1.G(eqd1.t+ eqd1.h, eqd1.theta+eqd1.ktheta3, eqd1.alpha+eqd1.kalpha3, eqd1.dtheta+eqd1.kdtheta3, eqd1.dalpha+eqd1.kdalpha3)*eqd1.h;
 
  
  eqd1.theta=eqd1.theta+(eqd1.ktheta1+2*eqd1.ktheta2+2*eqd1.ktheta3+eqd1.ktheta4)*0.166;
  eqd1.alpha=eqd1.alpha+(eqd1.kalpha1+2*eqd1.kalpha2+2*eqd1.kalpha3+eqd1.kalpha4)*0.166;
  eqd1.dtheta=eqd1.dtheta+(eqd1.kdtheta1+2*eqd1.kdtheta2+2*eqd1.kdtheta3+eqd1.kdtheta4)*0.166;
  eqd1.dalpha=eqd1.dalpha+(eqd1.kdalpha1+2*eqd1.kdalpha2+2*eqd1.kdalpha3+eqd1.kdalpha4)*0.166;
      */
       double k1=eqd1.F(eqd1.t,eqd1.dtheta,eqd1.dalpha,eqd1.theta,eqd1.alpha);
    double j1=eqd1.G(eqd1.t,eqd1.dtheta,eqd1.dalpha,eqd1.theta,eqd1.alpha);
    double b1=eqd1.dtheta;
    double p1=eqd1.dalpha;
    
    //2eme coefficient
    double k2=eqd1.F(eqd1.t,eqd1.dtheta+eqd1.h/2.*k1,eqd1.dalpha+eqd1.h/2.*j1,eqd1.theta+eqd1.h/2*b1,eqd1.alpha+eqd1.h/2.*p1);
    double j2=eqd1.G(eqd1.t,eqd1.dtheta+eqd1.h/2.*k1,eqd1.dalpha+eqd1.h/2.*j1,eqd1.theta+eqd1.h/2*b1,eqd1.alpha+eqd1.h/2.*p1);
    double b2=eqd1.dtheta+eqd1.h/2.*b1;
    double p2=eqd1.dalpha+eqd1.h/2.*p1;
    
    //3eme coefficient
    double k3=eqd1.F(eqd1.t,eqd1.dtheta+eqd1.h/2.*k2,eqd1.dalpha+eqd1.h/2.*j2,eqd1.theta+eqd1.h/2.*b2,eqd1.alpha+eqd1.h/2.*p2);
    double j3=eqd1.G(eqd1.t,eqd1.dtheta+eqd1.h/2.*k2,eqd1.dalpha+eqd1.h/2.*j2,eqd1.theta+eqd1.h/2.*b2,eqd1.alpha+eqd1.h/2.*p2);
    double b3=eqd1.dtheta+eqd1.h/2.*b2;
    double p3=eqd1.dalpha+eqd1.h/2.*p2;
    
    //4eme coefficient
    double k4=eqd1.F(eqd1.t,eqd1.dtheta+eqd1.h*k3,eqd1.dalpha+eqd1.h*j3,eqd1.theta+eqd1.h*b3,eqd1.alpha+eqd1.h*p3);
    double j4=eqd1.G(eqd1.t,eqd1.dtheta+eqd1.h*k3,eqd1.dalpha+eqd1.h*j3,eqd1.theta+eqd1.h*b3,eqd1.alpha+eqd1.h*p3);
    double b4=eqd1.dtheta+eqd1.h*b3;
    double p4=eqd1.dalpha+eqd1.h*p3;
    
    //valeur
    eqd1.dtheta=eqd1.dtheta+eqd1.h/6.*(k1+2*k2+2*k3+k4);
    eqd1.dalpha=eqd1.dalpha+eqd1.h/6.*(j1+2*j2+2*j3+j4);
    eqd1.theta = eqd1.theta+eqd1.h/6.*(b1+2*b2+2*b3+b4);
    eqd1.alpha =eqd1.alpha+eqd1.h/6.*(p1+2*p2+2*p3+p4);
    
  B.x=eqd1.l1*sin(eqd1.theta*M_PI/180);
  B.y=-eqd1.l1*cos(eqd1.theta*M_PI/180);

  C.x=eqd1.l1*sin(eqd1.theta*M_PI/180)+eqd1.l2*sin(eqd1.alpha*M_PI/180);
  C.y=-eqd1.l1*cos(eqd1.theta*M_PI/180)-eqd1.l2*cos(eqd1.alpha*M_PI/180);

  fichier<<"  "<<eqd1.t<<"  "<<eqd1.theta<<"  "<<eqd1.alpha<<" "<<eqd1.dtheta<<" "<<eqd1.dalpha<<endl;

  fichier1<<" "<<A.x<<" "<<A.y<<endl;
  fichier1<<" "<<B.x<<" "<<B.y<<endl;
  fichier1<<" "<<C.x<<" "<<C.y<<endl;
  fichier1<<" "<<endl;
  fichier1<<" "<<endl;
  
  cout<<" "<<eqd1.t<<endl;

  eqd1.t=eqd1.t+eqd1.h;
  
 }
  
    ofstream animation("animationpendule.gnu"); // création du fichier de sortie pour le gif

    animation<<"set title \"Animation du pendule au cours du temps\""<<endl;
    animation<<"set terminal gif animate delay "<<eqd1.h*20<<endl;
    animation<<"set output \"penduleanime.gif\"" <<endl;
    animation<<"set xrange ["<<-(eqd1.l1+eqd1.l2)<<":"<<(eqd1.l1+eqd1.l2)<<"]"<<endl; //Permet d'adapter la fenetre à la longeur des fils en x
   
    animation<<"set yrange ["<<-(eqd1.l1+eqd1.l2)<<":"<<(eqd1.l1+eqd1.l2)<<"]"<<endl; //Permet d'adapter la fenetre à la longeur des fils en y
   
 
   
    animation<<"set tics nomirror"<<endl; //Impossible d'avoir une image qui apparait 2 fois
    animation<<"set pointintervalbox 3"<<endl; //Définis un style de ligne
    animation<<"do for [i=1:999] {plot \"pendule_double_point.out\" index i using 1:2  title 'double-pendule' w linespoints ls 1}"<<endl;

    system("gnuplot animationpendule.gnu");
    
    
 return 0;
}
