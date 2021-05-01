#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

//On va passer par une classe

class eqd {
public:
  float theta, dtheta, alpha, dalpha, M1, M2,  t, h, T, g, l1, l2, ktheta1, kalpha1, kdtheta1, kdalpha1, ktheta2, kalpha2, kdtheta2, kdalpha2, ktheta3, kalpha3, kdtheta3, kdalpha3, ktheta4, kalpha4, kdtheta4, kdalpha4 ;

  float F(float t, float theta, float dtheta, float alpha, float dalpha) {
    return (pow(dtheta,2)*M2*l1*cos((alpha-theta)*M_PI/180)*sin((alpha-theta)*M_PI/180)+pow(dalpha,2)*M2*l2*sin((alpha-theta)*M_PI/180)-(M1+M2)*g*sin(theta*M_PI/180)+M2*cos((alpha-theta)*M_PI/180)*g*sin(alpha*M_PI/180))/((M1+M2)*l1-M2*l2*pow(cos((alpha-theta)*M_PI/180),2));
  }
  
  float G(float t, float theta, float dtheta, float alpha, float dalpha) {
    return (-pow(dalpha,2)*M2*l2*cos((alpha-theta)*M_PI/180)*sin((alpha-theta)*M_PI/180)+(M1+M2)*(g*sin(theta*M_PI/180)*cos((alpha-theta)*M_PI/180)-l1*pow(dtheta,2)*sin((alpha-theta)*M_PI/180)-g*sin(alpha*M_PI/180)))/((M1+M2)*l2-M2*l2*pow(cos((alpha-theta)*M_PI/180),2));
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
  point C;
  
  //On déclare la classe utilisé
  eqd eqd1;
  
  //On déclare nos variables:
  
  //PAS
  eqd1.h=0.1;
  
  eqd1.g=9.81;
  //Longeurs des deux fils
  eqd1.l1=1;
  eqd1.l2=1;
  
  //Masse des deux masses 
  eqd1.M1=1;
  eqd1.M2=1;
  
  //Conditions initiales pour le premier pendule
  eqd1.theta=45;
  eqd1.alpha=45;
  
  //Conditions initiale pour le deuxième pendule
  eqd1.dtheta=0;
  eqd1.dalpha=0;
  
  //Conditions de temps
  eqd1.t=0;
  eqd1.T=10;
  
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
  
  eqd1.ktheta1=eqd1.dtheta*eqd1.h;
  eqd1.kalpha1=eqd1.dalpha*eqd1.h;
  eqd1.kdtheta1=eqd1.F(eqd1.t, eqd1.theta, eqd1.alpha, eqd1.dtheta, eqd1.dalpha)*eqd1.h;
  eqd1.kdalpha1=eqd1.G(eqd1.t, eqd1.theta, eqd1.alpha, eqd1.dtheta, eqd1.dalpha)*eqd1.h;

  eqd1.ktheta2=(eqd1.dtheta+(1/2)*eqd1.kdtheta1)*eqd1.h;
  eqd1.kalpha2=(eqd1.dalpha+(1/2)*eqd1.kdalpha1)*eqd1.h;
  eqd1.kdtheta2=eqd1.F(eqd1.t+ eqd1.h/2, eqd1.theta+eqd1.ktheta1/2, eqd1.alpha+eqd1.kalpha1/2, eqd1.dtheta+eqd1.kdtheta1/2, eqd1.dalpha+eqd1.kdalpha1/2)*eqd1.h;
  eqd1.kdtheta2=eqd1.G(eqd1.t+ eqd1.h/2, eqd1.theta+eqd1.ktheta1/2, eqd1.alpha+eqd1.kalpha1/2, eqd1.dtheta+eqd1.kdtheta1/2, eqd1.dalpha+eqd1.kdalpha1/2)*eqd1.h;

  eqd1.ktheta3=(eqd1.dtheta+(1/2)*eqd1.kdtheta2)*eqd1.h;
  eqd1.kalpha3=(eqd1.dalpha+(1/2)*eqd1.kdalpha2)*eqd1.h;
  eqd1.kdtheta3=eqd1.F(eqd1.t+ eqd1.h/2, eqd1.theta+eqd1.ktheta2/2, eqd1.alpha+eqd1.kalpha2/2, eqd1.dtheta+eqd1.kdtheta2/2, eqd1.dalpha+eqd1.kdalpha2/2)*eqd1.h;
  eqd1.kdtheta3=eqd1.G(eqd1.t+ eqd1.h/2, eqd1.theta+eqd1.ktheta2/2, eqd1.alpha+eqd1.kalpha2/2, eqd1.dtheta+eqd1.kdtheta2/2, eqd1.dalpha+eqd1.kdalpha2/2)*eqd1.h;

  eqd1.ktheta4=(eqd1.dtheta+eqd1.kdtheta3)*eqd1.h;
  eqd1.kalpha4=(eqd1.dalpha+eqd1.kdalpha3)*eqd1.h;
  eqd1.kdtheta4=eqd1.F(eqd1.t+ eqd1.h, eqd1.theta+eqd1.ktheta3, eqd1.alpha+eqd1.kalpha3, eqd1.dtheta+eqd1.kdtheta3, eqd1.dalpha+eqd1.kdalpha3)*eqd1.h;
  eqd1.kdtheta4=eqd1.G(eqd1.t+ eqd1.h, eqd1.theta+eqd1.ktheta3, eqd1.alpha+eqd1.kalpha3, eqd1.dtheta+eqd1.kdtheta3, eqd1.dalpha+eqd1.kdalpha3)*eqd1.h;

  eqd1.theta=eqd1.theta+(1/6)*(eqd1.ktheta1+2*eqd1.ktheta2+2*eqd1.ktheta3+eqd1.ktheta4);
  eqd1.alpha=eqd1.alpha+(1/6)*(eqd1.kalpha1+2*eqd1.kalpha2+2*eqd1.kalpha3+eqd1.kalpha4);
  eqd1.dtheta=eqd1.dtheta+(1/6)*(eqd1.kdtheta1+2*eqd1.kdtheta2+2*eqd1.kdtheta3+eqd1.kdtheta4);
  eqd1.dalpha=eqd1.dalpha+(1/6)*(eqd1.kdalpha1+2*eqd1.kdalpha2+2*eqd1.kdalpha3+eqd1.kdalpha4);

   
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
  
 return 0;
}
