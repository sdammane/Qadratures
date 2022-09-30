#include <iostream>
#include "Lobatto.hpp"
#include <cmath>

//Le constructeur qui initialise une quadrature à n points
QuadratureLobatto::QuadratureLobatto(int n3){
  n2=n3;
}

//La fonction de changement de variable pour se ramener à
//l'intervlle [-1,1]
double QuadratureLobatto::g1(double t){

    return f1( ( (b2-a2)*t+(b2+a2) ) /2.);
}


//Surchage de l'operateur ()
double QuadratureLobatto:: operator() (double (*g2) (double),double a3,double b3){
  //Initialisation des bornes de l'intervalle
  a2=a3;
  b2=b3;
  //Initialiation de la fonction
  f1=g2;

  return Gauss_Lobatto() ;
}





double QuadratureLobatto::Gauss_Lobatto(){

  //Les poids associés aux points xi
  double w1,w2,w3,w4,w5;
  double x1,x2,x3,x4,x5;

  //Variable pour la quadrature de Gauss Legendre
  double quadrature_Lobatto;



  if(n2==2){
     w1=1;
     w2=1;

     x1= -1;
     x2= 1;


    quadrature_Lobatto= ( (b2-a2) / 2. ) * (w1 * g1(x1) + w2 * g1(x2) );
  }

  else if(n2==3){
     w1= 1./3.;
     w2= 4./3.;
     w3= 1./3.;

     x1= -1;
     x2=0;
     x3= 1;


    quadrature_Lobatto= ( (b2-a2) / 2. ) * (w1 * g1(x1) + w2 * g1(x2) + w3 * g1(x3) );
  }

  else if(n2==4){
     w1= 1. / 6.;
     w2= 5. / 6.;
     w3= 5. / 6.;
     w4= 1 / 6.;

     x1= -1 ;
     x2= -1./ sqrt(5);
     x3= 1./ sqrt(5);
     x4= 1;


    quadrature_Lobatto= ( (b2-a2) / 2. ) * (w1 * g1(x1) + w2 * g1(x2) + w3 * g1(x3) + w4 * g1(x4) );
  }

  else if(n2==5){
     w1= 1./10.;
     w2= 49./90.;
     w3= 32./45.;
     w4= 49./90.;
     w5= 1./10.;

     x1= -1;
     x2= -sqrt( 3./7.);
     x3= 0;
     x4= +sqrt( 3./7.);;
     x5= 1;


    quadrature_Lobatto= ( (b2-a2) / 2. ) * (w1 * g1(x1) + w2 * g1(x2) + w3 * g1(x3) + w4 * g1(x4) + w5 * g1(x5));
  }

  return quadrature_Lobatto;
}
