#include <iostream>
#include "Legendre.hpp"
#include <cmath>



//Le constructeur qui initialise une quadrature à n points
QuadratureLegendre::QuadratureLegendre(int n1){
  n=n1;
}

//La fonction de changement de variable pour se ramener à
//l'intervalle [-1,1]
double QuadratureLegendre::g(double t){

    return f( ( (b-a)*t+(b+a) ) /2.);
}


//Surchage de l'operateur ()
double QuadratureLegendre:: operator() (double (*g1) (double),double a1,double b1){
  //Initialisation des bornes de l'intervalle
  a=a1;
  b=b1;
  //Initialiation de la fonction
  f=g1;

  return Gauss_Legendre() ;
}



//La fonction qui calcule la quadrature de Gauss Legendre

double QuadratureLegendre::Gauss_Legendre(){

  //Les poids associés aux points xi
  double w1,w2,w3,w4,w5;
  double x1,x2,x3,x4,x5;

  //Variable pour la quadrature de Gauss Legendre
  double quadrature;



  if(n==2){
     w1=1;
     w2=1;

     x1= 1. / sqrt( 3 );
     x2= -1. / sqrt( 3 );


    quadrature= ( (b-a) / 2. ) * (w1 * g(x1) + w2 * g(x2) );
  }

  else if(n==3){
     w1= 5./9.;
     w2= 8./9.;
     w3= 5./9.;

     x1= -sqrt( 0.6 );
     x2=0;
     x3= sqrt( 0.6 );


    quadrature= ( (b-a) / 2. ) * (w1 * g(x1) + w2 * g(x2) + w3 * g(x3) );
  }

  else if(n==4){
     w1= ( 18 + sqrt( 30 )) / 36.;
     w2= ( 18 - sqrt( 30 )) / 36.;
     w3= ( 18 -sqrt( 30 )) / 36.;
     w4= ( 18 + sqrt( 30 )) / 36.;

     x1= -sqrt( 3. / 7. + (2 * sqrt(6) )/ (7*sqrt(5)) ) ;
     x2= -sqrt( 3. / 7. - (2 * sqrt(6) )/ (7*sqrt(5)) ) ;
     x3= +sqrt( 3. / 7. - (2 * sqrt(6) )/ (7*sqrt(5)) ) ;
     x4= +sqrt( 3. / 7. + (2 * sqrt(6) )/ (7*sqrt(5)) ) ;


    quadrature= ( (b-a) / 2. ) * (w1 * g(x1) + w2 * g(x2) + w3 * g(x3) + w4 * g(x4) );
  }

  else if(n==5){
     w1= ( 322 + 13 * sqrt(70) ) / 900.;
     w2= ( 322 - 13 * sqrt(70) ) / 900.;
     w3= 128. / 225.;
     w4= ( 322 - 13 * sqrt(70) ) / 900.;
     w5= ( 322 + 13 * sqrt(70) ) / 900.;

     x1= ( -1. / 3. ) * sqrt ( 5 + 2 * sqrt (10. / 7.));
     x2= ( -1. / 3. ) * sqrt ( 5 - 2 * sqrt (10. / 7.));
     x3= 0;
     x4= ( 1. / 3. ) * sqrt ( 5 - 2 * sqrt (10. / 7.));
     x5= ( 1. / 3. ) * sqrt ( 5 + 2 * sqrt (10. / 7.));


    quadrature= ( (b-a) / 2. ) * (w1 * g(x1) + w2 * g(x2) + w3 * g(x3) + w4 * g(x4) + w5 * g(x5));
  }

  return quadrature;
}
