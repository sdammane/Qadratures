#include <iostream>
#include "Legendre.hpp"
#include "Lobatto.hpp"
#include <cmath>
#include <vector>
#include <fstream>

double pi=3.141592654;


//Les fonctions à tester

double f1( double x ){
  return pow(x,4);
}

double f2( double x ){
  return 1./sqrt(x);
}

double f3( double x ){
  return (x*x) * log(x);
}

double f4( double x ){
  return cos(x) /sqrt( sin(x) );
}


//La fonction erreur de la méthode de quadrature de Gauss Legendre
double erreur(double soluAp,double soluExacte){
  return fabs( ((soluExacte)-soluAp) );
}


  int main(){

    //Initilisation des bornes de l'intégrale
    double  a1= 1;
    double  b1= 2;
    double  a2= 0;
    double  b2= 1;
    double  a3= 1;
    double  b3= exp(1);
    double  a4= 0;
    double  b4= pi/4.; //En utilisant le "M_PI" cela ne fonctionne pas la commande me disant M_PI was not declared it in this scope


    //Les résultats théoriques des intégrales
    double I1= 31./5.;
    double I2=  2.;
    double I3= ( (2 * exp(3)) + 1) / 9.;
    double I4= pow(2, 3./4.);

    //Les n points de la méthode de quadrature de Gauss et de Lobatto
    int n=2;

    //Les résultats de la quadrature de Legendre des sous intervalles
    double resultat1, resultat2, resultat3, resultat4;

    //Les résultats de la quadrature de Lobatto des sous intervalles
    double res1, res2, res3, res4;



    //Initialisation des sous intervalles
    double N=2; //Le nombre de sous intervalle de [a,b]

    double h1,h2,h3,h4; //Les pas de subdivision pour f1,f2,f3 et f4

    h1= double(b1-a1)/ N;
    h2= double(b2-a2)/ N;
    h3= double(b3-a3)/ N;
    h4= double(b4-a4)/ N;

    //Les bornes des sous intervalles
    //[xi, xi+1] pour f1
    //[yi, yi+1] pour f2
    //[zi, zi+1] pour f3
    //[ti, ti+1] pour f4
    double xi, xiP1,yi, yiP1,zi, ziP1,ti, tiP1;




    //La quadrature sur l'ensemble des sous intervalles pour f1,f2,f3 et f4
    double GaussLegendre1, GaussLegendre2, GaussLegendre3, GaussLegendre4;
    GaussLegendre1=0;
    GaussLegendre2=0;
    GaussLegendre3=0;
    GaussLegendre4=0;


    //La quadrature de Gauss Lobatto sur l'ensemble des sous intervalles
    double GaussLobatto1, GaussLobatto2, GaussLobatto3, GaussLobatto4;
    GaussLobatto1=0;
    GaussLobatto2=0;
    GaussLobatto3=0;
    GaussLobatto4=0;


    QuadratureLegendre Q(n);//Construction d'une quadrature de Legendre avec n point

    QuadratureLobatto L1(n);//Construction d'une quadrature de Lobatto avec n points


    std::cout<<"Les test unitaires sur les quadratures avec "<<n<<" points "<<std::endl;

    for( int i = 0; i<N ; i++ ){

      //Initialisation des sous intervalles des intégrales à évaluer
      xi =   a1+i*h1;
      xiP1 = a1+(i+1)*h1;

      yi =   a2+i*h2;
      yiP1 = a2+(i+1)*h2;

      zi =   a3+i*h3;
      ziP1 = a3+(i+1)*h3;

      ti =   a4+i*h4;
      tiP1 = a4+(i+1)*h4;

      //Initialisation des auadrature de Legendre des sous intervalles [xi,xi+1]
      resultat1 = Q(f1,xi,xiP1);
      resultat2 = Q(f2,yi,yiP1);
      resultat3 = Q(f3,zi,ziP1);
      resultat4 = Q(f4,ti,tiP1);


      //Initialisation des quadrature de Lobatto des sous intervalles [xi,xi+1]
      res1 = L1(f1,xi,xiP1);
      res2 = L1(f2,yi,yiP1);
      res3 = L1(f3,zi,ziP1);
      res4 = L1(f4,ti,tiP1);


      //Les résultats des quadratures de Legendre sur l'ensemble des sous intervalles
      GaussLegendre1 += resultat1;
      GaussLegendre2 += resultat2;
      GaussLegendre3 += resultat3;
      GaussLegendre4 += resultat4;

      //Les résultats des quadratures de Lobatto sur l'ensemble des sous intervalles
      GaussLobatto1 += res1;
      GaussLobatto2 += res2;
      GaussLobatto3 += res3;
      GaussLobatto4 += res4;

    }

    std::cout<<"n:"<<n<<" I1_Gauss:" <<GaussLegendre1<<"\t I2_Gauss:" <<GaussLegendre2<<"\t I3_Gauss:"<<GaussLegendre3<<"\t I4_Gauss:"<<GaussLegendre4<<std::endl;
    std::cout<<std::endl;

    std::cout<<"n:"<<n <<" I1_Lobatto:" <<GaussLobatto1<< "\t I2_Lobatto:" <<"\t"<< GaussLobatto2 <<"\t I3_Lobatto:"<<GaussLobatto3<<"\t I4_Lobatto:"<<GaussLobatto4<<std::endl;
    std::cout<<std::endl;

    std::cout<<"Les erreurs des quadratures avec "<<n<<" points"<<std::endl;
    std::cout<< "Er(I1_Gauss):" <<erreur(GaussLegendre1,I1)<< " Er(I2_Gauss):" <<erreur(GaussLegendre2,I2)<< "\t Er(I3_Gauss):" <<erreur(GaussLegendre3,I3)<< " Er(I4_Gauss):" <<erreur(GaussLegendre4,I4)<<std::endl;
    std::cout<<std::endl;

    std::cout<< "Er(I1_Lobatto):" <<erreur(GaussLobatto1,I1)<< " Er(I2_Lobatto):" <<erreur(GaussLobatto2,I2)<< "\t Er(I3_Lobatto):" <<erreur(GaussLobatto3,I3)<< " Er(I4_Lobatto):" <<erreur(GaussLobatto4,I4)<<std::endl;

    //Remmetre les variables à 0 pour une utilisation
    //des variables ultérieure
    GaussLegendre1=0;
    GaussLegendre2=0;
    GaussLegendre3=0;
    GaussLegendre4=0;

    GaussLobatto1 = 0;
    GaussLobatto2 = 0;
    GaussLobatto3 = 0;
    GaussLobatto4 = 0;




    //Creation des fichiers pour afficher les quadrature des sous intervalles,
    //des quadrature sur [a,b], des erreurs des quadratures


     std::ofstream oFile("quadrature_SousIntervalle_Legendre.dat");
     std::ofstream oFile1("quadrature_SousIntervalle_Lobatto.dat");

      std::ofstream oFile2("quadrature_Legendre.dat");
      std::ofstream oFile3("quadrature_Lobatto.dat");


      std::ofstream oFile4("err_Legendre.dat");
      std::ofstream oFile5("err_Lobatto.dat");




    //faire varier le nombre de points de 2 à 5
    //et le nombre des sous intervalles

    for( int n1 = 2 ; n1 <= 5; n1++ ){
      for( int i = 0; i<N ; i++ ){

        QuadratureLegendre Q2(n1);//Construction d'une quadrature de Legendre avec n points

        QuadratureLobatto L(n1);//Construction d'une quadrature de Lobatto avec n points

        //Initialisation des sous intervalles des intégrales à évaluer
        xi =   a1+i*h1;
        xiP1 = a1+(i+1)*h1;

        yi =   a2+i*h2;
        yiP1 = a2+(i+1)*h2;

        zi =   a3+i*h3;
        ziP1 = a3+(i+1)*h3;

        ti =   a4+i*h4;
        tiP1 = a4+(i+1)*h4;

        //Initialisation des auadrature de Legendre des sous intervalles [xi,xi+1]
        resultat1 = Q2(f1,xi,xiP1);
        resultat2 = Q2(f2,yi,yiP1);
        resultat3 = Q2(f3,zi,ziP1);
        resultat4 = Q2(f4,ti,tiP1);

        //Initialisation des quadrature de Lobatto des sous intervalles [xi,xi+1]
        res1 = L(f1,xi,xiP1);
        res2 = L(f2,yi,yiP1);
        res3 = L(f3,zi,ziP1);
        res4 = L(f4,ti,tiP1);


        //Les résultats des quadratures de Legendre sur l'ensemble des sous intervalles
        GaussLegendre1 += resultat1;
        GaussLegendre2 += resultat2;
        GaussLegendre3 += resultat3;
        GaussLegendre4 += resultat4;

        //Les résultats des quadratures de Lobatto sur l'ensemble des sous intervalles
        GaussLobatto1 += res1;
        GaussLobatto2 += res2;
        GaussLobatto3 += res3;
        GaussLobatto4 += res4;


        oFile << "n:" <<n1<<"\t" <<"I1_Gauss sur "<<"[ "<<xi<<","<<xiP1<<"]\t:"<< resultat1 <<"\t I2_Gauss"<<"["<<yi<<","<<yiP1<<"]\t:" << resultat2 <<"\t I3_Gauss"<<"["<<zi<<","<<ziP1<<"]\t:"<< resultat3 <<"\t I4_Gauss"<< "["<<ti<<","<<tiP1<<"]\t:"<< resultat4 << std::endl;
        oFile1 << "n:" <<n1<<"\t" <<"I1_Lobatto sur "<<"[ "<<xi<<","<<xiP1<<"]\t:"<< res1 << "\t I2_Lobatto"<<"["<<yi<<","<<yiP1<<"]\t:" << res2 << "\t I3_Lobatto"<< "["<<zi<<","<<ziP1<<"]\t:" << res3 << "\t I4_Lobatto"<< "["<<ti<<","<<tiP1<<"]\t:" << res4 << std::endl;
      }


      //Quadrature de Legendre sur l'ensemble des sous intervalles de [a,b]
      oFile2<<"n:"<<n1<<"\t I1_Gauss:" <<GaussLegendre1<<"\t I2_Gauss:" <<GaussLegendre2<<"\t I3_Gauss:"<<GaussLegendre3<<"\t I4_Gauss:"<<GaussLegendre4<<std::endl;

      //Les erreurs de quadrature de Legendre sur l'ensemble des sous intervalles de [a,b]
      oFile4<<"n:"<<n1<<"\t"<< "Er(I1_Gauss):\t" <<erreur(GaussLegendre1,I1)<< "\t Er(I2_Gauss):" <<erreur(GaussLegendre2,I2)<< "\t Er(I3_Gauss):" <<erreur(GaussLegendre3,I3)<< "\t Er(I4_Gauss):" <<erreur(GaussLegendre4,I4)<<std::endl;


      //Quadrature de Lobatto sur l'ensemble des sous intervalles de [a,b]
      oFile3<<"n:"<<n1<<"\t I1_Lobatto:" <<GaussLobatto1<< "\t I2_Lobatto:" <<GaussLobatto2<<"\t I3_Lobatto:"<<GaussLobatto3<<"\t I4_Lobatto:"<<GaussLobatto4<<std::endl;

      //Les erreurs de quadrature de Legendre sur l'ensemble des sous intervalles de [a,b]
      oFile5<<"n:"<<n1<<"\t"<< "Er(I1_Lobatto):\t" <<erreur(GaussLobatto1,I1)<< "\t Er(I3_Lobatto):" <<erreur(GaussLobatto3,I3)<<std::endl;



      //Remmetre les variables à 0 sinon le résultat sera fausssé //quel résultat
      GaussLegendre1=0;
      GaussLegendre2=0;
      GaussLegendre3=0;
      GaussLegendre4=0;

      GaussLobatto1=0;
      GaussLobatto2=0;
      GaussLobatto3=0;
      GaussLobatto4=0;
    }

    return 0;
  }
