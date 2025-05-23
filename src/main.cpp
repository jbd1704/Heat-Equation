#include <iostream>
#include "edp.h"
#include "graphique.h"

int main() {

    //test();
    int N = 1001;
    int TIME = 100;
    //! Premier cas : la barre
    Materiau cu(0,0,0);
    Materiau fe(0,0,0);
    Materiau ver(0,0,0);
    Materiau poly(0,0,0);

    cu.Materiau::setCuivre();
    fe.Materiau::setFer();
    ver.Materiau::setVerre();
    poly.Materiau::setPolystyrene();

    Edp barre(1, 16, 286.15, 353.15);            //! Les températures sont en Kelvin, en accord avec les coefficients des matériaux

    
    Matrix U_1_cuivre = barre.Edp::getSolution_1(N, TIME, cu);
    Matrix U_1_fer = barre.Edp::getSolution_1(N, TIME, fe);
    Matrix U_1_verre = barre.Edp::getSolution_1(N, TIME, ver);
    Matrix U_1_polystyrene = barre.Edp::getSolution_1(N, TIME, poly); 

    test(U_1_cuivre, "Cuivre");
    test(U_1_fer, "Fer");
    test(U_1_verre, "Verre");
    test(U_1_polystyrene, "Polystyrene");


   /* for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
        std::cout << "temps: "<< i << ", espace: "<< j << "=" << U_1_fer.Matrix::getPoint(j, i) << std::endl;
    }
    } 
*/

    //! Deuxième cas : la plaque

    Edp plaque(1, 16, 286.15, 353.15);


    Matrix U_2_cuivre = plaque.Edp::getSolution_2(TIME, N, cu);
    Matrix U_2_fer = plaque.Edp::getSolution_2(TIME, N, fe);
    Matrix U_2_verre = plaque.Edp::getSolution_2(TIME, N,ver);
    Matrix U_2_polystyrene = plaque.Edp::getSolution_2(TIME, N, poly);
    
    test2D( U_2_polystyrene, N, "Polystyrene");
    test2D( U_2_verre, N, "Verre");
    test2D( U_2_fer, N, "Fer");
    test2D( U_2_cuivre, N, "Cuivre");
    
    std::cout << "Nous n'avons pas reussi a coder la methode implicite pour la plaque" << std::endl;
    std::cout << "Et la methode explicite diverge... (on ne voit que les sommets des zones chauffantes)" << std::endl;
    
    

/*    for (int i = 0; i < TIME; i++){
        std::cout << U_2_cuivre.Matrix::getPoint(i, 600600) << std::endl;
    }*/

//    int a ;
//    std::cin >> a ;


    return 0;
}