#ifndef EDP_H
#define EDP_H
#include <iostream>
#include <math.h>
#include "matrix.h"
#include "materiau.h"

/*!
    * @file edp.h
    * @brief interface de la classe edp
    * @author Claire Tualle, Jaad Belhouari 
    */



class Edp {

    /*!
     * @brief déclaration des attributs de la classe
     */

private:
    double _L; //! longueur de l'objet
    double _tmax; //! durée de la propagation de la chaleur
    double _u_0; //! température initiale
    double _f; //! ajout de chaleur

public:

    /*!
     * @brief : definition des constructeurs de la classe edp
     */

    /*!
     * @brief constructeur par défaut
     * */

    Edp();
    /*!
     * @brief : constructeur
     * @param coefficients du matériau, dimension de l'objet, temps de la propagation, température initiale, source de chaleur
     * @return construction de l'équation différentielle représentative du problème
     */
    Edp(double L, double tmax, double u_0, double f);

    /*!
     * @brief méthodes implémentées dans la classe edp
     */

    /*!
     * @brief source de chaleur pour le problème à une dimension
     * @return la température créée au point x
     */
    double F_1(double x) const;

    /*!
     * @brief source de chaleur pour le problème à deux dimensions
     * @return la température créée au point (x, y)
     */
    double F_2(double x, double y) const;

    /*!
     * @brief getter de la solution du problème à une dimension en fonction du nombre de points choisis temporellement et spatialement
     * @return la solution de l'équation différentielle à une dimension
     */

    // Dans notre cas, N_temps = N_espace = N = 1001

    Matrix getSolution_1(int N_temps, int N_espace, Materiau param) const;

    /*!
    * @brief getter de la solution du problème à deux dimensions en fonction du nombre de points choisis temporellement et spatialement
    * @return la solution de l'équation différentielle à deux dimensions
    */
    Matrix getSolution_2(int N_temps, int N_espace, Materiau param) const;
};


#endif


