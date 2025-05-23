#ifndef MATERIAU_H
#define MATERIAU_H
#include <iostream>


/*!
    * @file materiau.h
    * @brief interface de la classe materiau
    * @author Claire Tualle, Jaad Belhouari 
    */


class Materiau {

    /*!
     * @brief déclaration des attributs de la classe materiau
     */

private:
    double _lambda; //! conductivité thermique
    double _rho; //! masse volumique 
    double _c; //! chaleur massique

public:

    /*!
     * @brief : definition des constructeurs de la classe materiau
     */

    /*!
     * @brief constructeur par défaut
     */
    Materiau();

    /*!
     * @brief : constructeur
     * @param conductivité thermique, masse volumique et chaleur massique
     * @return construction des coefficients représentatifs du matériau
     */
    Materiau(double lambda, double rho, double c);

    /*!
     * @brief méthodes implémentées dans la classe materiau
     */

    /*!
     * @brief le cuivre
     * @return les grandeurs physiques du cuivre
     */
    void setCuivre();

    /*!
     * @brief le fer
     * @return les grandeurs physiques du fer
     */
    void setFer();

    /*!
     * @brief le verre
     * @return les grandeurs physiques du verre
     */
    void setVerre();

    /*!
     * @brief le polystyrène
     * @return les grandeurs physiques du polystyèrne
     */
    void setPolystyrene();

    /*!
     * @brief getter de la conductivité thermique du matériau
     * @return la conductivité thermique du matériau
     */
    double getLambda() const;

    /*!
     * @brief getter de la masse volumique du matériau
     * @return la masse volumique du matériau
     */
    double getRho() const;

    /*!
     * @brief getter de la chaleur massique du matériau
     * @return la chaleur massique du materiau
     */
    double getC() const;

};


#endif


