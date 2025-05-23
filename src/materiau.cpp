#include "materiau.h"

/*!
 *
 * @param param
 * @param L
 * @param t_max
 * @param u_0
 * @param f
 * @param dim
 */

//constructeur

Materiau::Materiau() {
    this -> _lambda = 0;
    this -> _rho = 0;
    this -> _c = 0;
}

Materiau::Materiau(double lambda, double rho, double c) {
    _lambda = lambda;
    _rho = rho;
    _c = c;
}

/*!
 * @brief déclaration des méthodes
 */


void Materiau::setCuivre(){
    this->_lambda = 389;
    this->_rho = 8940;
    this->_c = 380;
}

void Materiau::setFer(){
    this->_lambda = 80.2;
    this->_rho = 7874;
    this->_c = 440;
}

void Materiau::setVerre(){
    this->_lambda = 1.2;
    this->_rho = 2530;
    this->_c = 840;
}

void Materiau::setPolystyrene(){
    this->_lambda = 0.1;
    this->_rho = 1040;
    this->_c = 1200;
}

double Materiau::getLambda() const {
    return this->_lambda;
}

double Materiau::getRho() const {
    return this->_rho;
}

double Materiau::getC() const {
    return this->_c;
}
