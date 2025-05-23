#include "edp.h"
#include "matrix.h"
#include "materiau.h"


/*!
 *
 * @param param
 * @param L
 * @param t_max
 * @param u_0
 * @param f
 */

//constructeur
Edp::Edp(double L, double tmax, double u_0, double f) {
    _L = L;
    _tmax = tmax;
    _u_0 = u_0;
    _f = f;
}

Edp::Edp(){
    this -> _L = 0;
    this -> _tmax = 0;
    this -> _u_0 = 0;
    this -> _f = 0;
}

//méthodes

double Edp::F_1(double x) const {
    double tmax = this->_tmax;
    double L = this->_L;
    double f = this->_f;
    double res = 0;
    if (x >= L / 10 && x <= 2 * L / 10) {
        res = tmax * pow(f,2);
    }
    if (x >= 5 * L / 10 && x <= 6 * L / 10) {
        res = 3 * tmax / 4 * pow(f,2);
    }
    return res;
}

double Edp::F_2(double x, double y) const {
    double tmax = this->_tmax;
    double L = this->_L;
    double f = this->_f;
    double res = 0;
    if (x >= L / 6 && x <= 2 * L / 6 && y >= L / 6 && y <= 2 * L / 6) {
        res = tmax * pow(f,2);
    }
    if (x >= 4 * L / 6 && x <= 5 * L / 6 && y >= L / 6 && y <= 2 * L / 6) {
        res = tmax * pow(f,2);
    }
    if (x >= L / 6 && x <= 2 * L / 6 && y >= 4 * L / 6 && y <= 5 * L / 6) {
        res = tmax * pow(f,2);
    }
    if (x >= 4 * L / 6 && x <= 5 * L / 6 && y >= 4 * L / 6 && y <= 5 * L / 6) {
        res = tmax * pow(f,2);
    }
    return res;
}

Matrix Edp::getSolution_1(int N_espace, int N_temps, Materiau param) const {
    double tmax = this->_tmax;
    double u_0 = this->_u_0;
    double L = this->_L;
    double dt = tmax / N_temps;
    double dx = L / N_espace;
    double beta = dt * param.getLambda() / (param.getRho() * param.getC() * pow(dx,2));
    Matrix U(N_temps , N_espace);                    //! La matrice résultante qui contiendra les températures
    Matrix B(N_espace, N_espace);                  //! Matrice tridiagonale
    double C[N_espace];                                   //! Vecteur colonne contenant les termes constants correspondants à la production de chaleur
    double U_j[N_espace];                                 //! Colonne j de U : les températures de tout point à l'instant j
    double U_j_next[N_espace];
    //! Remplissage des matrices permettant le calcul et contenant les données initiales
    for (int i = 0; i < N_espace; i++){
         U_j[i] = u_0;                //! Condition initiale
         C[i] = beta * pow(dx,2) * F_1(i * dx) / param.getLambda();
    }
    for (int i = 2; i < N_espace - 1; i++) { //! On ne prend ici pas en compte les bords
        B.Matrix::setPoint(1 + 2 * beta, i, i);
        B.Matrix::setPoint(-beta, i, i-1);
        B.Matrix::setPoint(-beta, i-1, i);
    }
    B.Matrix::setPoint(1 + beta, 1 , 1); //! Condition de Neumann
 
    U.setLine(0, U_j);

    double alpha[N_espace];
    double gamma[N_espace];
    double D[N_espace];
    for(int j = 1; j < N_temps; j++){
        for (int i = 0; i < N_espace; i ++){
            D[i] = U_j[i] - C[i];
        }
        U_j_next[N_espace - 1] = u_0;                           //! Condition de Dirichlet
        D[N_espace - 2] += beta*u_0; 
        alpha[1] = B.Matrix::getPoint(1,1);
        gamma[1] = D[1] / alpha[1];
        for (int i = 2; i < N_espace - 1; i++){             //! Algorithme de Thomas
            alpha[i] = B.Matrix::getPoint(i,i) - B.Matrix::getPoint(i,i-1)*B.Matrix::getPoint(i-1,i)/alpha[i-1];
            gamma[i] = (D[i] - B.Matrix::getPoint(i-1,i)*gamma[i-1])/alpha[i];
        }
        U_j_next[N_espace-2] = gamma[N_espace-2];
        for (int i = N_espace-3; i > 0; i--){
            U_j_next[i] = gamma[i] - B.Matrix::getPoint(i,i+1)*U_j_next[i+1]/alpha[i];
        }
        U_j_next[0] = U_j_next[1];                              //! Condition de Neumann
    
        U.setLine(j, U_j_next);


   // Vérification de la méthode de Thomas

    /*    double S=0.0;
        double foo;
        for(int i = 1; i < N_espace - 1; i ++)
        {
            foo = (1+2*beta)*U_j_next[i]-beta*(U_j_next[i+1]+U_j_next[i-1]);
            foo = foo - (U_j[i] - C[i]);
            S = S + foo*foo;
        }  

        std::cout << "Voici le résultat ("<< j << ") de S, S =  " << S << std::endl;
     */


   ///

        for (int i = 0; i < N_espace; i ++){
            U_j[i] = U_j_next[i];
        }
    }

    return U;
}



/*Méthode implicite : résolution du système par calcul matriciel -> abandonée car calcul de l'inverse trop coûteux
 *
 * Matrix Edp::getSolution_1(int N_temps, int N_espace, Materiau param) const {
    double tmax = this->_tmax;
    double u_0 = this->_u_0;
    double L = this->_L;
    double dt = tmax / N_temps;
    double dx = L / N_espace;
    double beta = dt * param.getLambda() / (param.getRho() * param.getC() * dx * *2);
    Matrix C = Matrix(N_temps, 1);
    Matrix U = Matrix(N_temps, N_espace);
    Matrix B = Matrix(N_temps, N_espace);
    double* U_j[N_espace];
    double* U_j_next[N_espace];
    //! Remplissage des matrices permettant le calcul et contenant les données initiales
    for (int i = 0; i < N_espace; i++) {
        C._matrix[i][0] = beta * dx * *2 * F_1(i * dx) / param.getLambda();        //! x = i*L/N = i*dx
        B._matrix[i][i] = 1 + 2 * beta;
        U_j[i] = u0;                //! Condition initiale
        for (int j = 0; j < N_temps, j++) {
            if (i == j - 1 || i == j + 1) {
                B._matrix[i][j] = -beta;
            }
        }
    }
    U.setColumn(0, U_j);
    //! Calcul de la température pour chaque instant et remplissage de U
    for (int J = 1; J < N_temps; J++) {
        U.setPoint(u0, L, J);       //! Condition de Dirichlet
        for (int i = 0; i < N_espace; i++) {
                double tmp = 0;
                double tmp_bis = 0;
                for (int k = 0; k < N_espace; k++){
                    tmp = tmp + inv_B._matrix[i][k] * U_j[k];         //! On calcule les coefficients du produit B^-1*U_j
                    tmp_bis = tmp_bis + inv_B._matrix[i][k] * C[k];   //! On calcule les coefficeints du produit B^-1*C
                }
                U_j_next[i] = tmp - tmp_bis;              //! La différence des deux correspond au terme i de U_j+1
            }
        U.setColumn(J, U_j_next);
        U_j = U_j_next;
    }
    return U;
}
*/


/* Méthode explicite -> pas ce qui est demandé
 *
 * Matrix Edp::getSolution_1(int N_temps, int N_espace, Materiau param) const{
    double t_max = this -> _t_max;
    double L = this -> _L;
    double u_0 = this -> _u_0;
    Matrix u = Matrix(N_temps, N_espace, 0);
    dt = t_max / N_temps;
    dx = L / N_espace;
    // Conditions initiales
    for (int i = 0; i < N_espace; i++){
        u[i][0] = u_0;
    }
    for (int j = 0; j <*//* N_temps; j++){
        u[N_espace][j] = u_0; 
    }
    // Remplissage de la matrice u
    for (int i = 1; i < N_espace - 1; i++){
        for (int j = 0; j < N_temps - 2; j++){
            u[i][j+1] = u[i][j] + dt/(param.getRho() * param.getC())*(param.getLambda()*((u[i+1][j] - 2*u[i][j] + u[i-1][j])/dx**2) + F_1(i*dx));      // x = i*L/n = i*dx
        }
    }
    return u;
};*/


/*
void Edp::MatSolve( double *diagG, double *idownG, double *iupG, , double *jdownG, double *jupG, double *D, double *U, int n) {

    double x,a;
    double dG(n*n);                         //! Diagonale de la matrice G
    double idG(n*n);                        //! G(i+1,j;i,j)
    double jdG(n*n);                        //! G(i,j+1;i,j)
    double iuG(n*n);                          //! G(i,j;i+1,j)
    double juG(n*n);                          //! G(i,j;i,j+1)

    //! Copie de la matrice creuse

    for(int i=0; i<n*n; i++){
        dG[i] = diagG[i];
        idG[i] = idownG[i];
        jdG[i] = jdownG[i];
        iuG[i] = iupG[i];
        juG[i] = jupG[i];
    }

    for(int di=1; di<n-1; di++){
        for(int dj=1; dj<n-1; dj++) {
            index = di*n+dj;
            x = dG[index];


        }
    }

}
*/

/*
Matrix Edp::getSolution_2(int N_temps, int N_espace, Materiau param) const {
    double tmax = this->_tmax;
    double u_0 = this->_u_0;
    double L = this->_L;
    double dt = tmax / N_temps;
    double dx = L / N_espace;                       //! Delta x = Delta y
    double mu = dt * param.getLambda() / (param.getRho() * param.getC() * pow(dx,2));
    Matrix U(N_temps , N_espace*N_espace);                    //! La matrice résultante qui contiendra les températures
    cMatrix G;                                                // Matrice creuse                        
    double *H;                                   //! Vecteur colonne contenant les termes constants correspondants à la production de chaleur
    double *U_k;                                 //! Colonne k de U : les températures de tout point à l'instant k
    double *U_k_next;
    double *D;
    double foo;
    int index;
    int Nm2 = N_espace-2;

    std::cout << "Mu = " << mu <<  std::endl;

    H = (double *) calloc(N_espace*N_espace, sizeof(double));
    U_k = (double *) calloc(N_espace*N_espace, sizeof(double));
    U_k_next = (double *) calloc(N_espace*N_espace, sizeof(double));
    D = (double *) calloc(Nm2*Nm2, sizeof(double));

    //! Remplissage des matrices permettant le calcul et contenant les données initiales
    for (int i = 0; i < N_espace; i++) {
        for (int j = 0; j < N_espace; j++) {
            index = i * N_espace + j; 
            H[index] = F_2(i * dx, j * dx) / (param.getRho() * param.getC());
            U_k[index] = u_0;                //! Condition initiale
        }
    }

    for (int i = 1; i < Nm2; i++) {
        for (int j = 1; j < Nm2; j++) {
            index = i * Nm2 + j;
            G.setPoint(1+4*mu,index,index);
            G.setPoint(-mu,index-1,index);
            G.setPoint(-mu,index,index-1);
            G.setPoint(-mu,index-Nm2,index);
            G.setPoint(-mu,index,index-Nm2);
        }
    }

    for (int i = 1; i < Nm2; i++) {
            index = i * Nm2 ;
            G.setPoint(1+3*mu,index,index);   //! Condition de Neumann
            G.setPoint(-mu,index-Nm2,index);
            G.setPoint(-mu,index,index-Nm2);

            index = i ;
            G.setPoint(1+3*mu,index,index);
            G.setPoint(-mu,index-1,index);
            G.setPoint(-mu,index,index-1);
    }

    G.setPoint(1+2*mu,0,0);  //! Condition de Neumann

    std::cout << "condition initiale : " << u_0  << std::endl;
    U.setLine(0, U_k);
    
    for(int ntime=1; ntime<N_temps; ntime++){

        //! Calcul du modèle explicite

        for (int i = 1; i < N_espace-1; i++) {
                for (int j = 1; j < N_espace-1; j++) {
                    index = i * N_espace + j; 
                    foo = (1-4*mu)*U_k[index];
                    foo += mu*(U_k[index-1]+U_k[index+1]+U_k[index-N_espace]+U_k[index+N_espace]) + H[index];         
                    U_k_next[index] = foo;
                }
            }

        //! Condition de Dirichlet

        for (int i = 0; i < N_espace; i++) { 
            index = (N_espace-1)*N_espace + i;
            U_k_next[index] = u_0;
            index = i*N_espace + N_espace-1;
            U_k_next[index] = u_0;
        }
        //! Condition de Neumann

        for (int i = 1; i < N_espace; i++) { 
            index = i;
            U_k_next[index] = U_k_next[index + N_espace];
            index = i*N_espace; 
            U_k_next[index] = U_k_next[index + 1];
        }
        U_k_next[0] = U_k_next[N_espace];

        //! Enregistrement du résultat et copie
        U.setLine(ntime, U_k_next);
        
        for (int i = 0; i < N_espace*N_espace; i++) {
            U_k[i] = U_k_next[i];
        }
    }

    free(U_k_next);
    free(U_k);
    free(H);

    return U;
 }*/


 Matrix Edp::getSolution_2(int N_temps, int N_espace, Materiau param) const {
    double tmax = this->_tmax;
    double u_0 = this->_u_0;
    double L = this->_L;
    double dt = tmax / N_temps;
    double dx = L / N_espace;                       //! Delta x = Delta y
    double mu = dt * param.getLambda() / (param.getRho() * param.getC() * pow(dx,2));
    Matrix U(N_temps , N_espace*N_espace);                    //! La matrice résultante qui contiendra les températures
    double *H;                                   //! Vecteur colonne contenant les termes constants correspondants à la production de chaleur
    double *U_k;                                 //! Colonne k de U : les températures de tout point à l'instant k
    double *U_k_next;
    double foo;
    int index;

    std::cout << "Mu = " << mu <<  std::endl;

    H = (double *) calloc(N_espace*N_espace, sizeof(double));
    U_k = (double *) calloc(N_espace*N_espace, sizeof(double));
    U_k_next = (double *) calloc(N_espace*N_espace, sizeof(double));

    //! Remplissage des matrices permettant le calcul et contenant les données initiales
    for (int i = 0; i < N_espace; i++) {
        for (int j = 0; j < N_espace; j++) {
            index = i * N_espace + j; 
            H[index] = F_2(i * dx, j * dx) / (param.getRho() * param.getC());
            U_k[index] = u_0;                //! Condition initiale
        }
    }
 
    std::cout << "condition initiale : " << u_0  << std::endl;
    U.setLine(0, U_k);
    
    for(int ntime=1; ntime<N_temps; ntime++){

        //! Calcul du modèle explicite

        for (int i = 1; i < N_espace-1; i++) {
                for (int j = 1; j < N_espace-1; j++) {
                    index = i * N_espace + j; 
                    foo = (1-4*mu)*U_k[index];
                    foo += mu*(U_k[index-1]+U_k[index+1]+U_k[index-N_espace]+U_k[index+N_espace]) + H[index];         
                    U_k_next[index] = foo;
                }
            }

        //! Condition de Dirichlet

        for (int i = 0; i < N_espace; i++) { 
            index = (N_espace-1)*N_espace + i;
            U_k_next[index] = u_0;
            index = i*N_espace + N_espace-1;
            U_k_next[index] = u_0;
     }
        //! Condition de Neumann

        for (int i = 1; i < N_espace; i++) { 
            index = i;
            U_k_next[index] = U_k_next[index + N_espace];
            index = i*N_espace; 
            U_k_next[index] = U_k_next[index + 1];
        }
        U_k_next[0] = U_k_next[N_espace];

        //! Enregistrement du résultat et copie
        U.setLine(ntime, U_k_next);
        
        for (int i = 0; i < N_espace*N_espace; i++) {
            U_k[i] = U_k_next[i];
        }
    }

    free(U_k_next);
    free(U_k);
    free(H);

    return U;
 }

