#include "matrix.h"


//constructeur

Matrix::Matrix(){
    this -> _n = 0;
    this -> _m = 0;
    this -> _matrix = NULL;
}

Matrix::Matrix(int n, int m) {
    this -> _m = m;
    this -> _n = n;
    this -> _matrix = new double*[n];
    for (int i = 0; i < n; i++) {
        this -> _matrix[i] = new double[m];
        for (int j = 0; j < m; j++) {
            this -> _matrix[i][j] = 0;
        }
    }
}

Matrix::Matrix(int n, int m, double **value) {
    this->_n = n;
    this->_m = m;
    this->_matrix = new double *[m];
    for (int i = 0; i < n; i++) {
        this -> _matrix[i] = new double[n];
        for (int j = 0; j < m; j++) {
            this->_matrix[i][j] = value[i][j];
        }
    }
}


Matrix::Matrix(const Matrix &matrice) {
    this->_n = matrice._n;
    this->_m = matrice._m;
    this->_matrix = new double *[_m];
    for (int i = 0; i < _n; i++) {
        this -> _matrix[i] = new double[_n];
        for (int j = 0; j < _m; j++) {
            this -> _matrix[i][j] = matrice._matrix[i][j];
        }
    }
}

//destructeur

Matrix::~Matrix() {
    int n = this->_n;
    for (int i = 0; i < n; i++) {
        delete[] _matrix[i];
    }
    delete[] _matrix;
}


//méthodes
int Matrix::getN() const {
    return this->_n;
}

int Matrix::getM() const {
    return this->_m;
}

double Matrix::getPoint(int i, int j) {
    return this->_matrix[i][j];
}

int Matrix::getGraphical(int i, int jx, double Mx, double Mn, int n) {
    int m = this->_m;
    int j = (m*jx)/n;
    double x = this->_matrix[i][j];
    return (int)(n*(x-Mn)/(Mx-Mn));
}

void Matrix::getColor(int i, int jx, int kx, double Mx, double Mn, int Esp, int n, int& R, int& G, int& B) {
   
    int j = (Esp*jx)/n;
    int k = (Esp*kx)/n;
    double x = this->_matrix[i][j*Esp+k];

    x = 255*(x-Mn)/(Mx-Mn);
    R = (int)(x);
    G = 0;
    B = (int)(255-x);

    return;
}


void Matrix::setPoint(double value, int i, int j) {
    this->_matrix[i][j] = value;
}

double* Matrix::getColumn(Matrix matrice, int j) {
    int n = matrice.getN();
    double* vect = NULL;
    for (int i = 0; i < n; i++) {
        vect[i] = matrice.getPoint(i,j);
    }
    return vect;
}

void Matrix::setColumn(int j, double* vect) {
    int n = this -> _n;
    for (int i = 0; i < n; i++) {
        this->_matrix[i][j] = vect[i];
    }
}

double* Matrix::getLine(const Matrix matrice, int i) {
    return matrice._matrix[i];
}

void Matrix::setLine(int j, double* vect) {
    int m = this -> _m;
    for (int i = 0; i < m; i++) {
        this->_matrix[j][i] = vect[i];
    }
}

double Matrix::min() const {
    int n = this->_n;
    int m = this->_m;
    double minVal = this->_matrix[0][0];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (this->_matrix[i][j] < minVal) {
                minVal = this->_matrix[i][j];
            }
        }
    }

    return minVal;
}

double Matrix::max() const {
    int n = this->_n;
    int m = this->_m;
    double maxVal = this->_matrix[0][0];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (this->_matrix[i][j] > maxVal) {
                maxVal = this->_matrix[i][j];
            }
        }
    }
    return maxVal;
}


bool Matrix::operator==(const Matrix &matrice) const {
    int m = this -> _m;
    int n = this -> _n;
    bool boolean = true;
    if (m != matrice._m || n != matrice._n) {
        throw std::invalid_argument("Erreur : les dimensions des matrices posent problème");
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (this->_matrix[i][j] != matrice._matrix[i][j]) {
                boolean = false;
            }
        }
    }
    return boolean;
}


bool Matrix::operator!=(const Matrix &matrice) const {
    int m = this->_m;
    int n = this->_n;
    bool boolean = true;
    if (m != matrice._m || n != matrice._n) {
        throw std::invalid_argument("Erreur : les dimensions des matrices posent problème");
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (this->_matrix[i][j] == matrice._matrix[i][j]) {
                boolean = false;
            }
        }
    }
    return boolean;
}

Matrix Matrix::operator=(const Matrix &matrice) {
    int m = this->_m;
    int n = this->_n;
    if (m != matrice._m || n != matrice._n) {
        throw std::invalid_argument("Erreur : les dimensions des matrices posent problème");
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            this->_matrix[i][j] = matrice._matrix[i][j];
        }
    }
    return *this;
}



Matrix &Matrix::operator+(const Matrix &matrice) {
    int m = this->_m;
    int n = this->_n;
    if (m != matrice._m || n != matrice._n) {
        throw std::invalid_argument("Erreur : les dimensions des matrices posent problème");
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
           this -> _matrix[i][j] = this -> _matrix[i][j] + matrice._matrix[i][j];
        }
    }
    return *this;
}



Matrix &Matrix::operator-(const Matrix &matrice) {
    int m = this->_m;
    int n = this->_n;
    if (m != matrice._m || n != matrice._n) {
        throw std::invalid_argument("Erreur : les dimensions des matrices posent problème");
    }
    this -> _matrix = new double*[n];
    for (int i = 0; i < n; i++) {
        this -> _matrix[i] = new double[m];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            this -> _matrix[i][j] = this->_matrix[i][j] - matrice._matrix[i][j];
        }
    }
    return *this;

}
Matrix &Matrix::operator*(double value) {
    int m = this->_m;
    int n = this->_n;
    this -> _matrix = new double*[n];
    for (int i = 0; i < n; i++) {
        this -> _matrix[i] = new double[m];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            this -> _matrix[i][j] = value * this->_matrix[i][j];
        }
    }
    return *this;
}


Matrix &Matrix::operator/(double value) {
    int m = this->_m;
    int n = this->_n;
    this -> _matrix = new double*[n];
    for (int i = 0; i < n; i++) {
        this -> _matrix[i] = new double[m];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            this -> _matrix[i][j] = this->_matrix[i][j] / value;
        }
    }
    return *this;
}

