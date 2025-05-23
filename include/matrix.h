#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <stdexcept>

/*!
    * @file matrix.h
    * @brief interface de la classe matrix
    * @author Claire Tualle, Jaad Belhouari
    */

    class Matrix {

        /*!
         * @brief déclaration des attributs de la classe
         */

    private:
        int _m; //! nombre de colonnes de la matrice
        int _n; //! nombre de lignes de la matrice
        double **_matrix; //! définition d'un tableau 2D

    public:

        /*!
         * @brief : definition des constructeurs de la classe Matrix
         */

        /*!
         * @brief constructeur par défaut
         */
        Matrix();

        /*!
         * @brief : constructeur
         * @param nombre de lignes n, nombre de colonnes m
         * @return construction d'une matrice de taille n,m  initialisée à 0
         */
        Matrix(int n, int m);

        /*!
        * @brief : constructeur
        * @param nombre de lignes n, nombre de colonnes m, tableau à 2 dimensions value
        * @return construction d'une matrice de taille n,m avec les valeurs du tableau value
        */
        Matrix(int n, int m, double **value);

        /*!
         * @brief: constructeur de copie
         * @param matrice
         * @return construction d'une copie de la matrice matrice
         */
        Matrix(const Matrix &matrice);

        /*!
         * @brief destructeur
         * @return "détruit" la matrice
         */

        ~Matrix();

        /*!
         * @brief méthodes implémentées dans la classe Matrix
         */

        /*!
         * @brief getter du nombre de lignes de la matrice
         * @return le nombre de lignes de la matrice
         */
        int getN() const;

        /*!
         * @brief getter du nombre de colonnes de la matrice
         * @return le nombre de colonnes de la matrice
         */
        int getM() const;

        /*!
         * @brief getter du coefficient (i,j) de la matrice
         * @param ligne i, colonne j
         * @return le coefficient (i,j) de la matrice
         */
        double getPoint(int i, int j);

        /*!
         * @brief getter du coefficient (i,j) de la matrice en coordonnées graphique
         * @param ligne i, colonne j, Mx MAx et Mn Min, n taille fenêtre
         * @return le coefficient (i,j) de la matrice
         */
        int getGraphical(int i, int j, double Mx, double Mn, int n);

        /*!
         * @brief getter du coefficient (i,j) de la matrice : rend des fausses couleurs R,G,B
         * @param ligne i, coord j,k , Mx MAx et Mn Min, Esp taille matrice spatiale, n taille fenêtre, RGB
         * @return le coefficient (i,j) de la matrice
         */
        void getColor(int i, int jx, int kx, double Mx, double Mn, int Esp, int n, int& R, int& G, int& B);

        /*!
         * @brief ajout d'une valeur dans la matrice dans à un emplacement précis
         * @param ligne i, colonne j, valeur value
         * @return applique à la matrice la valeur value à la ligne i et la colonne j
         */
        void setPoint(double value, int i, int j);

        /*!
         * @brief méthode retournant la colonne j
         * @param matrice, colonne j
         * @return vecteur correspondant à la colonne j
         */
        double *getColumn(const Matrix matrice, int j);


        /*!
         * @brief méthode attribuant à la colonne j d'une matrice un vecteur donné
         * @param indice j, vecteur colonne
         * @return matrice modifiée
         */
        void setColumn(int j, double *vect);

        /*!
         * @brief méthode retournant la ligne i
         * @param matrice, ligne i
         * @return vecteur correspondant à la ligne i
         */

        double *getLine(const Matrix matrice, int i);

        /*!
         * @brief méthode attribuant à la ligne j d'une matrice un vecteur donné
         * @param indice j, vecteur ligne
         * @return matrice modifiée
         */
        void setLine(int j, double *vect);

         /*!
          * @brief Obtient le minimum des coefficients de la matrice.
          * @return La valeur minimale dans la matrice.
          */
        double min() const;

        /**
         * @brief Obtient le maximum des coefficients de la matrice.
         * @return La valeur maximale dans la matrice.
         */
        double max() const;


        /*!
         * @brief opérateurs internes
         */

        /*!
         * @brief opérateurs internes
         */

        /*!
         * @brief opérateur test égalité
         * @param matrice
         * @return teste si la matrice entrée en argument et la matrice this sont égales en vérifiant qu'elles sont bien de même taille */
        bool operator==(const Matrix &matrice) const;

        /*!
         * @brief opérateur test différence
         * @param matrice
         * @return teste si la matrice entrée en argument et la matrice this sont différentes en vérifiant qu'elles sont bien de même taille */

        bool operator!=(const Matrix &matrice) const;

        /*!
         * @brief opérateur égalité
         * @param matrice
         * @return affecte à la matrice this les valeurs de la matrice matrice
         */
        Matrix operator=(const Matrix &matrice);

        /*!
         * @brief opérateur addition
         * @param matrice
         * @return ajoute à la matrice this les valeurs de la matrice matrice
         */
        Matrix &operator+(const Matrix &matrice);

        /*!
         * @brief opérateur soustraction
         * @param matrice
         * @return soustrait à la matrice this les valeurs de la matrice matrice
         */
        Matrix &operator-(const Matrix &matrice);

        /*!
         * @brief opérateur multiplication scalaire
         * @param value
         * @return multiplication de la valeur value à la matrice this
         */
        Matrix &operator*(double value);

        /*!
         * @brief opérateur division scalaire
         * @param value
         * @return division de la valeur value à la matrice this
         */
        Matrix &operator/(double value);
    };


#endif