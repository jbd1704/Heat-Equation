#ifndef GRAPHIQUE_H_INCLUDED
#define GRAPHIQUE_H_INCLUDED

/**
 * @file graphique.h
 * @brief Définition des fonctions pour le dessin graphique.
 */

#include <stdint.h>
#include "matrix.h"
#include <SDL2/SDL.h>

/**
 * @brief Dessine une ligne entre deux points avec SDL_Renderer.
 * @param renderer Le SDL_Renderer sur lequel dessiner la ligne.
 * @param ix La coordonnée x du point de départ.
 * @param iy La coordonnée y du point de départ.
 * @param fx La coordonnée x du point d'arrivée.
 * @param fy La coordonnée y du point d'arrivée.
 */
    void sdlLine(SDL_Renderer *r, int xi, int yi, int xf, int yf);


/**
 * @brief Fonction de test pour le dessin.
 * 
 * Cette fonction initialise une fenêtre SDL, dessine une série de points rouges
 * sur la diagonale et gère les événements jusqu'à ce que la fenêtre soit fermée.
 * 
 * @return EXIT_SUCCESS si le programme s'est exécuté avec succès.
 */
int test(Matrix& M, const char *sz);

/**
 * @brief Fonction de test pour le dessin 2D.
 * 
 * Cette fonction initialise une fenêtre SDL, dessine une série de points rouges
 * sur la diagonale et gère les événements jusqu'à ce que la fenêtre soit fermée.
 * 
 * @return EXIT_SUCCESS si le programme s'est exécuté avec succès.
 */

int test2D(Matrix& U, int nEsp, const char *sz);

/**
 * @brief Dessine un cercle.
 * 
 * Cette fonction dessine un cercle avec le renderer SDL donné.
 * 
 * @param renderer Le renderer SDL où dessiner le cercle.
 * @param centreX La coordonnée X du centre du cercle.
 * @param centreY La coordonnée Y du centre du cercle.
 * @param radius Le rayon du cercle à dessiner.
 */
void DrawCircle(SDL_Renderer *renderer, int32_t centreX, int32_t centreY, int32_t radius);

#endif /* GRAPHIQUE_H_INCLUDED */
