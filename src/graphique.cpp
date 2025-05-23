#include "graphique.h"
#include <iostream>
#include <unistd.h>    

#define WINDOW_WIDTH 500
#define MARGIN 30


//int SDL_Line(int ix,int iy,fx,fy)

void sdlLine(SDL_Renderer *r, int xi, int yi, int xf, int yf)
{
    int x,y,ix,iy,z;
    double de,e=0.0;

    if (xi>xf)
    {
        z = xf;
        xf = xi;
        xi = z;
    }

   if (yi>yf)
    {
        z = yf;
        yf = yi;
        yi = z;
    }

    double dx = xf-xi;
    double dy = yf-yi;

    if(dx>dy)
    {
        ix = 1;
        iy = 0;
        de = dy/dx;
    }
    else
    {
        ix = 0;
        iy = 1;
        de = dx/dy;
    }


    x = xi;
    y = yi;

    while(x<xf)
    {
        SDL_RenderDrawPoint(r, x, y);
        if (e>1) {
            x+= iy;
            y+= ix;
            e--;
        }
        else {
            x+= ix;
            y+= iy;        
            e+= de;
        }
    }
}

int test(Matrix& U, const char *sz) {
    SDL_Event event;
    SDL_Renderer *renderer;
    SDL_Window *window;
    int i,j,aj,ak,bj,bk;

    int n = U.getN();
    int size = WINDOW_WIDTH - 2*MARGIN;

    double Mn = U.min();
    double Mx = U.max();

    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(WINDOW_WIDTH, WINDOW_WIDTH, 0, &window,
&renderer);


    for(i=0;i<n;i++){
        
        std::cout << sz <<" --> temps: "<< i << std::endl;
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 0);
        SDL_RenderClear(renderer);

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);

        bj = MARGIN;
        bk = MARGIN + U.getGraphical(i,0,Mx,Mn,size);

        for (j = 1; j < size; ++j)  {
                aj = MARGIN + j;
                ak = MARGIN + U.getGraphical(i,j,Mx,Mn,size);
                sdlLine(renderer, bj, bk, aj, ak);
                bj = aj;
                bk = ak;
        }

         //    SDL_RenderDrawPoint(renderer, MARGIN + j, MARGIN + U.getGraphical(i,j,Mx,Mn,size));
    
        sdlLine(renderer, MARGIN, MARGIN+size, MARGIN+size, MARGIN+size);
        sdlLine(renderer, MARGIN, MARGIN, MARGIN+1, MARGIN+size);

        SDL_RenderPresent(renderer);
 
        usleep(100000); //! 100ms de pause 
 
     }
std::cout << "Fermer la fenetre graphique pour continuer" << std::endl;
while (1) {
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT)
            break;
    }
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return EXIT_SUCCESS;
}

int test2D(Matrix& U, int nEsp, const char *sz) {
    SDL_Event event;
    SDL_Renderer *renderer;
    SDL_Window *window;
    int i,j,k;
    int R,G,B;

    int n = U.getN();
    int size = WINDOW_WIDTH - 2*MARGIN;

    double Mn = U.min();
    double Mx = U.max();

    std::cout << Mn << "<" << Mx << std::endl;

    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(WINDOW_WIDTH, WINDOW_WIDTH, 0, &window,
&renderer);


    for(i=0;i<n;i++){
        std::cout << sz <<" --> temps: "<< i << std::endl;
        
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
        SDL_RenderClear(renderer);

        for (j = 0; j < size; j++)
            for (k = 0; k < size; k++) {
                 U.getColor(i,j,k,Mx,Mn,nEsp,size,R,G,B);

                 SDL_SetRenderDrawColor(renderer, R, G, B, 255);
             SDL_RenderDrawPoint(renderer, MARGIN + j, MARGIN + k);
            }
        SDL_RenderPresent(renderer);
 
        usleep(1000); //! 100ms de pause 
 
     }

std::cout << "Fermer la fenetre graphique pour continuer" << std::endl;

while (1) {
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT)
            break;
    }
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    return EXIT_SUCCESS;
}


void DrawCircle(SDL_Renderer * renderer, int32_t centreX, int32_t
centreY, int32_t radius)
{
    const int32_t diameter = (radius * 2);

    int32_t x = (radius - 1);
    int32_t y = 0;
    int32_t tx = 1;
    int32_t ty = 1;
    int32_t error = (tx - diameter);

    while (x >= y)
    {
      //  Each of the following renders an octant of the circle
      SDL_RenderDrawPoint(renderer, centreX + x, centreY - y);
      SDL_RenderDrawPoint(renderer, centreX + x, centreY + y);
      SDL_RenderDrawPoint(renderer, centreX - x, centreY - y);
      SDL_RenderDrawPoint(renderer, centreX - x, centreY + y);
      SDL_RenderDrawPoint(renderer, centreX + y, centreY - x);
      SDL_RenderDrawPoint(renderer, centreX + y, centreY + x);
      SDL_RenderDrawPoint(renderer, centreX - y, centreY - x);
      SDL_RenderDrawPoint(renderer, centreX - y, centreY + x);

      if (error <= 0)
      {
          ++y;
          error += ty;
          ty += 2;
      }

      if (error > 0)
      {
          --x;
          tx += 2;
          error += (tx - diameter);
      }
    }
}