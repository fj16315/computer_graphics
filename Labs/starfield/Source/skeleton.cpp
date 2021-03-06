#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false


/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
int t;
vector<vec3> stars(1000);
//vector<vec3> newstars(1000);

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(screen* screen);
void Interpolate(float a, float b, vector<float>& result);
void Interpolate(vec3 a, vec3 b, vector<vec3>& result);

int main( int argc, char* argv[] )
{

  /*vector<float> result(10);
  Interpolate(5, 14, result);
  for (size_t i = 0; i < result.size(); i++) {
    cout << result[i] << " ";
  }

  vector<vec3> resultvec(4);
  vec3 a(1, 4, 9.2);
  vec3 b(4, 1, 9.8);
  Interpolate(a, b, resultvec);
  for (size_t i = 0; i < resultvec.size(); i++) {
    cout << "( " << resultvec[i].x << ", " << resultvec[i].y << ", " << resultvec[i].z << " )";
  }*/

  //vector<vec3> stars(1000);
  for (int i = 0; i < 1000; i++) {
    float r = float(rand())/float(RAND_MAX);
    stars[i].x = (2*r)-1;
    r = float(rand())/float(RAND_MAX);
    stars[i].y = (2*r)-1;
    r = float(rand())/float(RAND_MAX);
    stars[i].z = r;
  }

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  t = SDL_GetTicks();	/*Set start value for timer.*/

  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
  while( NoQuitMessageSDL() )
    {
      Draw(screen);
      Update();
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  //vec3 colour(1, 1, 1);
  vec3 colour;
  /*for(int i=0; i<screen->height; i++)
    {
      //uint32_t x = rand() % screen->width;
      //uint32_t y = rand() % screen->height;
      //uint32_t x = i % screen->width;
      uint32_t y = i;
      for (int j = 0; j < screen->width; j++) {
        uint32_t x = j;
        PutPixelSDL(screen, x, y, colour);
      }
      //PutPixelSDL(screen, x, y, colour);
    }*/
  /*vec3 topLeft(1, 0, 0);
  vec3 topRight(0, 0, 1);
  vec3 bottomRight(0, 1, 0);
  vec3 bottomLeft(1, 1, 0);

  vector<vec3> leftSide(SCREEN_HEIGHT);
  vector<vec3> rightSide(SCREEN_HEIGHT);
  vector<vec3> row(SCREEN_WIDTH);

  Interpolate(topLeft, bottomLeft, leftSide);
  Interpolate(topRight, bottomRight, rightSide);

  for (int i = 0; i < screen->height; i++) {
    Interpolate(leftSide[i], rightSide[i], row);
    for (size_t j = 0; j < row.size(); j++) {
      PutPixelSDL(screen, j, i, row[j]);
    }
  }*/

  for (size_t i = 0; i < stars.size(); i++) {
    uint32_t u = ((screen->height*stars[i].x)/(2*stars[i].z)) + (screen->width / 2);
    uint32_t v = ((screen->height*stars[i].y)/(2*stars[i].z)) + (screen->width / 2);
    colour = 0.2f * vec3(1, 1, 1) / (stars[i].z*stars[i].z);
    PutPixelSDL(screen, u, v, colour);
  }

}

/*Place updates of parameters here*/
void Update()
{
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  /*Good idea to remove this*/
  //std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/

  for (size_t i = 0; i < stars.size(); i++) {
    stars[i].z -= dt/32;
    if (stars[i].z <= 0) {
      stars[i].z += 1;
    }
    if (stars[i].z > 1) {
      stars[i].z -=1;
    }
  }
}

void Interpolate(float a, float b, vector<float>& result) {
  result[0] = a;
  if (result.size() != 1) {
    float incr = (b-a)/(result.size()-1);
    for (size_t i = 0; i < (result.size()-1); i++) {
      result[i+1] = result[i]+incr;
    }
  }
}

void Interpolate(vec3 a, vec3 b, vector<vec3>& result) {
  vector<float> xvalues(result.size());
  vector<float> yvalues(result.size());
  vector<float> zvalues(result.size());
  Interpolate(a.x, b.x, xvalues);
  Interpolate(a.y, b.y, yvalues);
  Interpolate(a.z, b.z, zvalues);
  for (size_t i = 0; i < result.size(); i++) {
    result[i].x = xvalues[i];
    result[i].y = yvalues[i];
    result[i].z = zvalues[i];
  }
}
