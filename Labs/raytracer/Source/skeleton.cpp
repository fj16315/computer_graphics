#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <math.h>
#include "limits.h"
#include <vector>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;

#define SCREEN_WIDTH 1920
#define SCREEN_HEIGHT 1080
#define FULLSCREEN_MODE false
#define PI 3.14159265
#define TRIANGLE 1
#define SPHERE 2
#define DIFFUSE 3
#define REFLECT 4
#define REFRACT 5

vector<Triangle> triangles;
vec3 defaultCameraPos = vec3(0, 0, -2.001);
mat3 R;
float yaw = 0;
vec3 lightPos(0, -0.5, -0.7);
vec3 lightColour = 14.f * vec3(1, 1, 1);
vec3 indirectLight = 0.5f*vec3(1, 1, 1);
int maxDepth = 5;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

struct Intersection
{
  vec3 position;
  float distance;
  int surfaceIndex;
  int type;
  vec3 colour;
};

struct Sphere
{
  vec3 centre;
  float radius;
  vec3 colour;
  vec3 normal;
  int type;
};

struct Surface
{
  int shapeIdentifier;
  vec3 v0;
  vec3 v1;
  vec3 v2;
  vec3 centre;
  float radius;
  vec3 colour;
  vec3 normal;
  int type;
};

vector<Sphere> spheres;
vector<Surface> surfaces;

bool Update();
void Draw(screen* screen);
bool ClosestIntersection(vec3 start, vec3 dir, const vector<Surface>& surfaces, Intersection& closestIntersection);
vec3 DirectLight(const Intersection& i);
void LoadSurfaces();
void LoadSpheres();
bool ComputeReflection(vec3 start, vec3 dir, const vector<Surface>& surfaces, Intersection& closestIntersection, int depth);
vec3 ComputeRefraction(vec3 I, vec3 N, float ior);
float Fresnel(vec3 I, vec3 N, float ior);

int main( int argc, char* argv[] )
{
  LoadSurfaces();

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  while ( Update())
    {
      Draw(screen);
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

  float focalLength = screen->height/2;
  vec3 cameraPos = defaultCameraPos*R;
  for (int i = 0; i < screen->height; i++)
  {
    for (int j = 0; j < screen->width; j++)
    {
      Intersection intersection;
      Intersection intersectionLight;
      vec3 dir = vec3(j-(screen->width/2), i-(screen->height/2), focalLength);
      dir = dir*R;
      if (ComputeReflection(cameraPos, dir, surfaces, intersection, 1))
      {
        //PutPixelSDL(screen, j, i, surfaces[intersection.surfaceIndex].colour);
        //PutPixelSDL(screen, j, i, DirectLight(intersection));
        //PutPixelSDL(screen, j, i, triangles[intersection.triangleIndex].color*DirectLight(intersection));
        PutPixelSDL(screen, j, i, intersection.colour*(DirectLight(intersection)+indirectLight));
      }
      else
      {
        PutPixelSDL(screen, j, i, vec3(0, 0, 0));
      }
    }
  }
}

/*Place updates of parameters here*/
bool Update()
{
  float radians;

  //static int t = SDL_GetTicks();
  /* Compute frame time */
  //int t2 = SDL_GetTicks();
  //float dt = float(t2-t);
  //t = t2;

  SDL_Event e;
  while(SDL_PollEvent(&e))
    {
      if (e.type == SDL_QUIT)
	{
	  return false;
	}
      else
	if (e.type == SDL_KEYDOWN)
	  {
	    int key_code = e.key.keysym.sym;
	    switch(key_code)
	      {
	      case SDLK_UP:
		/* Move camera forward */
		defaultCameraPos.z += 0.5;
		break;
	      case SDLK_DOWN:
		/* Move camera backwards */
		defaultCameraPos.z -= 0.5;
		break;
	      case SDLK_LEFT:
		/* Move camera left */
		yaw -= 10;
		yaw = fmod(yaw, 360);
		radians = (yaw*PI)/180;
		R = mat3(vec3(cos(radians), 0, -sin(radians)), vec3(0, 1, 0), vec3(sin(radians), 0, cos(radians)));
		break;
	      case SDLK_RIGHT:
		/* Move camera right */
		yaw += 10;
		yaw = fmod(yaw, 360);
		radians = (yaw*PI)/180;
		R = mat3(vec3(cos(radians), 0, -sin(radians)), vec3(0, 1, 0), vec3(sin(radians), 0, cos(radians)));
		break;
	      case SDLK_ESCAPE:
		/* Move camera quit */
		return false;
		break;
              case SDLK_w:
		/* Move light forward */
		lightPos.z += 0.1;
		break;
              case SDLK_s:
		/* Move light backward */
		lightPos.z -= 0.1;
		break;
              case SDLK_a:
		/* Move light left */
		lightPos.x -= 0.1;
		break;
              case SDLK_d:
		/* Move light right */
		lightPos.x += 0.1;
		break;
              case SDLK_q:
		/* Move light down */
		lightPos.y += 0.1;
		break;
              case SDLK_e:
		/* Move light up */
		lightPos.y -= 0.1;
		break;
	      }
	  }
    }
  return true;
}

bool ClosestIntersection(vec3 start, vec3 dir, const vector<Surface>& surfaces, Intersection& closestIntersection)
{
  float furthest = std::numeric_limits<float>::max();
  closestIntersection.distance = furthest;

  for (uint i = 0; i < surfaces.size(); i++)
  {
    if (surfaces[i].shapeIdentifier == TRIANGLE) {
      vec3 v0 = surfaces[i].v0;
      vec3 v1 = surfaces[i].v1;
      vec3 v2 = surfaces[i].v2;

      vec3 e1 = v1-v0;
      vec3 e2 = v2-v0;
      vec3 b = start-v0;

      mat3 A(-dir, e1, e2);
      vec3 x = inverse(A)*b;

      if ((x.y+x.z<1) && (x.x>0) && (x.y>=0) && (x.z>=0))
      {
	vec3 r = start+(x.x*dir);
	float distance = length(x.x*dir);

        if (distance < closestIntersection.distance)
        {
          closestIntersection.distance = distance;
          closestIntersection.position = r;
          closestIntersection.surfaceIndex = i;
	  closestIntersection.colour = surfaces[i].colour;
        }
      }
    }

    else if (surfaces[i].shapeIdentifier == SPHERE) {
      vec3 centre = surfaces[i].centre;
      float radius = surfaces[i].radius;
      vec3 L = centre-start;
      float tca = dot(normalize(dir), L);

      if (tca >= 0) {
	float d = sqrt(dot(L, L)-pow(tca, 2));

	if ((d >= 0) && (d <= radius)) {
          float thc = sqrt(pow(radius, 2)-pow(d, 2));
          float t0 = tca-thc;
          vec3 P = start+(t0*normalize(dir));
	  float distance = length(P-start);

	  if (distance < closestIntersection.distance) {
	    closestIntersection.distance = distance;
	    closestIntersection.position = P;
	    closestIntersection.surfaceIndex = i;
	    closestIntersection.colour = surfaces[i].colour;
	  }
	}
      }
    }
    if (i == surfaces.size()-1)
    {
      if (closestIntersection.distance == furthest)
      {
        return false;
      }
      else
      {
        return true;
      }
    }
  }
  return false;
}

vec3 DirectLight(const Intersection& i)
{
  vec3 ray = lightPos-i.position;
  Intersection intersectionLight;
  vec3 normal = surfaces[i.surfaceIndex].normal;

  if (surfaces[i.surfaceIndex].shapeIdentifier == SPHERE) {
    normal = i.position-surfaces[i.surfaceIndex].centre;
  }

  normal = normalize(normal);

  vec3 start = i.position+vec3(normal.x*0.0001, normal.y*0.0001, normal.z*0.0001);
  ClosestIntersection(start, ray, surfaces, intersectionLight);
  float r = length(ray);
  float scalar = dot(normalize(ray), normal);

  if (scalar < 0)
  {
    scalar = 0;
  }

  vec3 directLight;
  if ((acos(scalar)*(180/PI) > 90) || (intersectionLight.distance < r))
  {
    directLight = vec3(0, 0, 0);
  }
  else
  {
    float coef = scalar/(4*PI*pow(r, 2));
    directLight = lightColour*coef;
  }
  return directLight;
}

void LoadSurfaces() {
  Surface surface;
  LoadTestModel(triangles);
  LoadSpheres();

  for (uint i = 0; i < /*triangles.size()*/30; i++) {
    surface.shapeIdentifier = TRIANGLE;
    surface.v0 = vec3(triangles[i].v0.x, triangles[i].v0.y, triangles[i].v0.z);
    surface.v1 = vec3(triangles[i].v1.x, triangles[i].v1.y, triangles[i].v1.z);
    surface.v2 = vec3(triangles[i].v2.x, triangles[i].v2.y, triangles[i].v2.z);
    surface.colour = triangles[i].color;
    //surface.type = triangles[i].type;
    surface.type = DIFFUSE;
    surface.normal = vec3(triangles[i].normal.x, triangles[i].normal.y, triangles[i].normal.z);
    surfaces.push_back(surface);
  }

  for (uint i = 0; i < spheres.size(); i++) {
    surface.shapeIdentifier = SPHERE;
    surface.centre = spheres[i].centre;
    surface.radius = spheres[i].radius;
    surface.colour = spheres[i].colour;
    surface.normal = spheres[i].normal;
    surface.type = spheres[i].type;
    surfaces.push_back(surface);
  }
}

void LoadSpheres() {
  Sphere sphere;
  sphere.centre = vec3(-0.5, 0.6, -0.5);
  sphere.radius = 0.4;
  sphere.colour = vec3(1, 0, 0);
  sphere.type = REFLECT;
  spheres.push_back(sphere);
}

bool ComputeReflection(vec3 start, vec3 dir, const vector<Surface>& surfaces, Intersection& closestIntersection, int depth) {
  if (depth > maxDepth) {
    return false;
  }

  if (ClosestIntersection(start, dir, surfaces, closestIntersection)) {
    Surface currentSurface = surfaces[closestIntersection.surfaceIndex];
    if (currentSurface.type == REFLECT) {
      vec3 indicentRay = closestIntersection.position-start;
      vec3 normal = currentSurface.normal;

      if (currentSurface.shapeIdentifier == SPHERE) {
	normal = closestIntersection.position-currentSurface.centre;
      }

      float scalar = dot(indicentRay, normalize(normal));
      vec3 coef = 2*scalar*normalize(normal);
      vec3 reflectedRay = indicentRay-coef;

      return ComputeReflection(closestIntersection.position, reflectedRay, surfaces, closestIntersection, depth+1);
    }
    else if (currentSurface.type == DIFFUSE) {
      return true;
    }
    /*else if (currentSurface.type == REFRACT) {
      vec3 normal = currentSurface.normal;

      if (currentSurface.shapeIdentifier == SPHERE) {
	normal = closestIntersection.position-currentSurface.centre;
      }
      normal = normalize(normal);

      float kr = Fresnel(dir, normal, 1.5);
      bool outside = 0;
      if (dot(dir, normal) < 0) {
	outside = 1;
      }

      Intersection intersectionRefract;
      if (kr < 1) {
	vec3 refractDir = ComputeRefraction(dir, normal, 1.5);
	vec3 newStart;
	if (outside) {
	  newStart = start-vec3(0.001*normal.x, 0.001*normal.y, 0.001*normal.z);
	}
	else {
	  newStart = start+vec3(0.001*normal.x, 0.001*normal.y, 0.001*normal.z);
	}
	ComputeReflection(newStart, refractDir, surfaces, intersectionRefract, depth+1);
      }

      float scalar = dot(dir, normalize(normal));
      vec3 coef = 2*scalar*normalize(normal);
      vec3 reflectedRay = dir-coef;
      vec3 newStart2;
      if (outside) {
	newStart2 = start+vec3(0.001*normal.x, 0.001*normal.y, 0.001*normal.z);
      }
      else {
	newStart2 = start-vec3(0.001*normal.x, 0.001*normal.y, 0.001*normal.z);
      }

      Intersection intersectionReflect;
      ComputeReflection(newStart2, reflectedRay, surfaces, intersectionReflect, depth+1);

      vec3 reflectColour = vec3(intersectionReflect.colour.x*kr, intersectionReflect.colour.y*kr, intersectionReflect.colour.z*kr);
      vec3 refractColour = vec3(intersectionRefract.colour.x*(1-kr), intersectionRefract.colour.y*(1-kr), intersectionRefract.colour.z*(1-kr));
      closestIntersection.colour = reflectColour+refractColour;

      vec3 refractedRay = ComputeRefraction(dir, normal, 1.5);
      return ComputeReflection(closestIntersection.position, refractedRay, surfaces, closestIntersection, depth+1);
      //return true;
    }*/
  }
  return false;
}

float Fresnel(vec3 I, vec3 N, float ior) {
  float kr;
  float cosi = dot(I, N);
  float etai = 1;
  float etat = ior;
  if (cosi > 0) {
    swap(etai, etat);
  }
  float sint = etai/(etat*sqrt(max(0.f, 1-pow(cosi, 2))));
  if (sint >= 1) {
    kr = 1;
  }
  else {
    float cost = sqrt(max(0.f, 1-pow(sint, 2)));
    cosi = fabs(cosi);
    float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
    float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
    kr = (Rs * Rs + Rp * Rp) / 2;
  }
  return kr;
}

vec3 ComputeRefraction(vec3 I, vec3 N, float ior) {
  float cosi = dot(I, N);
  float etai = 1;
  float etat = ior;
  if (cosi < 0) {
    cosi = -cosi;
  }
  else {
    swap(etai, etat);
    N = -N;
  }
  float eta = etai/etat;
  float k = 1-(pow(eta, 2)*(1-pow(cosi, 2)));
  vec3 dir = (eta*I)+(((eta*cosi)-sqrt(k))*N);
  return dir;
}