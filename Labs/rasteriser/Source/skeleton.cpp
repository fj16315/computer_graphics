#include <iostream>
#include <glm/glm.hpp>
#include "SDL2/SDL.h"
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <math.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::ivec2;
using glm::vec2;

SDL_Event event;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false
#define PI 3.14159265

vector<Triangle> triangles;
float yaw = 0;
vec4 cameraPos(0, 0, 3, 1);
mat4 R;
vec3 currentColour = vec3(1, 1, 1);
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec4 lightPos(0, -0.5, 2.3, 1);
vec3 lightPower = 14.f*vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f*vec3(1, 1, 1);
vec4 currentNormal;
vec3 currentReflectance;
vec3 colourBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
int stencilBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
float n = 0.001;
float f = 5;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

struct Pixel {
  int x;
  int y;
  float zinv;
  vec4 pos3d;
};

struct Vertex {
  vec4 position;
};

struct Surface {
  vec4 normal;
  vector<Vertex> vertices;
  bool scalar;
};

bool Update();
void Draw(screen* screen);
void TransformationMatrix(mat4& M);
void DrawPolygon(screen* screen, const vector<Vertex>& vertices);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void VertexShader(const Vertex& v, Pixel& p);
void PixelShader(screen* screen, const Pixel& p);
void ComputeShadows(const vector<Vertex>& vertices, screen* screen);
void DrawShadows(screen* screen);
void ClipAtWEqualsZero(const vector<Vertex>& vertices, vector<Vertex>& in);
void ClipAtWEqualsPlane(const vector<Vertex>& vertices, int plane, vector<Vertex>& in);
void Clip(const vector<Vertex>& vertices, vector<Vertex>& in);

int main( int argc, char* argv[] )
{
  LoadTestModel(triangles);
  cout << "Hello World!" << endl;

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
  for (int y = 0; y < SCREEN_HEIGHT; y++) {
    for (int x = 0; x < SCREEN_WIDTH; x++) {
      depthBuffer[y][x] = 0;
      colourBuffer[y][x] = vec3(0, 0, 0);
      stencilBuffer[y][x] = 0;
    }
  }

  yaw = fmod(yaw, 360);
  float radians = (yaw*PI)/180;
  mat4 rotation = mat4(vec4(cos(radians), 0, -sin(radians), 0), vec4(0, 1, 0, 0), vec4(sin(radians), 0, cos(radians), 0), vec4(0, 0, 0, 1));

  /*triangles[0].v0 = vec4(-0.75, 0.75, -0.5, 1);
  triangles[0].v1 = vec4(0.75, 0.75, -0.5, 1);
  triangles[0].v2 = vec4(0, 0.75, 0.5, 1);
  vec4 line1 = triangles[0].v1-triangles[0].v0;
  vec3 line13 = vec3(line1.x, line1.y, line1.z);
  vec4 line2 = triangles[0].v2-triangles[0].v0;
  vec3 line23 = vec3(line2.x, line2.y, line2.z);
  vec3 cross = glm::cross(line13, line23);
  triangles[0].normal = normalize(vec4(cross.x, cross.y, cross.z, 1));

  triangles[1].v0 = vec4(-0.25, 0.25, -0.25, 1);
  triangles[1].v1 = vec4(0.25, 0.25, -0.25, 1);
  triangles[1].v2 = vec4(0, 0.25, 0.25, 1);
  line1 = triangles[1].v1-triangles[1].v0;
  line13 = vec3(line1.x, line1.y, line1.z);
  line2 = triangles[1].v2-triangles[1].v0;
  line23 = vec3(line2.x, line2.y, line2.z);
  cross = glm::cross(line13, line23);
  triangles[0].normal = normalize(vec4(cross.x, cross.y, cross.z, 1));*/

  mat4 M;
  TransformationMatrix(M);
  mat4 projectionMatrix = mat4(vec4(0.8, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, (f+n)/(f-n), 1), vec4(0, 0, -((2*f*n)/(f-n)), 0));

  for (uint i = 0; i < triangles.size(); i++) {
    currentColour = triangles[i].color;
    vector<Vertex> vertices(3);

    vertices[0].position = M*triangles[i].v0;
    vertices[0].position = projectionMatrix*vertices[0].position;
    vertices[1].position = M*triangles[i].v1;
    vertices[1].position = projectionMatrix*vertices[1].position;
    vertices[2].position = M*triangles[i].v2;
    vertices[2].position = projectionMatrix*vertices[2].position;

    currentNormal = rotation*triangles[i].normal;
    currentReflectance = triangles[i].color;

    vector<Vertex> in(0);
    Clip(vertices, in);
    if (in.size() != 0) {
      DrawPolygon(screen, in);
    }
  }

  /*for (uint i = 0; i < triangles.size(); i++) {
    vector<Vertex> vertices(3);

    vertices[0].position = M*triangles[i].v0;
    vertices[1].position = M*triangles[i].v1;
    vertices[2].position = M*triangles[i].v2;

    currentNormal = rotation*triangles[i].normal;

    //cout << "computing shadows" << endl;
    //if ((dot(vertices[0].position, currentNormal) <= 0) || (dot(vertices[1].position, currentNormal) <= 0) || (dot(vertices[2].position, currentNormal) <= 0)) {
      //ComputeShadows(vertices, screen);
    //}
  }*/

  //DrawShadows(screen);
}

/*Place updates of parameters here*/
bool Update()
{
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
    cameraPos.z -= 0.5;
    lightPos.z -= 0.5;
		break;
	      case SDLK_DOWN:
		/* Move camera backwards */
    cameraPos.z += 0.5;
    lightPos.z += 0.5;
		break;
	      case SDLK_LEFT:
		/* Move camera left */
    yaw -= 10;
		break;
	      case SDLK_RIGHT:
		/* Move camera right */
    yaw += 10;
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
	      case SDLK_ESCAPE:
		/* Move camera quit */
		return false;
	      }
	  }
    }
  return true;
}

void TransformationMatrix(mat4& M) {
  yaw = fmod(yaw, 360);
  float radians = (yaw*PI)/180;

  M = mat4(vec4(cos(radians), 0, -sin(radians), 0), vec4(0, 1, 0, 0), vec4(sin(radians), 0, cos(radians), 0), cameraPos);
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result) {
  int n = result.size();
  int currentX = a.x;
  int currentY = a.y;
  int dx = b.x-a.x;
  int dy = b.y-a.y;
  int flippedX = 1;
  int flippedY = 1;
  int transfer = 0;
  int d = 0;
  int twodxdy = 0;
  int twodydx = 0;
  if (dx < 0) {
    flippedX = -1;
    transfer = a.x;
    a.x = b.x;
    b.x = transfer;
    dx = -dx;
  }
  if (dy < 0) {
    flippedY = -1;
    transfer = a.y;
    a.y = b.y;
    b.y = transfer;
    dy = -dy;
  }
  if (dx > dy) {
    twodydx = 2*(dy-dx);
    d = (2*dy)-dx;
  } else {
    twodxdy = 2*(dx-dy);
    d = (2*dx)-dy;
  }
  float stepZinv = float(b.zinv-a.zinv)/float(max(n-1, 1));
  vec4 stepPos = (b.pos3d-a.pos3d)/float(max(n-1, 1));
  float currentZinv = a.zinv;
  vec4 currentPos = a.pos3d;

  for (int i = 0; i < n; i++) {
    result[i].x = currentX;
    result[i].y = currentY;
    if (dy == 0) {
      currentX += flippedX;
    } else if (dx == 0) {
      currentY += flippedY;
    } else if (dx > dy) {
      currentX += flippedX;
      if (d < 0) {
	d += 2*dy;
      } else {
	currentY += flippedY;
	d += twodydx;
      }
    } else {
      currentY += flippedY;
      if (d < 0) {
	d += 2*dx;
      } else {
	currentX += flippedX;
	d += twodxdy;
      }
    }
    result[i].zinv = currentZinv;
    result[i].pos3d = currentPos;
    currentZinv += stepZinv;
    currentPos += stepPos;
  }
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
  int maxY = vertexPixels[0].y;
  int minY = vertexPixels[0].y;

  for (uint i = 1; i < vertexPixels.size(); i++) {
    maxY = max(maxY, vertexPixels[i].y);
    minY = min(minY, vertexPixels[i].y);
  }

  int rows = maxY-minY+1;
  leftPixels.resize(rows);
  rightPixels.resize(rows);

  for (int i = 0; i < rows; i++) {
    leftPixels[i].x = numeric_limits<int>::max();
    rightPixels[i].x = -numeric_limits<int>::max();
  }

  for (uint i = 0; i < vertexPixels.size(); i++) {
    Pixel delta(vertexPixels[0]);
    delta.x = glm::abs(vertexPixels[i].x-vertexPixels[(i+1)%vertexPixels.size()].x);
    delta.y = glm::abs(vertexPixels[i].y-vertexPixels[(i+1)%vertexPixels.size()].y);
    int pixels = glm::max(delta.x, delta.y)+1;
    vector<Pixel> line(pixels);
    Interpolate(vertexPixels[i], vertexPixels[(i+1)%vertexPixels.size()], line);

    for (int j = 0; j < pixels; j++) {
      int diff = line[j].y - minY;
      if (line[j].x <= leftPixels[diff].x) {
        leftPixels[diff] = line[j];
      }

      if (line[j].x >= rightPixels[diff].x) {
        rightPixels[diff] = line[j];
      }
    }
  }
}

void DrawPolygonRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels) {
  for (uint i = 0; i < leftPixels.size(); i++) {
    vector<Pixel> row(rightPixels[i].x-leftPixels[i].x+1);
    Interpolate(leftPixels[i], rightPixels[i], row);
    for (uint j = 0; j < row.size(); j++) {
      PixelShader(screen, row[j]);
    }
  }
}

void VertexShader(const Vertex& v, Pixel& p) {
  vec4 perspective = v.position/v.position.w;

  mat4 viewPort = mat4(vec4((SCREEN_WIDTH-1)/2, 0, 0, 0), vec4(0, (SCREEN_HEIGHT-1)/2, 0, 0), vec4(0, 0, 0.5, 0), vec4((SCREEN_WIDTH-1)/2, (SCREEN_HEIGHT-1)/2, 0.5, 1));
  perspective = viewPort*perspective;

  p.zinv = 1/v.position.w;
  p.x = perspective.x;
  p.y = perspective.y;
  p.pos3d.x = v.position.x/v.position.w;
  p.pos3d.y = v.position.y/v.position.w;
  p.pos3d.z = 1;
  p.pos3d.w = 1;
}

void DrawPolygon(screen* screen, const vector<Vertex>& vertices) {
  int v = vertices.size();
  vector<Pixel> projectedVertices(v);

  for (int i = 0; i < v; i++) {
    VertexShader(vertices[i], projectedVertices[i]);
  }

  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows(projectedVertices, leftPixels, rightPixels);
  DrawPolygonRows(screen, leftPixels, rightPixels);
}

void PixelShader(screen* screen, const Pixel& p) {
  if (p.zinv > depthBuffer[p.y][p.x]) {
    vec4 position = vec4((p.pos3d.x*1.25)/p.zinv, p.pos3d.y/p.zinv, 1/p.zinv, 1);
    vec4 pixelToLight = lightPos-position;
    float scalar = dot(normalize(pixelToLight), currentNormal);
    if (scalar < 0) {
      scalar = 0;
    }
    vec3 top = lightPower*scalar;
    float bottom = 4*PI*pow(length(pixelToLight), 2);
    vec3 D = top/bottom;
    vec3 R = currentReflectance*(D+indirectLightPowerPerArea);
    depthBuffer[p.y][p.x] = p.zinv;
    PutPixelSDL(screen, p.x, p.y, R);
    //colourBuffer[p.y][p.x] = R;
  }
}

void ComputeShadows(const vector<Vertex>& vertices, screen* screen) {
  int n = vertices.size();
  bool alreadyTested[SCREEN_HEIGHT][SCREEN_WIDTH];
  for (int y = 0; y < SCREEN_HEIGHT; y++) {
    for (int x = 0; x < SCREEN_WIDTH; x++) {
      alreadyTested[y][x] = 0;
    }
  }
  //mat4 M;
  //TransformationMatrix(M);
  vector<Surface> surfaces(n+2);
  vector<vec4> rayVertices(0);
  //int max = numeric_limits<float>::infinity();
  int max = 100;
  //mat4 P = mat4(vec4(0.8, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 200/199, 1), vec4(0, 0, -(400/199), 0));

  for (int i = 0; i < n; i++) {
    surfaces[i].vertices.resize(4);
    vec4 currentVertex = vertices[i].position;
    vec4 nextVertex = vertices[(i+1)%n].position;
    vec4 ray1 = currentVertex-lightPos;
    vec4 rayExtend1 = vec4(max*ray1.x, max*ray1.y, max*ray1.z, 1);
    surfaces[i].vertices[3].position = currentVertex+rayExtend1;

    if (surfaces[i].vertices[3].position.z < 1) {
      float t = (1-currentVertex.z)/(surfaces[i].vertices[3].position.z-currentVertex.z);
      //cout << "t: " << t << endl;
      vec4 I = currentVertex+(t*(surfaces[i].vertices[3].position-currentVertex));
      surfaces[i].vertices[3].position = I;
    }

    if (surfaces[i].vertices[3].position.z > 20) {
      float t = (20-currentVertex.z)/(surfaces[i].vertices[3].position.z-currentVertex.z);
      //cout << "t: " << t << endl;
      vec4 I = currentVertex+(t*(surfaces[i].vertices[3].position-currentVertex));
      surfaces[i].vertices[3].position = I;
    }

    if (surfaces[i].vertices[3].position.x < -20) {
      float t = (-20-currentVertex.x)/(surfaces[i].vertices[3].position.x-currentVertex.x);
      //cout << "t: " << t << endl;
      vec4 I = currentVertex+(t*(surfaces[i].vertices[3].position-currentVertex));
      surfaces[i].vertices[3].position = I;
    }

    if (surfaces[i].vertices[3].position.x > 20) {
      float t = (20-currentVertex.x)/(surfaces[i].vertices[3].position.x-currentVertex.x);
      //cout << "t: " << t << endl;
      vec4 I = currentVertex+(t*(surfaces[i].vertices[3].position-currentVertex));
      surfaces[i].vertices[3].position = I;
    }

    if (surfaces[i].vertices[3].position.y < -20) {
      float t = (-20-currentVertex.y)/(surfaces[i].vertices[3].position.y-currentVertex.y);
      //cout << "t: " << t << endl;
      vec4 I = currentVertex+(t*(surfaces[i].vertices[3].position-currentVertex));
      surfaces[i].vertices[3].position = I;
    }

    if (surfaces[i].vertices[3].position.y > 20) {
      float t = (20-currentVertex.y)/(surfaces[i].vertices[3].position.y-currentVertex.y);
      //cout << "t: " << t << endl;
      vec4 I = currentVertex+(t*(surfaces[i].vertices[3].position-currentVertex));
      surfaces[i].vertices[3].position = I;
    }
    rayVertices.push_back(surfaces[i].vertices[3].position);

    vec4 ray2 = nextVertex-lightPos;
    vec4 rayExtend2 = vec4(max*ray2.x, max*ray2.y, max*ray2.z, 1);
    surfaces[i].vertices[2].position = nextVertex+rayExtend2;

    if (surfaces[i].vertices[2].position.z < 1) {
      float t = (1-nextVertex.z)/(surfaces[i].vertices[2].position.z-nextVertex.z);
      //cout << "t: " << t << endl;
      vec4 I = nextVertex+(t*(surfaces[i].vertices[2].position-nextVertex));
      surfaces[i].vertices[2].position = I;
    }

    if (surfaces[i].vertices[2].position.z > 20) {
      float t = (20-nextVertex.z)/(surfaces[i].vertices[2].position.z-nextVertex.z);
      //cout << "t: " << t << endl;
      vec4 I = nextVertex+(t*(surfaces[i].vertices[2].position-nextVertex));
      surfaces[i].vertices[2].position = I;
    }
    if (surfaces[i].vertices[2].position.x < -20) {
      float t = (-20-nextVertex.x)/(surfaces[i].vertices[2].position.x-nextVertex.x);
      //cout << "t: " << t << endl;
      vec4 I = nextVertex+(t*(surfaces[i].vertices[2].position-nextVertex));
      surfaces[i].vertices[2].position = I;
    }
    if (surfaces[i].vertices[2].position.x > 20) {
      float t = (20-nextVertex.x)/(surfaces[i].vertices[2].position.x-nextVertex.x);
      //cout << "t: " << t << endl;
      vec4 I = nextVertex+(t*(surfaces[i].vertices[2].position-nextVertex));
      surfaces[i].vertices[2].position = I;
    }
    if (surfaces[i].vertices[2].position.y < -20) {
      float t = (-20-nextVertex.y)/(surfaces[i].vertices[2].position.y-nextVertex.y);
      //cout << "t: " << t << endl;
      vec4 I = nextVertex+(t*(surfaces[i].vertices[2].position-nextVertex));
      surfaces[i].vertices[2].position = I;
    }
    if (surfaces[i].vertices[2].position.y > 20) {
      float t = (20-nextVertex.y)/(surfaces[i].vertices[2].position.y-nextVertex.y);
      //cout << "t: " << t << endl;
      vec4 I = nextVertex+(t*(surfaces[i].vertices[2].position-nextVertex));
      surfaces[i].vertices[2].position = I;
    }

    surfaces[i].vertices[0].position = currentVertex;
    surfaces[i].vertices[1].position = nextVertex;
    //cout << "vertices:" << endl;
    //cout << endl;
    //cout << "vertex[0]: (" << surfaces[i].vertices[0].position.x << "," << surfaces[i].vertices[0].position.y << "," << surfaces[i].vertices[0].position.z << ")" << endl;
    //cout << "vertex[1]: (" << surfaces[i].vertices[1].position.x << "," << surfaces[i].vertices[1].position.y << "," << surfaces[i].vertices[1].position.z << ")" << endl;
    //cout << "vertex[2]: (" << surfaces[i].vertices[2].position.x << "," << surfaces[i].vertices[2].position.y << "," << surfaces[i].vertices[2].position.z << ")" << endl;
    //cout << "vertex[3]: (" << surfaces[i].vertices[3].position.x << "," << surfaces[i].vertices[3].position.y << "," << surfaces[i].vertices[3].position.z << ")" << endl;
    vec4 line1 = nextVertex-currentVertex;
    vec4 line2 = rayExtend1;
    vec3 line13 = vec3(line1.x, line1.y, line1.z);
    vec3 line23 = vec3(line2.x, line2.y, line2.z);
    vec3 normal3 = cross(line13, line23);
    vec4 normal = normalize(vec4(normal3.x, normal3.y, normal3.z, 1));
    surfaces[i].scalar = (dot(normal, -currentVertex) <= 0)/* || (dot(normal, -surfaces[i].vertices[1].position) <= 0) || (dot(normal, -surfaces[i].vertices[2].position) <= 0) || (dot(normal, -surfaces[i].vertices[3].position) <= 0)*/;
    //cout << "scalar: " << surfaces[i].scalar << endl;
    //vector<Vertex> in(0);
    //Clip(surfaces[i].vertices, in);
    //surfaces[i].vertices = in;
  }

  surfaces[n].normal = currentNormal;
  //surfaces[n].scalar = (dot(currentNormal, -vertices[0].position) <= 0) || (dot(currentNormal, -vertices[1].position) <= 0) || (dot(currentNormal, -vertices[2].position) <= 0);
  surfaces[n].scalar = 1;
  surfaces[n].vertices.resize(3);
  for (int i = 0; i < 3; i++) {
    surfaces[n].vertices[i].position = vertices[i].position;
    //cout << "vertex[0]: (" << surfaces[i].vertices[0].position.x << "," << surfaces[i].vertices[0].position.y << "," << surfaces[i].vertices[v].position.z << ")" << endl;
  }

  surfaces[n+1].normal = currentNormal;
  //surfaces[n+1].scalar = (dot(currentNormal, -vertices[0].position) <= 0) || (dot(currentNormal, -vertices[1].position) <= 0) || (dot(currentNormal, -vertices[2].position) <= 0);
  surfaces[n+1].scalar = 1;
  surfaces[n+1].vertices.resize(3);
  for (int i = 0; i < 3; i++) {
    surfaces[n+1].vertices[i].position = rayVertices[i];
    //cout << "vertex[0]: (" << surfaces[i].vertices[0].position.x << "," << surfaces[i].vertices[0].position.y << "," << surfaces[i].vertices[v].position.z << ")" << endl;
  }

  for (uint i = 0; i < surfaces.size(); i++) {
    vector<Pixel> projectedVertices(surfaces[i].vertices.size());
    //cout << "scalar: " << surfaces[i].scalar << endl;
    if (!surfaces[i].scalar) {
      for (uint v = 0; v < surfaces[i].vertices.size(); v++) {
	//cout << "in: (" << surfaces[i].vertices[v].position.x*1.25 << "," << surfaces[i].vertices[v].position.y << "," << surfaces[i].vertices[v].position.w << ")" << endl;
	//cout << "before vertexShader" << endl;
	VertexShader(surfaces[i].vertices[v], projectedVertices[v]);
	//cout << "pixel[" << v << "]: (" << projectedVertices[v].x << "," << projectedVertices[v].y << ")" << endl;
      }
      vector<Pixel> leftPixels;
      vector<Pixel> rightPixels;
      ComputePolygonRows(projectedVertices, leftPixels, rightPixels);
      //cout << "computed polygon rows" << endl;
      for (uint j = 0; j < leftPixels.size(); j++) {
	vector<Pixel> row(rightPixels[j].x-leftPixels[j].x+1);
	Interpolate(leftPixels[j], rightPixels[j], row);
	for (uint k = 0; k < row.size(); k++) {
	  if ((row[k].x >= 0) && (row[k].x < SCREEN_WIDTH) && (row[k].y >= 0) && (row[k].y < SCREEN_HEIGHT)) {
	    if ((row[k].zinv <= depthBuffer[row[k].y][row[k].x]) && (!alreadyTested[row[k].y][row[k].x])) {
	      stencilBuffer[row[k].y][row[k].x]++;
	      alreadyTested[row[k].y][row[k].x] = 1;
	      //PutPixelSDL(screen, row[k].x, row[k].y, vec3(0, 0, 1));
	      //alreadyTested[row[k].y][row[k].x] = 1;
	    }
	  }
	}
      }
    }
  }

  for (int y = 0; y < SCREEN_HEIGHT; y++) {
    for (int x = 0; x < SCREEN_WIDTH; x++) {
      alreadyTested[y][x] = 0;
    }
  }

  for (uint i = 0; i < surfaces.size(); i++) {
    vector<Pixel> projectedVertices(surfaces[i].vertices.size());
    //cout << "scalar: " << surfaces[i].scalar << endl;
    if (surfaces[i].scalar) {
      for (uint v = 0; v < surfaces[i].vertices.size(); v++) {
	VertexShader(surfaces[i].vertices[v], projectedVertices[v]);
      }
      vector<Pixel> leftPixels;
      vector<Pixel> rightPixels;
      ComputePolygonRows(projectedVertices, leftPixels, rightPixels);
      for (uint j = 0; j < leftPixels.size(); j++) {
	vector<Pixel> row(rightPixels[j].x-leftPixels[j].x+1);
	Interpolate(leftPixels[j], rightPixels[j], row);
	for (uint k = 0; k < row.size(); k++) {
	  if ((row[k].x >= 0) && (row[k].x < SCREEN_WIDTH) && (row[k].y >= 0) && (row[k].y < SCREEN_HEIGHT)) {
	    if ((row[k].zinv <= depthBuffer[row[k].y][row[k].x]) && (stencilBuffer[row[k].y][row[k].x] > 0) && (!alreadyTested[row[k].y][row[k].x])) {
	      stencilBuffer[row[k].y][row[k].x]--;
	      //PutPixelSDL(screen, row[k].x, row[k].y, vec3(0, 0, 1));
	      alreadyTested[row[k].y][row[k].x] = 1;
	    }
	  }
	}
      }
    }
  }
}

void DrawShadows(screen* screen) {
  for (int y = 0; y < SCREEN_HEIGHT; y++) {
    for (int x = 0; x < SCREEN_WIDTH; x++) {
      if (stencilBuffer[y][x] == 0) {
	PutPixelSDL(screen, x, y, colourBuffer[y][x]);
      }
    }
  }
}

void Clip(const vector<Vertex>& vertices, vector<Vertex>& in) {
  vector<Vertex> inZero(0);
  vector<Vertex> inX(0);
  vector<Vertex> inY(0);
  ClipAtWEqualsZero(vertices, inZero);
  if (inZero.size() != 0) {
    ClipAtWEqualsPlane(inZero, 1, inX);
    if (inX.size() != 0) {
      ClipAtWEqualsPlane(inX, 2, inY);
      if (inY.size() != 0) {
	ClipAtWEqualsPlane(inY, 3, in);
	/*for (uint i = 0; i < in.size(); i++) {
	  cout << "in: (" << in[i].position.x << "," << in[i].position.y << "," << in[i].position.z << "," << in[i].position.w << ")" << endl;
	}*/
      }
    }
  }
}

void ClipAtWEqualsZero(const vector<Vertex>& vertices, vector<Vertex>& in) {
  int n = vertices.size();
  Vertex currentVertex;
  Vertex prevVertex = vertices[n-1];
  int currentDot;
  int prevDot;
  float t;
  vec4 I;
  float wPlane = 0.00001;
  Vertex inVertex = vertices[0];

  for (int i = 0; i < n; i++) {
    currentVertex = vertices[i];
    //cout << "currentVertex: (" << currentVertex.position.x << "," << currentVertex.position.y << "," << currentVertex.position.z << "," << currentVertex.position.w << ")" << endl;
    //cout << "prevVertex: (" << prevVertex.position.x << "," << prevVertex.position.y << "," << prevVertex.position.z << "," << prevVertex.position.w << ")" << endl;

    if (currentVertex.position.w < wPlane) {
      currentDot = -1;
    } else {
      currentDot = 1;
    }

    if (prevVertex.position.w < wPlane) {
      prevDot = -1;
    } else {
      prevDot = 1;
    }

    if (prevDot*currentDot < 0) {
      t = abs((wPlane-prevVertex.position.w)/(prevVertex.position.w-currentVertex.position.w));
      I = prevVertex.position+(t*(currentVertex.position-prevVertex.position));
      inVertex.position = I;
      in.push_back(inVertex);
      //cout << "intersection" << endl;
      //cout << "in: (" << in.back().position.x << "," << in.back().position.y << "," << in.back().position.z << "," << in.back().position.w << ")" << endl;
    }

    if (currentDot > 0) {
      in.push_back(currentVertex);
      //cout << "vertex in" << endl;
      //cout << "in: (" << in.back().position.x << "," << in.back().position.y << "," << in.back().position.z << "," << in.back().position.w << ")" << endl;
    }

    prevVertex = currentVertex;
  }
}

void ClipAtWEqualsPlane(const vector<Vertex>& vertices, int plane, vector<Vertex>& in) {
  int n = vertices.size();
  Vertex currentVertex;
  Vertex prevVertex = vertices[n-1];
  float currentAxis;
  float prevAxis;
  int currentDot;
  int prevDot;
  float t;
  vec4 I;
  Vertex inVertex = vertices[0];
  vector<Vertex> inPos(0);

  for (int i = 0; i < n; i++) {
    currentVertex = vertices[i];
    //cout << "currentVertex: (" << currentVertex.position.x << "," << currentVertex.position.y << "," << currentVertex.position.z << "," << currentVertex.position.w << ")" << endl;

    if (plane == 1) {
      currentAxis = currentVertex.position.x;
      prevAxis = prevVertex.position.x;
    } else if (plane == 2) {
      currentAxis = currentVertex.position.y;
      prevAxis = prevVertex.position.y;
    } else if (plane == 3) {
      currentAxis = currentVertex.position.z;
      prevAxis = prevVertex.position.z;
    }

    if (currentAxis <= currentVertex.position.w) {
      currentDot = 1;
    } else {
      currentDot = -1;
    }

    if (prevAxis <= prevVertex.position.w) {
      prevDot = 1;
    } else {
      prevDot = -1;
    }

    if (prevDot*currentDot < 0) {
      t = (prevVertex.position.w-prevAxis)/((prevVertex.position.w-prevAxis)-(currentVertex.position.w-currentAxis));
      I = prevVertex.position+(t*(currentVertex.position-prevVertex.position));
      inVertex.position = I;
      inPos.push_back(inVertex);
      //cout << "in: (" << inPos.back().position.x << "," << inPos.back().position.y << "," << inPos.back().position.z << "," << inPos.back().position.w << ")" << endl;
    }

    if (currentDot > 0) {
      inPos.push_back(currentVertex);
      //cout << "in: (" << inPos.back().position.x << "," << inPos.back().position.y << "," << inPos.back().position.z << "," << inPos.back().position.w << ")" << endl;
    }

    prevVertex = currentVertex;
  }

  n = inPos.size();
  if (n != 0) {
    prevVertex = inPos[n-1];
    for (int i = 0; i < n; i++) {
      currentVertex = inPos[i];
      //cout << "currentVertex: (" << currentVertex.position.x << "," << currentVertex.position.y << "," << currentVertex.position.z << "," << currentVertex.position.w << ")" << endl;

      if (plane == 1) {
	currentAxis = currentVertex.position.x;
	prevAxis = prevVertex.position.x;
      } else if (plane == 2) {
	currentAxis = currentVertex.position.y;
	prevAxis = prevVertex.position.y;
      } else if (plane == 3) {
	currentAxis = currentVertex.position.z;
	prevAxis = prevVertex.position.z;
      }

      if (-currentAxis <= currentVertex.position.w) {
	currentDot = 1;
      } else {
	currentDot = -1;
      }

      if (-prevAxis <= prevVertex.position.w) {
	prevDot = 1;
      } else {
	prevDot = -1;
      }

      if (prevDot*currentDot < 0) {
	t = (prevVertex.position.w+prevAxis)/((prevVertex.position.w+prevAxis)-(currentVertex.position.w+currentAxis));
	I = prevVertex.position+(t*(currentVertex.position-prevVertex.position));
	inVertex.position = I;
	in.push_back(inVertex);
	//cout << "in: (" << in.back().position.x << "," << in.back().position.y << "," << in.back().position.z << "," << in.back().position.w << ")" << endl;
      }

      if (currentDot > 0) {
	in.push_back(currentVertex);
	//cout << "in: (" << in.back().position.x << "," << in.back().position.y << "," << in.back().position.z << "," << in.back().position.w << ")" << endl;
      }

      prevVertex = currentVertex;
    }
  }
}