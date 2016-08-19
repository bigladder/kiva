/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef PIXELCOUNTER_H_
#define PIXELCOUNTER_H_

#define GLEW_STATIC
#include <GL/glew.h>

#include "Geometry.hpp"

#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>
#include <time.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

namespace Kiva {

class PixelCounter
{
public:

  PixelCounter(int size = 512, int numQueries = 8, bool showWindow = false);
  ~PixelCounter();


  double getAreaRatio(double orientation, double azimuth, double altitude, std::vector<Polygon3> shading, std::vector<Polygon3> shaded, int index);
  int retrievePixelCount(int index);
  bool queryComplete(int index);

private:

  glm::vec3 sun;
  glm::mat4 viewMatrix;
  glm::mat4 projectionMatrix;

  int size;
    int numQueries;
  bool showWindow;
    GLFWwindow* window;
    GLuint vao[2];
    GLuint vbo[2];
    GLuint ebo[2];
    std::vector<GLuint> queries, fbo, rbo1, rbo2;


  void drawRectangles(std::vector<Polygon3> polygons, float color, GLuint &vao, GLuint &vbo, GLuint &ebo);


};

}

#endif /* PIXELCOUNTER_H_ */
