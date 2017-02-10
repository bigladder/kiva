/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#include "PixelCounter.hpp"

namespace Kiva {

// Shader sources
const GLchar* vertexSource =
    "#version 120\n"
    "attribute vec3 position;"
  "attribute vec3 color;"
  "varying vec3 Color;"
    "uniform mat4 projectionMatrix;"
  "uniform mat4 viewMatrix;"
  "uniform mat4 model;"
    "void main() {"
  "   Color = color;"
    "   gl_Position = projectionMatrix*viewMatrix*vec4(position, 1.0);"
    "}";
const GLchar* fragmentSource =
    "#version 120\n"
  "varying vec3 Color;"
    "void main() {"
    "   gl_FragColor = vec4(Color, 1.0);"
    "}";

static const double PI = 4.0*atan(1.0);

PixelCounter::PixelCounter(int size, int numQueries, bool showWindow) :
    size(size), numQueries(numQueries), showWindow(showWindow)
{

    if (!glfwInit())
        throw;

    if (!showWindow)
    {
      glfwWindowHint(GLFW_VISIBLE, GL_FALSE);
      window = glfwCreateWindow(1, 1, "OpenGL", NULL, NULL);
    }
    else
      window = glfwCreateWindow(size, size, "OpenGL", NULL, NULL);

    // Make the window's context current
    glfwMakeContextCurrent(window);

    // Setup GLEW
    glewExperimental = GL_TRUE;
    glewInit();

    if (!glewIsSupported("GL_VERSION_2_1"))
      throw;
    if (!glewIsSupported("GL_EXT_framebuffer_object"))
      throw;

    if (!window)
    {
        glfwTerminate();
        throw;
    }

    // Set number of asynchronous queries to the number of available cores
    queries.resize(numQueries);
    fbo.resize(numQueries);
    rbo1.resize(numQueries);
    rbo2.resize(numQueries);

    // Setup framebuffer object for off-screen rendering
    if (!showWindow)
    {
        glGenFramebuffersEXT(fbo.size(), &fbo[0]);
        glGenRenderbuffersEXT(rbo1.size(), &rbo1[0]);
        glGenRenderbuffersEXT(rbo2.size(), &rbo2[0]);

        for (std::size_t q = 0; q < numQueries; q++)
        {
            glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo[q]);

            glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rbo1[q]);
            glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGBA8, size, size);
            glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_RENDERBUFFER_EXT, rbo1[q]);

            glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rbo2[q]);
            glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24, size, size);
            glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, rbo2[q]);

            GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
            if (status != GL_FRAMEBUFFER_COMPLETE_EXT)
              throw;

        }
    }

    glEnable(GL_DEPTH_TEST);
    glGenQueries(queries.size(), &queries[0]);

    if (!showWindow)
      glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

}

PixelCounter::~PixelCounter()
{
    glDeleteVertexArrays(1, &vao[0]);
    glDeleteVertexArrays(1, &vao[1]);
    glDeleteBuffers(1, &vbo[0]);
    glDeleteBuffers(1, &vbo[1]);
    glDeleteBuffers(1, &ebo[0]);
    glDeleteBuffers(1, &ebo[1]);
    glDeleteFramebuffers(numQueries, &fbo[0]);
    glDeleteRenderbuffers(numQueries, &rbo1[0]);
    glDeleteRenderbuffers(numQueries, &rbo2[0]);
    glfwTerminate();
}

double PixelCounter::getAreaRatio(double orientation, double azimuth, double altitude, std::vector<Polygon3> shading, std::vector<Polygon3> shaded, int index)
{

  float hypotenuse = cos(altitude);
  float sunX = hypotenuse*sin(azimuth-orientation);
  float sunY = hypotenuse*cos(azimuth-orientation);
  float sunZ = sin(altitude);

  glm::vec3 sunPosition(sunX, sunY, sunZ);
  sun = sunPosition;

    viewMatrix = glm::lookAt(
         sun,
           glm::vec3(0.0f, 0.0f, 0.0f),
           glm::vec3(0.0f, 0.0f, 1.0f)
       );

  // This function draws the geometry and also begins the query for the pixel count

  float left = FLT_MAX;
  float right = -left;
  float bottom = left;
  float top = right;
  float near = left;
  float far = right;

  double totalArea = 0;

  // Find the extents of the surface in eye space
  for (std::size_t p = 0; p < shaded.size(); p++)
  {
    for (std::size_t v = 0; v < shaded[p].outer().size(); v++)
    {
      float x = shaded[p].outer()[v].get<0>();
      float y = shaded[p].outer()[v].get<1>();
      float z = shaded[p].outer()[v].get<2>();

      glm::vec4 point(x,y,z,0);
      glm::vec4 trans = viewMatrix * point;
      trans.z = -trans.z;
      if (trans.x < left)   left = trans.x;
      if (trans.x > right)  right = trans.x;
      if (trans.y < bottom) bottom = trans.y;
      if (trans.y > top)    top = trans.y;
      if (trans.z < near)
        near = trans.z;
      if (trans.z > far)
        far = trans.z;

    }

    totalArea += boost::geometry::distance(shaded[p].outer()[0],shaded[p].outer()[1])
          *boost::geometry::distance(shaded[p].outer()[1],shaded[p].outer()[2]);
  }

  // Near & far extents can be set by the shading surfaces, as well.
  for (std::size_t p = 0; p < shading.size(); p++)
  {
    for (std::size_t v = 0; v < shading[p].outer().size(); v++)
    {
      float x = shading[p].outer()[v].get<0>();
      float y = shading[p].outer()[v].get<1>();
      float z = shading[p].outer()[v].get<2>();

      glm::vec4 point(x,y,z,0);
      glm::vec4 trans = viewMatrix * point;
      //trans.z = -trans.z;
      if (trans.z < near)
        near = trans.z;
      if (trans.z > far)
        far = trans.z;

    }
  }

  // Grow horizontal extents of view frustrum by one pixel on each side
  float deltaX = (right - left)/size;
  left -= deltaX;
  right += deltaX;

  // Grow vertical extents of view frustrum by one pixel on each side
  float deltaY = (top - bottom)/size;
  bottom -= deltaY;
  top += deltaY;

  // TODO: Detect depth
  float depth = far - near;
  // Grow depth extents of view frustrum by 4 times building size
  float deltaZ = depth*10;
  near -= deltaZ;
  far += deltaZ;

  double areaPerPixel = (right - left)*(top - bottom)/(size*size);

  projectionMatrix = glm::ortho(left, right, bottom, top, near, far);

  // Render
    if (showWindow)
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    else
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo[index]);



    if (showWindow)
  {
    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window))
    {
      // Render here //
      // Clear the screen to black
        glViewport( 0, 0, size, size);

      glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      // Draw
      glDepthFunc(GL_LESS);
      drawRectangles(shading, 0.5, vao[0], vbo[0], ebo[0]);

      glDepthFunc(GL_LEQUAL);
      drawRectangles(shaded, 0.0, vao[1], vbo[1], ebo[1]);

      // Swap front and back buffers //
      glfwSwapBuffers(window);

      // Poll for and process events //
      glfwPollEvents();

      if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
    }
    glfwSetWindowShouldClose(window, GL_FALSE);
  }

    glViewport( 0, 0, size, size);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDepthFunc(GL_LESS);
    drawRectangles(shading, 0.5, vao[0], vbo[0], ebo[0]);

    glBeginQuery(GL_SAMPLES_PASSED, queries[index]);

    glDepthFunc(GL_LEQUAL);
    drawRectangles(shaded, 0.0, vao[1], vbo[1], ebo[1]);

    glEndQuery(GL_SAMPLES_PASSED);

    return areaPerPixel/totalArea;
}

int PixelCounter::retrievePixelCount(int index)
{
  GLint result;

  glGetQueryObjectiv(queries[index], GL_QUERY_RESULT, &result);

  return result;

}

bool PixelCounter::queryComplete(int index)
{
  GLint result;

  glGetQueryObjectiv(queries[index], GL_QUERY_RESULT_AVAILABLE, &result);

  return result;

}


void PixelCounter::drawRectangles(std::vector<Polygon3> polygons, float color, GLuint &vao, GLuint &vbo, GLuint &ebo)
{
    float vertices[polygons.size()*4*6];

    int vNum = 0;
  for (std::size_t p = 0; p < polygons.size(); p++)
    {
      for (std::size_t v = 0; v < polygons[p].outer().size(); v++)
      {
        vertices[vNum] = polygons[p].outer()[v].get<0>();
        vNum++;
        vertices[vNum] = polygons[p].outer()[v].get<1>();
        vNum++;
        vertices[vNum] = polygons[p].outer()[v].get<2>();
        vNum++;
        vertices[vNum] = color;
        vNum++;
        vertices[vNum] = color;
        vNum++;
        vertices[vNum] = color;
        vNum++;
      }
    }

  GLsizei bufferSize = sizeof(vertices)/sizeof(vertices[0]);

  // Create Vertex Array Object
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Create a Vertex Buffer Object and copy the vertex data to it
    glGenBuffers(1, &vbo);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, bufferSize*sizeof(float), vertices, GL_STATIC_DRAW);

    // Create an element array
    glGenBuffers(1, &ebo);

    GLuint numVertices = bufferSize/6;

    GLuint numElements = numVertices/4*6;

    GLuint elements[numElements];

    for (GLuint r = 0; r < numVertices/4; r++)
    {
      elements[r*6+0] = r*4+0;
      elements[r*6+1] = r*4+1;
      elements[r*6+2] = r*4+2;
      elements[r*6+3] = r*4+2;
      elements[r*6+4] = r*4+3;
      elements[r*6+5] = r*4+0;
    }

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);

    GLuint shaderProgram = glCreateProgram();

    // Create and compile the vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexSource, NULL);
    glCompileShader(vertexShader);
    glAttachShader(shaderProgram, vertexShader);

    // Create and compile the fragment shader
    if (showWindow)
    {
        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
        glCompileShader(fragmentShader);
        glAttachShader(shaderProgram, fragmentShader);
    }

    // Link the vertex and fragment shader into a shader program
    glLinkProgram(shaderProgram);
    glUseProgram(shaderProgram);

    // Set up view matrix and projection matrix
    GLint uniView = glGetUniformLocation(shaderProgram, "viewMatrix");
    glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(viewMatrix));

    GLint uniProj = glGetUniformLocation(shaderProgram, "projectionMatrix");
    glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(projectionMatrix));

    // Specify the layout of the vertex data
    GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
    glEnableVertexAttribArray(posAttrib);
    glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);

    GLint colAttrib = glGetAttribLocation(shaderProgram, "color");
    glEnableVertexAttribArray(colAttrib);
    glVertexAttribPointer(colAttrib, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)(3*sizeof(float)));

    //GLint status;
    //glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &status);

    //char buffer[512];
    //glGetShaderInfoLog(vertexShader, 512, NULL, buffer);

    //std::cout << buffer << std::endl;

    glDrawElements(GL_TRIANGLES, numElements, GL_UNSIGNED_INT, 0);

}
}
