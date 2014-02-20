/* PixelCounter.h is part of Kiva (Written by Neal Kruis)
 * Copyright (C) 2012-2013 Big Ladder Software <info@bigladdersoftware.com>
 * All rights reserved.
 *
 * Kiva is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kiva is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kiva.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PIXELCOUNTER_H_
#define PIXELCOUNTER_H_

#define GLEW_STATIC
#include <GL/glew.h>

#include "Geometry.h"

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


#endif /* PIXELCOUNTER_H_ */
