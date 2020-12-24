#pragma once
#include <vector>
#include "glm/glm.hpp"
#include "Sphere.h"
#include "Triangle.h"

class Scene
{
public:
	std::vector<Sphere> spheres;
	std::vector<Triangle> triangles;
	// std::vector<Light> lights;
};