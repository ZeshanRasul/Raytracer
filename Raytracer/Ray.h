#pragma once
#include "glm/glm.hpp"

using namespace glm;

class Ray
{
public:
	Ray(vec3 origin, vec3 direction)
		:
		origin(origin),
		direction(direction)
	{}

	~Ray()
	{}

	vec3 origin;
	vec3 direction;
};