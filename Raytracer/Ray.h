#pragma once
#include "glm/glm.hpp"

using namespace glm;

class Ray
{
public:
	Ray(vec4 origin, vec3 direction)
		:
		origin(origin),
		direction(direction)
	{
	}

	vec4 origin;
	vec3 direction;
private:
};