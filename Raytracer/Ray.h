#pragma once
#include "glm/glm.hpp"

using namespace glm;

class Ray
{
public:
	Ray(vec4 origin, vec4 direction)
		:
		origin(origin),
		direction(direction)
	{}

	~Ray()
	{}

	vec4 origin;
	vec4 direction;
};