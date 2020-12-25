#pragma once
#include "glm/glm.hpp"

using namespace glm;

class PointLight
{
public:
	PointLight(vec3 position, vec3 colour)
		:
		position(position),
		colour(colour)
	{}

	vec3 position;
	vec3 colour;
};