#pragma once
#include "glm/glm.hpp"

using namespace glm;

class Light
{
public:
	Light(vec3 position, vec3 colour)
		:
		position(position),
		colour(colour)
	{}

	vec3 position;
	vec3 colour;
};