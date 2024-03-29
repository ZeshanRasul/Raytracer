#pragma once
#include "glm/glm.hpp"

using namespace glm;

class DirectionalLight
{
public:
	DirectionalLight(vec3 direction, vec3 colour)
		:
		direction(direction),
		colour(colour)
	{}

	vec3 direction;
	vec3 colour;
};