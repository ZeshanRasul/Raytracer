#pragma once
#include "glm/glm.hpp"

using namespace glm;

class Triangle
{
public:
	Triangle(vec3 vertex0, vec3 vertex1, vec3 vertex2, vec3 diffuse, vec3 specular, vec3 emission, float shininess)
		:
		vertex0(vertex0),
		vertex1(vertex1),
		vertex2(vertex2),
		diffuse(diffuse),
		specular(specular),
		emission(emission),
		shininess(shininess)
	{}

	vec3 vertex0, vertex1, vertex2;
	vec3 diffuse, specular, emission;
	float shininess;

};