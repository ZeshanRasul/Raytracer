#pragma once
#include "glm/glm.hpp"
#include "Object.h"

using namespace glm;

class Triangle : public Object
{
public:
	Triangle(vec3 vertex0, vec3 vertex1, vec3 vertex2, vec3 diffuse, vec3 specular, vec3 emission, float shininess, vec3 ambient)
		:
		vertex0(vertex0),
		vertex1(vertex1),
		vertex2(vertex2),
		diffuse(diffuse),
		specular(specular),
		emission(emission),
		shininess(shininess),
		ambient(ambient)
	{}

	vec3 vertex0, vertex1, vertex2;
	vec3 normal;
	vec3 diffuse, specular, emission, ambient;
	float shininess;

};