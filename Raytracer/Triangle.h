#pragma once
#include "glm/glm.hpp"
#include "Object.h"

using namespace glm;

class Triangle : public Object
{
public:
	Triangle(vec3 vertex0, vec3 vertex1, vec3 vertex2, vec3 diffuse, vec3 specular, vec3 emission, float shininess, vec3 ambient)
		:
		vertex0(vec4(vertex0, 1)),
		vertex1(vec4(vertex1, 1)),
		vertex2(vec4(vertex2, 1)),
		diffuse(diffuse),
		specular(specular),
		emission(emission),
		shininess(shininess),
		ambient(ambient)

	{
		normalA = normalize(cross((vertex2 - vertex0), (vertex1 - vertex0)));
		normalB = normalize(cross((vertex0 - vertex1), (vertex2 - vertex1)));
		normalC = normalize(cross((vertex1 - vertex2), (vertex0 - vertex2)));
	}

	Triangle()
	{}

	void SetNormal(vec3 newNormal)
	{
		normal = newNormal;
	}


	vec4 vertex0, vertex1, vertex2;
	vec3 normalA;
	vec3 normalB;
	vec3 normalC;
	vec3 normal;
	vec3 diffuse, specular, emission, ambient;
	float shininess;

};