#pragma once
#include "glm/glm.hpp"
#include "Object.h"

using namespace glm;

class Sphere : public Object
{
public:
	Sphere(vec3 center, float radius, vec3 diffuse, vec3 specular, vec3 emission, float shininess, vec3 ambient)
		:
		center(center),
		radius(radius),
		diffuse(diffuse),
		specular(specular),
		emission(emission),
		shininess(shininess),
		ambient(ambient)
	{}
	
	Sphere()
	{}

	~Sphere()
	{};

	void SetNormal(vec3 newNormal)
	{
		normal = newNormal;
	}

	vec3 center;
	float radius;
	vec3 normal;
	vec3 diffuse, specular, emission, ambient;
	float shininess;
};