#pragma once
#include "glm/glm.hpp"

using namespace glm;

class Sphere
{
public:
	Sphere(vec3 center, float radius, vec3 diffuse, vec3 specular, vec3 emission, float shininess)
		:
		center(center),
		radius(radius),
		diffuse(diffuse),
		specular(specular),
		emission(emission),
		shininess(shininess)
	{}
	
	~Sphere()
	{};

	vec3 center;
	float radius;
	vec3 diffuse, specular, emission;
	float shininess;
};