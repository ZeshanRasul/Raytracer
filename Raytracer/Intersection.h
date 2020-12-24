#pragma once
#include "glm/glm.hpp"

using namespace glm;

class Intersection
{
public:
	Intersection(bool didHit,
		vec3 intersectionPoint,
		vec3 hitObjectDiffuse,
		vec3 hitObjectSpecular,
		vec3 hitObjectEmission,
		float hitObjectShininess,
		vec3 hitObjectAmbient,
		vec3 hitObjectNormal)
		:
		didHit(didHit),
		intersectionPoint(intersectionPoint),
		hitObjectDiffuse(hitObjectDiffuse),
		hitObjectSpecular(hitObjectSpecular),
		hitObjectEmission(hitObjectEmission),
		hitObjectShininess(hitObjectShininess),
		hitObjectAmbient(hitObjectAmbient),
		hitObjectNormal(hitObjectNormal)
	{		
	}

	~Intersection()
	{
	}

	bool didHit = false;
	vec3 intersectionPoint;
	vec3 hitObjectDiffuse;
	vec3 hitObjectSpecular;
	vec3 hitObjectEmission;
	vec3 hitObjectAmbient;
	float hitObjectShininess;
	vec3 hitObjectNormal;
};