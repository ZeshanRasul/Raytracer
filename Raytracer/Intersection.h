#pragma once
#include "glm/glm.hpp"

using namespace glm;

class Intersection
{
public:
	Intersection(bool didHit,
		vec4 intersectionPoint,
		vec3 hitObjectDiffuse,
		vec3 hitObjectSpecular,
		vec3 hitObjectEmission,
		float hitObjectShininess,
		vec3 hitObjectAmbient,
		vec4 hitObjectNormal,
		bool hitObjectIsSphere)
		:
		didHit(didHit),
		intersectionPoint(intersectionPoint),
		hitObjectDiffuse(hitObjectDiffuse),
		hitObjectSpecular(hitObjectSpecular),
		hitObjectEmission(hitObjectEmission),
		hitObjectShininess(hitObjectShininess),
		hitObjectAmbient(hitObjectAmbient),
		hitObjectNormal(hitObjectNormal),
		hitObjectIsSphere(hitObjectIsSphere)
	{		
	}

	~Intersection()
	{
	}

	bool didHit = false;
	vec4 intersectionPoint;
	vec3 hitObjectDiffuse;
	vec3 hitObjectSpecular;
	vec3 hitObjectEmission;
	vec3 hitObjectAmbient;
	float hitObjectShininess;
	vec4 hitObjectNormal;
	bool hitObjectIsSphere;
};