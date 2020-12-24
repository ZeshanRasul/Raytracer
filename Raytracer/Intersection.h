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
				 float hitObjectShininess)
		:
		didHit(didHit),
		intersectionPoint(intersectionPoint),
		hitObjectDiffuse(hitObjectDiffuse),
		hitObjectSpecular(hitObjectSpecular),
		hitObjectEmission(hitObjectEmission),
		hitObjectShininess(hitObjectShininess)
	{		
	}

	~Intersection()
	{
	}

	bool didHit;
	vec4 intersectionPoint;
	vec3 hitObjectDiffuse;
	vec3 hitObjectSpecular;
	vec3 hitObjectEmission;
	float hitObjectShininess;
};