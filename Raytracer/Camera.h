#pragma once
#include "glm/glm.hpp"

using namespace glm;

class Camera
{
public:
	Camera(vec3 eyePos, vec3 center, vec3 up)
		:
		eyePos(eyePos),
		center(center),
		up(up)
	{
	}

	~Camera()
	{

	}

	vec3 eyePos, center, up;
};