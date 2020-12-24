#include "FreeImage.h"
#include <string>
#include "Ray.h"
#include "Camera.h"
#include "Intersection.h"
#include "Scene.h"
#include "Sphere.h"
#include "Triangle.h"
#include <vector>

std::vector <mat4> modelViewStack;

void pushMatrix(mat4 mat)
{
	modelViewStack.push_back(mat);
}

void popMatrix(mat4& mat)
{
	if (modelViewStack.size())
	{
		mat = modelViewStack.back();
		modelViewStack.pop_back();
	}
	else
	{
		mat = mat4(1.0);
	}
}

mat4 lookAt(const vec3& eye, const vec3& center, const vec3& up)
{
	vec3 a = eye - center;
	vec3 w = normalize(a);

	vec3 crossUpW = glm::cross(up, w);
	vec3 u = glm::normalize(crossUpW);
	//u = -u;
	vec3 v = cross(w, u);

	mat4 lookAtMatrix = mat4(u.x, u.y, u.z, -glm::dot(u, eye),
							 v.x, v.y, v.z, -glm::dot(v, eye),
							 w.x, w.y, w.z, -glm::dot(w, eye),
							 0, 0, 0, 1);

	return glm::transpose(lookAtMatrix);

}


Ray ShootRay(Camera cam, int i, int j, int width, int height)
{
	// Create coordinate system from camera eye - center and up
	vec3 a = cam.eyePos - cam.center;
	vec3 b = cam.up;

	vec3 w = normalize(a);
	vec3 u = normalize(cross(b, w));
	vec3 v = cross(w, u);

	// Create new ray
	float alpha = float(tan(1 / 2) * ((i - (width / 2)) / (width / 2)));
	float beta = tan(cam.fovY / 2) * (((height / 2) - j) / (height / 2));

	vec3 direction = normalize((alpha * u) + (beta * v) - w);
	vec3 origin = cam.eyePos;
	Ray ray(vec4(origin, 1), vec4(direction, 0));

	// Return ray
	return ray;

}

float CheckSphereIntersection(Sphere sphere, Ray ray)
{
	// Check discriminant
	// If 2 positive roots take smaller
	// If both roots the same tangent 
	// If 1 positive 1 negative choose positive
	// If complex roots no intersection
	float a = dot(ray.direction, ray.direction);
	float b = 2 * (dot(ray.direction, (ray.origin - sphere.center)));
	float c = dot((ray.origin - sphere.center), (ray.origin - sphere.center)) - (sphere.radius * sphere.radius);

	float returnVal;
	int roots = 0;

	float discriminant = (b * b) - (4 * a * c);
	if (discriminant < 0)
	{
		return 0.0f;
	} 
	else if (discriminant > 0)
	{
		roots = 2;
	}

	if (roots == 2)
	{
		float t = (-b + sqrt(((b * b) - (4 * a * c))) / (2 * a));
		float t1 = (-b - sqrt(((b * b) - (4 * a * c))) / (2 * a));

		if (t > 0 && t1 > 0)
		{
			t < t1 ?  returnVal = t : returnVal = t1;
		}

		if (t > 0 && t1 < 0)
		{
			returnVal = t;
		}

		if (t < 0 && t1 > 0)
		{
			returnVal = t1;
		}

		vec4 intersecPoint = ray.origin + (ray.direction * returnVal);
		sphere.normal = normalize(intersecPoint - sphere.center);

		return returnVal;
	} 
	else if (roots == 1)
	{
		float t = (-b + sqrt(((b * b) - (4 * a * c))) / (2 * a));
		returnVal = t;

		vec4 intersecPoint = ray.origin + (ray.direction * returnVal);
		sphere.normal = normalize(intersecPoint - sphere.center);

		return returnVal;
	}

}

Intersection FindIntersection(Scene scene, Ray ray)
{
	// For each sphere in scene and for each triangle in scene
	// i < scene->spheres.length() test for intersection
	// call sphereIntersection and triangleIntersection respectively
	// return closest intersection
	
	float minDist = INFINITY;
	//Object hitObject;
	bool didHit = false;
	float t = minDist;
	vec3 hitObjectDiffuse = vec3(0, 0, 0);
	vec3 hitObjectSpecular = vec3(0, 0, 0);
	vec3 hitObjectEmission = vec3(0, 0, 0);
	float hitObjectShininess = 0.0f;


	for (int i = 0; i < scene.spheres.size(); i++)
	{
		t = CheckSphereIntersection(scene.spheres[i], ray);

		if (t < minDist && t > 0)
		{
		//	hitObject = scene.spheres[i];
			hitObjectDiffuse = scene.spheres[i].diffuse;
			hitObjectSpecular = scene.spheres[i].specular;
			hitObjectEmission = scene.spheres[i].emission;
			hitObjectShininess = scene.spheres[i].shininess;
			
			minDist = t;
			didHit = true;

		}
	}

	// Same loop as above for triangles
	

	// Calculate intersection point
	vec4 intersectionPoint = ray.origin + ray.direction * t;

	return Intersection(didHit, intersectionPoint, hitObjectDiffuse, hitObjectSpecular, hitObjectEmission, hitObjectShininess);
}

vec3 FindColour(Intersection intersection)
{
	return intersection.hitObjectDiffuse;
}

int main()
{
	const int width = 128;
	const int height = 128;
	unsigned char pixels[width * height * 3] = { 0 };
	std::string outputFilename = "Raytracer.png";

	vec3 eyePosition = vec3(0, 0, 4);
	vec3 center = vec3(0, 0, 0);
	vec3 up = vec3(0, 1, 0);
	float fovY = 30;
	// Create new Camera with default values 
	Camera camera(eyePosition, center, up, fovY);

	// Create new Scene and add Sphere and then Triangle
	Scene scene;

	Sphere sphere0(vec4(0, 0, 0, 1), 0.1f, vec3(0.67, 0.33, 0.93), vec3(0.2, 0.2, 0.2), vec3(0.1, 0.1, 0.1), 20.0f);
	//Create Triangle here
	scene.spheres.push_back(sphere0);

	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			// Shoot Ray
			Ray ray = ShootRay(camera, i, j, width, height);
			Intersection intersection = FindIntersection(scene, ray);
			vec3 colour = FindColour(intersection);
			unsigned char col[] = { 0, 0, 0 };
			col[0] = colour[0] * 255;
			col[1] = colour[1] * 255;
			col[2] = colour[2] * 255;
			memcpy(&pixels[((j * width) + i)  * 3], &col, 3);
		}
	}

	FreeImage_Initialise();
	FIBITMAP* img = FreeImage_ConvertFromRawBits((BYTE*)pixels, width, height, 3 * width, 24, 0xFF0000, 0x00FF00, 0x0000FF, false);
	FreeImage_Save(FIF_PNG, img, outputFilename.c_str(), 0);
	FreeImage_DeInitialise();
}
