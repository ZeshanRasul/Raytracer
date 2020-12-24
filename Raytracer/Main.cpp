#include "FreeImage.h"
#include <string>
#include "Ray.h"
#include "Camera.h"
#include "Intersection.h"
#include "Scene.h"
#include "Sphere.h"
#include "Triangle.h"

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

}

Intersection FindIntersection(Scene scene, Ray ray)
{
	// For each sphere in scene and for each triangle in scene
	// i < scene->spheres.length() test for intersection
	// call sphereIntersection and triangleIntersection respectively
	// return closest intersection
	
	float minDist = INFINITY;
	Object hitObject;
	bool didHit = false;
	float t;
	vec3 hitObjectDiffuse, hitObjectSpecular, hitObjectEmission;
	float hitObjectShininess;


	for (int i = 0; i < scene.spheres.size(); i++)
	{
		t = CheckSphereIntersection(scene.spheres[i], ray);

		if (t < minDist && t > 0)
		{
			hitObject = scene.spheres[i];
			hitObjectDiffuse = scene.spheres[i].diffuse;
			hitObjectSpecular = scene.spheres[i].specular;
			hitObjectEmission = scene.spheres[i].emission;
			hitObjectShininess = scene.spheres[i].shininess;
		}
	}

	// Same loop as above for triangles
	
	if (minDist < INFINITY)
	{
		didHit = true;
	}

	// Calculate intersection point
	vec4 intersectionPoint = ray.origin + ray.direction * t;

	return Intersection(didHit, intersectionPoint, hitObjectDiffuse, hitObjectSpecular, hitObjectEmission, hitObjectShininess);
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

	Sphere sphere0(vec4(0, 0, 0, 1), 0.666, vec3(0.67, 0.33, 0.93), vec3(0.2, 0.2, 0.2), vec3(0.1, 0.1, 0.1), 20.0f);
	//Create Triangle here

	scene.spheres.push_back(sphere0);

	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			// Shoot Ray
			ShootRay(camera, i, j, width, height);
		}
	}


	FreeImage_Initialise();
	FIBITMAP* img = FreeImage_ConvertFromRawBits(pixels, width, height, 3 * width, 24, 0xFF0000, 0x00FF00, 0x0000FF, false);
	FreeImage_Save(FIF_PNG, img, outputFilename.c_str(), 0);
}
