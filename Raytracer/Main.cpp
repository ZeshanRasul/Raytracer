#include "FreeImage.h"
#include <string>
#include "Ray.h"
#include "Camera.h"
#include "Intersection.h"

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

Ray ShootRay(Camera cam, int i, int j, int width, int height)
{
	// Create coordinate system from camera eye - center and up
	vec3 a = cam.eyePos - cam.center;
	vec3 b = cam.up;

	vec3 w = normalize(a);
	vec3 u = normalize(cross(b, w));
	vec3 v = cross(w, u);

	// Create new ray
	float alpha = tan(1 / 2) * ((i - (width / 2)) / (width / 2));
	float beta = tan(cam.fovY / 2) * (((height / 2) - j) / (height / 2));

	vec3 direction = normalize((alpha * u) + (beta * v) - w);
	vec3 origin = cam.eyePos;
	Ray ray(vec4(origin, 1), vec4(direction, 0));

	// Return ray
	return ray;

}

Intersection FindIntersection(/*Scene scene,*/ Ray ray)
{
	// For each sphere in scene and for each triangle in scene
	// i < scene->spheres.length() test for intersection
	// call sphereIntersection and triangleIntersection respectively
}