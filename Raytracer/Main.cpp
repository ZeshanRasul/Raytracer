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

std::string viewMode = "sphere";

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
	vec3 b = up;

	vec3 w = normalize(a);
	vec3 u = normalize(cross(b, w));
	vec3 v = cross(w, u);

	mat4 lookAtMatrix = mat4(u.x, u.y, u.z, -glm::dot(u, eye),
							 v.x, v.y, v.z, -glm::dot(v, eye),
							 w.x, w.y, w.z, -glm::dot(w, eye),
							 0, 0, 0, 1);

	return transpose(lookAtMatrix);

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
//	float fovX = 2 * (1 / tan(cam.fovY / -2));
	float fovX = 2 * atan((height * tan(cam.fovY/ 2)) / width);
	float alpha = (tan(cam.fovY / 2) * ((float(i) - (float(width) / 2)) / (float(width) / 2)));
	float beta = tan(fovX / 2) * (((float(height) / 2) - float(j)) / (float(height) / 2));

	vec3 direction = normalize((alpha * u) + (beta * v) - w);
	vec3 origin = cam.eyePos;
	Ray ray(origin, direction);

	// Return ray
	return ray;

}

/*
float CheckSphereIntersection(Sphere sphere, Ray ray)
{

	//vec4 newCenter = lookAt(camera.eyePos, camera.center, camera.up) * sphere.center;
	vec3 newCenter = sphere.center;

	//float t = -(dot(ray.direction, (ray.origin - sphere.center))) + sqrt(((dot(ray.direction, (ray.origin - sphere.center))) * (dot(ray.direction, (ray.origin - sphere.center)))) - ((normalize(ray.origin - sphere.center) * normalize(ray.origin - sphere.center)) - (sphere.radius * sphere.radius)));

//	ray.origin = lookAt(camera.eyePos, camera.center, camera.up) * ray.origin;
//	ray.direction = lookAt(camera.eyePos, camera.center, camera.up) * ray.direction;

	// Check discriminant
	// If 2 positive roots take smaller
	// If both roots the same tangent 
	// If 1 positive 1 negative choose positive
	// If complex roots no intersection
	
	float a = dot(ray.direction, ray.direction);
	float b = 2 * (dot(ray.direction, (ray.origin - newCenter)));
	float c = dot((ray.origin - newCenter), (ray.origin - newCenter)) - (sphere.radius * sphere.radius);

	float returnVal = INFINITY;
	int roots = 0;

//	float discriminant = (b * b) - (4 * a * c);
	float discriminant = (dot(ray.direction, (ray.origin - newCenter))) * (dot(ray.direction, (ray.origin - newCenter))) - ((dot(ray.origin - newCenter, ray.origin - newCenter)) - (sphere.radius * sphere.radius));
	if (discriminant < 0)
	{
		return INFINITY;
	} 
	else if (discriminant > 0)
	{
		roots = 2;
	}
	else if (discriminant == 0)
	{
		roots = 1;
	}

	if (roots == 2)
	{
		float t = (-b + sqrt(((b * b) - (4 * a * c))) / (2 * a));
		float t1 = (-b - sqrt(((b * b) - (4 * a * c))) / (2 * a));
	//	float t = (-(dot(ray.direction, (ray.origin - sphere.center)))) + sqrt(((dot(ray.direction, (ray.origin - sphere.center))) * (dot(ray.direction, (ray.origin - sphere.center)))) - (dot(ray.origin - sphere.center, ray.origin - sphere.center)) - (sphere.radius * sphere.radius));
	//	float t1 = (-(dot(ray.direction, (ray.origin - sphere.center)))) - sqrt(((dot(ray.direction, (ray.origin - sphere.center))) * (dot(ray.direction, (ray.origin - sphere.center)))) - (dot(ray.origin - sphere.center, ray.origin - sphere.center)) - (sphere.radius * sphere.radius));

		if (t > 0 && t1 > 0)
		{
			t < t1 ?  returnVal = t : returnVal = t1;
		}

		else if (t > 0 && t1 < 0)
		{
			returnVal = t;
		}

		else if (t < 0 && t1 > 0)
		{
			returnVal = t1;
		}

		vec3 intersecPoint = ray.origin + (ray.direction * returnVal);
		sphere.SetNormal(normalize(intersecPoint - sphere.center));

		return returnVal;
	} 
	else if (roots == 1)
	{
		float t = -(dot(ray.direction, (ray.origin - sphere.center))) + sqrt(((dot(ray.direction, (ray.origin - sphere.center))) * (dot(ray.direction, (ray.origin - sphere.center)))) - (dot(ray.origin - sphere.center, ray.origin - sphere.center)) - (sphere.radius * sphere.radius));
		returnVal = t;

		vec3 intersecPoint = ray.origin + (ray.direction * returnVal);
		sphere.SetNormal(normalize(intersecPoint - sphere.center));

		return returnVal;
	}
	
	return INFINITY;


}
*/


float CheckSphereIntersection(Sphere sphere, Ray ray)
{
	vec3 oc = ray.origin - sphere.center;
	float a = dot(ray.direction, ray.direction);
	float b = 2 * dot(oc, ray.direction);
	float c = dot(oc, oc) - (sphere.radius * sphere.radius);
	float discriminant = (b * b) - (4 * a * c);
	if (discriminant < 0)
	{
		return -1.00f;
	}
	else
	{
		return (-b - (sqrt(discriminant)) / (2.00f * a));
	}

	// TODO Complete this
}

vec3 CheckTriangleIntersection(Triangle triangle, Ray ray)
{
	// Find plane normal
	/*
	// Find point on Ray Plane intersection
	float t = (dot(triangle.vertex0, triangle.normal) - dot(ray.origin, triangle.normal)) / dot(ray.direction, triangle.normal);
	
	// Check if point of Ray Plane intersection is within triangle
	vec3 c = (triangle.vertex1 - triangle.vertex0);
	vec3 b = (triangle.vertex2 - triangle.vertex0);
	vec3 h = ((ray.origin + (ray.direction * t)) - triangle.vertex0);

	float beta = ((b.x * h.y) - (b.y * h.x)) / ((b.x * c.y) - (b.y * c.x));
	float gamma = ((h.x * c.y) - (h.y * c.x)) / ((b.x * c.y) - (b.y * c.x));

	if (beta < 0 || gamma < 0 || beta + gamma > 1)
	{
		t = INFINITY;
		return t;
	}
	*/
	
	vec3 e1 = triangle.vertex1 - triangle.vertex0;
	vec3 e2 = triangle.vertex2 - triangle.vertex0;
	vec3 s = ray.origin - triangle.vertex0;
	vec3 s1 = cross(ray.direction, e2);
	vec3 s2 = cross(s, e1);

	vec3 tBetaGamma = (1 / dot(s1, e1)) * vec3(dot(s2, e2), dot(s1, s), dot(s2, ray.direction));

	float t = tBetaGamma.x;
	float beta = tBetaGamma.y;
	float gamma = tBetaGamma.z;
	float alpha = 1 - beta - gamma;
	

	vec3 q = ray.origin + (t * ray.direction);
	
	vec3 bMinusA = triangle.vertex1 - triangle.vertex0;
	vec3 cMinusA = triangle.vertex2 - triangle.vertex0;
	
	vec3 cMinusB = triangle.vertex2 - triangle.vertex1;
	vec3 qMinusB = q - triangle.vertex1;
	
	vec3 aMinusC = triangle.vertex0 - triangle.vertex2;
	vec3 qMinusC = q - triangle.vertex2;

	vec3 qMinusA = q - triangle.vertex0;

	float denominator = (dot(cross(bMinusA, cMinusA), triangle.normalA));

	float alpha2 = (dot(cross(cMinusB, qMinusB), triangle.normalA)) / denominator;
	float beta2 = (dot(cross(aMinusC, qMinusC), triangle.normalA)) / denominator;
	float gamma2 = (dot(cross(bMinusA, qMinusA), triangle.normalA)) / denominator;

	if (beta2 < 0 || gamma2 < 0 || beta2 + gamma2 > 1 || beta2 > 1 || gamma2 > 1)
	{
		vec3 nullT(INFINITY, 0, 0);
		return nullT;
	}
	
	vec3 tBetaGamma2(t, beta, gamma);

	triangle.SetNormal(normalize((alpha2 * triangle.normalA) + (beta2 * triangle.normalB) + (gamma2 * triangle.normalC)));
	return tBetaGamma2;
}

Intersection FindIntersection(Scene scene, Ray ray)
{
	// For each sphere in scene and for each triangle in scene
	// i < scene->spheres.length() test for intersection
	// call sphereIntersection and triangleIntersection respectively
	// return closest intersection
	
	float minDist = INFINITY;
	Sphere hitSphere;
	Triangle hitTriangle;
	bool didHit = false;
	float t = minDist;
	float tSphere = minDist;
	float tTriangle = minDist;
	vec3 hitObjectDiffuse = vec3(0, 0, 0);
	vec3 hitObjectSpecular = vec3(0, 0, 0);
	vec3 hitObjectEmission = vec3(0, 0, 0);
	vec3 hitObjectAmbient = vec3(0, 0, 0);
	vec3 hitObjectNormal = vec3(0, 0, 0);
	float hitObjectShininess = 0.0f;
	vec3 tBetaGamma2;
	vec3 intersectionPoint(INFINITY, INFINITY, INFINITY);
	vec3 center(0, 0, 0);
	bool hitObjectIsSphere = false;
	vec3 reflectionNormal(0, 0, 0);

	for (int i = 0; i < scene.spheres.size(); i++)
	{
		tSphere = CheckSphereIntersection(scene.spheres[i], ray);

		if (tSphere < minDist && tSphere > 0)
		{
			hitObjectDiffuse = scene.spheres[i].diffuse;
			hitObjectSpecular = scene.spheres[i].specular;
			hitObjectEmission = scene.spheres[i].emission;
			hitObjectAmbient = scene.spheres[i].ambient;
			hitObjectShininess = scene.spheres[i].shininess;

			hitObjectNormal = normalize(((ray.origin + (normalize(ray.direction) * tSphere)) - scene.spheres[i].center) / scene.spheres[i].radius);
		//	hitObjectNormal.x =  -(ray.origin.x + (normalize(ray.direction).x * tSphere)) - scene.spheres[i].center.x;
		//	hitObjectNormal = ((ray.origin + (ray.direction * tSphere)) - scene.spheres[i].center) / scene.spheres[i].radius;

		//	intersectionPoint = ray.origin + (normalize(ray.direction) * tSphere) ;
			intersectionPoint = ray.origin + (ray.direction * tSphere);
		//	intersectionPoint = ray.origin + (ray.direction * tSphere) + ray.origin;

		//	hitObjectNormal.x =  (ray.origin.x + (normalize(ray.direction).x * tSphere)) - scene.spheres[i].center.x + ray.origin.x;
		//	hitObjectNormal.y =  (ray.origin.y + (normalize(ray.direction).y * tSphere)) - scene.spheres[i].center.y + ray.origin.y;
		//	hitObjectNormal.z =  (ray.origin.z + (normalize(ray.direction).z * tSphere)) - scene.spheres[i].center.z + ray.origin.z;
			
		//	hitObjectNormal.x =  ((ray.origin.x + (normalize(ray.direction).x * tSphere)) - scene.spheres[i].center.x) / scene.spheres[i].radius;
		//	hitObjectNormal.y =  ((ray.origin.y + (normalize(ray.direction).y * tSphere)) - scene.spheres[i].center.y) / scene.spheres[i].radius;
		//	hitObjectNormal.z =  ((ray.origin.z + (normalize(ray.direction).z * tSphere)) - scene.spheres[i].center.z) / scene.spheres[i].radius;
			
			//hitObjectNormal.y = -hitObjectNormal.y;
			
			// check if ray dot normal is > 0 
			// if true ray is inside sphere, set normal to inverse normal
			// if false ray is outside sphere, set normal to normal;
			if (dot(ray.direction, hitObjectNormal) > 0)
			{
			//	hitObjectNormal = -hitObjectNormal;
			}
		//	hitObjectNormal.x = -hitObjectNormal.x;
		//	hitObjectNormal.z = -hitObjectNormal.z;
		//	hitObjectNormal = normalize(hitObjectNormal);
		//  TODO
			if (hitObjectNormal.y > 0.5)
			{
				ray.origin = ray.origin + vec3(0, 0, 0);
				hitSphere = scene.spheres[i];

				hitObjectNormal = hitObjectNormal + vec3(0, 0, 0);
			}
			reflectionNormal = normalize((ray.origin + (ray.direction * tSphere)) - scene.spheres[i].center + ray.origin);
			hitObjectIsSphere = true;

			center = scene.spheres[i].center;
			intersectionPoint = ray.origin + (normalize(ray.direction) * tSphere);

			if (intersectionPoint != vec3(INFINITY, INFINITY, INFINITY))
			{
				intersectionPoint = intersectionPoint + vec3(0, 0, 0);
			}

			minDist = tSphere;
			didHit = true;
			t = tSphere;
		}

	}

	// Same loop as above for triangles
	
	for (int i = 0; i < scene.triangles.size(); i++)
	{
		tBetaGamma2 = CheckTriangleIntersection(scene.triangles[i], ray);
		tTriangle = tBetaGamma2.x;
		float beta2 = tBetaGamma2.y;
		float gamma2 = tBetaGamma2.z;
		float alpha2 = 1 - beta2 - gamma2;

		if (tTriangle < minDist && tTriangle > 0)
		{
			hitObjectDiffuse = scene.triangles[i].diffuse;
			hitObjectSpecular = scene.triangles[i].specular;
			hitObjectEmission = scene.triangles[i].emission;
			hitObjectAmbient = scene.triangles[i].ambient;
			hitObjectShininess = scene.triangles[i].shininess;
			hitObjectIsSphere = false;


			hitObjectNormal = normalize((alpha2 * scene.triangles[i].normalA) + (beta2 * scene.triangles[i].normalB) + (gamma2 * scene.triangles[i].normalC));

			intersectionPoint = ray.origin + ray.direction * tTriangle;

			minDist = tTriangle;
			didHit = true;
			t = tTriangle;

		}
	}

	return Intersection(didHit, intersectionPoint, hitObjectDiffuse, hitObjectSpecular, hitObjectEmission, hitObjectShininess, hitObjectAmbient, hitObjectNormal, hitObjectIsSphere, center, reflectionNormal);
}

Ray ShootMirrorRay(Intersection intersection, Ray ray)
{
	vec3 direction = normalize(normalize(ray.direction) - ((2 * dot(normalize(ray.direction), intersection.hitObjectNormal) * intersection.hitObjectNormal)));
	vec3 origin = intersection.intersectionPoint + (direction * (0.00001f));

	Ray mirrorRay(origin, direction);
	return mirrorRay;
}

Ray ShootShadowRay(Intersection intersection, vec3 lightDirection)
{
	vec3 direction =  normalize(- lightDirection);
	vec3 origin = intersection.intersectionPoint + (direction * (0.00001f));

	Ray ray(origin, direction);
	return ray;
}

vec3 ComputePointLighting(Intersection intersection, PointLight light, Camera camera)
{
	vec3 finalColour;
	vec3 eyeDirection = normalize(camera.eyePos - intersection.intersectionPoint);
	vec3 directionPoint = normalize(light.position - intersection.intersectionPoint);

	float nDotL = dot(intersection.hitObjectNormal, directionPoint);
	// No attenuation for now
	vec3 lambert = intersection.hitObjectDiffuse * light.colour * max(nDotL, 0.0f);

	vec3 halfVec = normalize(light.position + eyeDirection);
	float nDotH = dot(intersection.hitObjectNormal, halfVec);
	vec3 phong = intersection.hitObjectSpecular * light.colour * pow(max(nDotH, 0.0f), intersection.hitObjectShininess);

	return finalColour = lambert + phong;
}


vec3 ComputeDirectionalLighting(Intersection intersection, DirectionalLight light, Camera camera, Scene scene)
{
	vec3 finalColour;
	vec3 eyeDirection = normalize(camera.eyePos - intersection.intersectionPoint);
	vec3 normalizedLightDirection = normalize(light.direction);

	float nDotL = dot(intersection.hitObjectNormal, normalizedLightDirection);
	// No attenuation for now
	vec3 lambert = intersection.hitObjectDiffuse * light.colour * max(nDotL, 0.0f);

	vec3 halfVec = normalize(normalizedLightDirection + eyeDirection);
	float nDotH = dot(intersection.hitObjectNormal, halfVec);
	vec3 phong = intersection.hitObjectSpecular * light.colour * pow(max(nDotH, 0.0f), intersection.hitObjectShininess);

	return finalColour = lambert + phong;
}

vec3 FindColour(Intersection intersection, Scene scene, Camera camera, Ray mirrorRay, Ray primaryRay)
{
	if (intersection.didHit == true)
	{
		vec3 finalColour(0, 0, 0); 
		vec3 col1(0, 0, 0);
		vec3 col2(0, 0, 0);

		for (int i = 0; i < scene.pointLights.size(); i++)
		{
			Ray ray = ShootShadowRay(intersection, -normalize(scene.pointLights[i].position - intersection.intersectionPoint));
			Intersection shadowIntersection = FindIntersection(scene, ray);
			if (shadowIntersection.didHit != true)
			{
				vec3 colour = ComputePointLighting(intersection, scene.pointLights[i], camera);
				col1 = col1 + colour;
			}
		}

		for (int i = 0; i < scene.dirLights.size(); i++)
		{
			Ray ray = ShootShadowRay(intersection, scene.dirLights[i].direction);
			Intersection shadowIntersection = FindIntersection(scene, ray);
		//	TODO FIX SHADOW RAYS
		//	if (shadowIntersection.didHit != true)
		//	{
				vec3 colour = ComputeDirectionalLighting(intersection, scene.dirLights[i], camera, scene);
				col2 = col2 + colour;
		//	}

		}

		finalColour = col1 + col2 + intersection.hitObjectAmbient + intersection.hitObjectEmission;

		if (mirrorRay.bounces < 5 && intersection.hitObjectSpecular != vec3(0, 0 ,0))
		{
			Ray previousMirrorRay;
			if (mirrorRay.bounces > 1)
			{
				 previousMirrorRay = mirrorRay;
			}
			else
			{
				previousMirrorRay = primaryRay;
			}

			vec3 tempNormal = intersection.hitObjectNormal;

			if (dot(previousMirrorRay.direction, intersection.hitObjectNormal) > 0)
			{
				tempNormal = -tempNormal;
				//	mirrorRay.direction = -mirrorRay.direction;
			}

			Intersection mirrorIntersection = FindIntersection(scene, mirrorRay);


			mirrorRay.direction.x = normalize(previousMirrorRay.direction).x - (2.0f * dot(normalize(previousMirrorRay.direction).x, tempNormal.x) * tempNormal.x);
			mirrorRay.direction.y = normalize(previousMirrorRay.direction).y - (2.0f * dot(normalize(previousMirrorRay.direction).y, tempNormal.y) * tempNormal.y);
			mirrorRay.direction.z = normalize(previousMirrorRay.direction).z - (2.0f * dot(normalize(previousMirrorRay.direction).z, tempNormal.z) * tempNormal.z);
			mirrorRay.direction = normalize(mirrorRay.direction);

			
			if (dot(mirrorRay.direction, intersection.hitObjectNormal) > 0)
			{
				mirrorRay.origin = intersection.intersectionPoint + (intersection.hitObjectNormal * 0.00001f);
			}
			else
			{
				mirrorRay.origin = intersection.intersectionPoint - (intersection.hitObjectNormal * 0.00001f);
			}

		//	Intersection mirrorIntersection = FindIntersection(scene, mirrorRay);
			vec3 mirrorColour;
			mirrorRay.bounces = mirrorRay.bounces + 1;
			vec3 reflectivity;			
		
			vec3 tempNormal2(intersection.hitObjectNormal.x, -intersection.hitObjectNormal.y, -intersection.hitObjectNormal.z);
			if (mirrorIntersection.didHit && intersection.didHit && dot(mirrorRay.direction, intersection.hitObjectNormal) < 0)
			{
				mirrorColour = FindColour(mirrorIntersection, scene, camera, mirrorRay, previousMirrorRay);
				reflectivity =  intersection.hitObjectSpecular * mirrorColour;

			}
			else
			{
 				reflectivity = intersection.hitObjectSpecular * vec3(1, 0, 0);
			}

			finalColour = reflectivity;
		}

		return finalColour;

	}
	else
	{
		return vec3(0, 0.33, 1);
	}
}

int main()
{
	const int width = 1280;
	const int height = 960;
	unsigned char* pixels = new unsigned char [width * height * 3];
	std::string outputFilename = "Raytracer.png";
	vec3 eyePosition;

	if (viewMode == "sphere")
	{
		eyePosition = vec3(0, 0, 15);
	}
	if (viewMode == "cube")
	{
		eyePosition = vec3(3, 0, 3);
	}

	vec3 center = vec3(0.0f, 0.0f, 0.0f);
	vec3 up = vec3(0, 1, 0);
	float fovY = radians(60.0f);
	// Create new Camera with default values 
	Camera camera(eyePosition, center, up, fovY);

	// Create new Scene and add Sphere and then Triangle
	Scene scene;

	if (viewMode == "sphere")
	{
		Sphere sphere0(vec3(0.00f, 0.00f, 0.00f), 1.0f, vec3(1.0f, 1.0f, 0.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.15f, 0.05f, 0.0f), 10.0f, vec3(0.1, 0.1, 0.1));
		Sphere sphere1(vec3(-3.00f, 0.00f, 0.00f), 1.0f, vec3(0.00f, 0.5, 0.5), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 10.0f, vec3(0.1, 0.1, 0.1));
		Sphere sphere2(vec3(3.00f, 0.00f, 0.00f), 1.00f, vec3(0.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 10.0f, vec3(0.1, 0.1, 0.1));
		Sphere sphere3(vec3(0.00f, 0.00f, -1.50f), 1.00f, vec3(0.0f, 1.0f, 0.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 10.0f, vec3(0.1, 0.1, 0.1));
		Sphere sphere4(vec3(0.00f, 0.00f, 1.50f), 1.00f, vec3(0.0f, 1.0f, 0.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 10.0f, vec3(0.1, 0.1, 0.1));
		Sphere sphere5(vec3(0.00f, 3.00f, 0.00f), 1.00f, vec3(1.0f, 1.0f, 0.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 10.0f, vec3(0.1, 0.1, 0.1));
		Sphere sphere6(vec3(0.00f, -3.00f, 0.00f), 1.00f, vec3(1.0f, 0.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 10.0f, vec3(0.1, 0.1, 0.1));
		
		scene.spheres.push_back(sphere0);
		scene.spheres.push_back(sphere1);
		scene.spheres.push_back(sphere2);
	//	scene.spheres.push_back(sphere3);
	//	scene.spheres.push_back(sphere4);
		scene.spheres.push_back(sphere5);
		scene.spheres.push_back(sphere6);
		
		DirectionalLight lightDir8(vec3(-1, 1, 1), vec3(0.5f, 0.5f, 0.5f));
		scene.dirLights.push_back(lightDir8);
		
		DirectionalLight lightDir9(vec3(0, 0, 1), vec3(0.5f, 0.5f, 0.5f));
	//	scene.dirLights.push_back(lightDir8);

		PointLight spherePoint0(vec3(-4, 0, -4), vec3(0.0f, 0.6f, 0.7f));
	//	scene.pointLights.push_back(spherePoint0);
		DirectionalLight lightDir2(vec3(1, 0, 1), vec3(0.0f, 0.6f, 0.7f));
	//	scene.dirLights.push_back(lightDir2);

		float triWidth = 10.0f;
		float triHeight = 10.00f;
		float triDepth = 20.0f;
		float triCenter = 0.00f;

		vec3 vert0(-triWidth , -triHeight , -triDepth );
		vec3 vert1(-triWidth , +triHeight , -triDepth );
		vec3 vert2(+triWidth , +triHeight , -triDepth );
		vec3 vert3(+triWidth , -triHeight , -triDepth );
		vec3 vert4(-triWidth , -triHeight , +triDepth );
		vec3 vert5(-triWidth , +triHeight , +triDepth );
		vec3 vert6(+triWidth , +triHeight , +triDepth );
		vec3 vert7(+triWidth , -triHeight , +triDepth );

		/*
		Triangle tri0(vert0, vert4, vert7, vec3(1.0f, 0.0f, 0.0f), vec3(1.00, 1.00, 1.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1, 0.1, 0.1));
		Triangle tri1(vert0, vert7, vert3, vec3(1.0f, 0.0f, 0.f), vec3(1.00, 1.00, 1.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri2(vert1, vert5, vert6, vec3(0.0f, 1.0f, 0.f), vec3(1.00, 1.00, 1.00f), vec3(1.0f, 1.00f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri3(vert1, vert6, vert2, vec3(0.0f, 1.0f, 0.f), vec3(1.00, 1.00, 1.00f), vec3(1.00f, 1.00f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri4(vert3, vert2, vert6, vec3(0.0f, 0.0f, 1.f), vec3(1.00, 1.00, 1.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri5(vert3, vert6, vert7, vec3(0.0f, 0.0f, 1.f), vec3(1.00, 1.00, 1.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri6(vert0, vert5, vert1, vec3(1.0f, 1.0f, 0.f), vec3(1.00, 1.00, 1.00f), vec3(0.00f, 1.00f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri7(vert0, vert4, vert5, vec3(1.0f, 1.0f, 0.f), vec3(1.00, 1.00, 1.00f), vec3(0.00f, 1.00f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri8(vert0, vert1, vert2, vec3(0.0f, 1.0f, 0.f), vec3(1.00, 1.00, 1.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri9(vert0, vert2, vert3, vec3(0.0f, 1.0f, 0.f), vec3(1.00, 1.00, 1.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri10(vert4, vert7, vert6, vec3(0.0f, 1.0f, 1.f), vec3(1.00, 1.00, 1.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri11(vert4, vert6, vert5, vec3(0.0f, 1.0f, 1.f), vec3(1.00, 1.00, 1.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		*/

		Triangle tri0(vert0, vert4, vert7, vec3(1.0f, 0.0f, 0.0f), vec3(0.00, 0.00, 0.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1, 0.1, 0.1));
		Triangle tri1(vert0, vert7, vert3, vec3(1.0f, 0.0f, 0.f), vec3(0.00, 0.00, 0.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri2(vert1, vert5, vert6, vec3(0.0f, 1.0f, 0.f), vec3(0.00, 0.00, 0.00f), vec3(1.0f, 1.00f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri3(vert1, vert6, vert2, vec3(0.0f, 1.0f, 0.f), vec3(0.00, 0.00, 0.00f), vec3(1.00f, 1.00f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri4(vert3, vert2, vert6, vec3(0.0f, 0.0f, 1.f), vec3(0.00, 0.00, 0.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri5(vert3, vert6, vert7, vec3(0.0f, 0.0f, 1.f), vec3(0.00, 0.00, 0.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri6(vert0, vert5, vert1, vec3(1.0f, 1.0f, 0.f), vec3(0.00, 0.00, 0.00f), vec3(0.00f, 1.00f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri7(vert0, vert4, vert5, vec3(1.0f, 1.0f, 0.f), vec3(0.00, 0.00, 0.00f), vec3(0.00f, 1.00f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri8(vert0, vert1, vert2, vec3(0.0f, 1.0f, 0.f), vec3(0.00, 0.00, 0.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri9(vert0, vert2, vert3, vec3(0.0f, 1.0f, 0.f), vec3(0.00, 0.00, 0.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri10(vert4, vert7, vert6, vec3(0.0f, 1.0f, 1.f), vec3(0.00, 0.00, 0.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));
		Triangle tri11(vert4, vert6, vert5, vec3(0.0f, 1.0f, 1.f), vec3(0.00, 0.00, 0.00f), vec3(0.15f, 0.05f, 0.0f), 1.00f, vec3(0.1f, 0.1f, 0.1f));

		
		// -Y
		scene.triangles.push_back(tri0);
		scene.triangles.push_back(tri1);

		// +Y
		scene.triangles.push_back(tri2);
		scene.triangles.push_back(tri3);

		// +X
		scene.triangles.push_back(tri4);
		scene.triangles.push_back(tri5);

		// -X
		scene.triangles.push_back(tri6);
		scene.triangles.push_back(tri7);
		// -Z
		scene.triangles.push_back(tri8);
		scene.triangles.push_back(tri9);

		// +Z
	//	scene.triangles.push_back(tri10);
	//	scene.triangles.push_back(tri11);
		
	}

	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			// Shoot Ray
			Ray ray = ShootRay(camera, i, j, width, height);
			Intersection intersection = FindIntersection(scene, ray);
		//	Ray mirrorRay = ShootMirrorRay(intersection, ray);
			Ray mirrorRay;

		//TODO
		//	vec3 tempNormal = intersection.reflectionNormal;
			vec3 tempNormal = intersection.hitObjectNormal;
		//	tempNormal.y = -intersection.hitObjectNormal.y;
		//	tempNormal.x = -intersection.hitObjectNormal.x;
		//	tempNormal.z = intersection.hitObjectNormal.z;
		
			// TODO, check when this code executes
			if (dot(ray.direction, intersection.hitObjectNormal) > 0)
			{
			//	tempNormal = -tempNormal;
			}

		//	tempNormal.x = 1;
		//	tempNormal.y = -1;
		//	tempNormal.z = 1;
			Ray tempRay(ray.origin, ray.direction);
			tempRay.direction = normalize(tempRay.direction);
			if (intersection.didHit)
			{
			//	mirrorRay.direction = normalize(normalize(tempRay.direction) - (2.0f * dot(normalize(tempRay.direction), tempNormal) * tempNormal));
				mirrorRay.direction.x = tempRay.direction.x - (2.0f * dot(tempRay.direction.x, tempNormal.x) * tempNormal.x);
				mirrorRay.direction.y = tempRay.direction.y - (2.0f * dot(tempRay.direction.y, tempNormal.y) * tempNormal.y);
				mirrorRay.direction.z = tempRay.direction.z - (2.0f * dot(tempRay.direction.z, tempNormal.z) * tempNormal.z);
			//	mirrorRay.direction.x = 1;
			//	mirrorRay.direction.y = 1;
			//	mirrorRay.direction.y = -mirrorRay.direction.y;
			//	mirrorRay.direction.z = 0;
				// TODO this section
			//	mirrorRay.direction.z = -mirrorRay.direction.z;
				mirrorRay.direction = normalize(mirrorRay.direction);
				if (intersection.didHit == true)
				{
					//	TODO Wrap everything in this condition
					tempNormal = tempNormal + vec3(0, 0, 0);
				}
			//	mirrorRay.direction = -mirrorRay.direction;
				if (dot(mirrorRay.direction, intersection.hitObjectNormal) > 0)
				{
					mirrorRay.origin = intersection.intersectionPoint + (mirrorRay.direction * 0.00001f)  ;
				}
				else
				{
					mirrorRay.origin.x = intersection.intersectionPoint.x - (mirrorRay.direction.x * 0.00001f)  ;
					mirrorRay.origin.y = intersection.intersectionPoint.y - (mirrorRay.direction.y * 0.00001f)  ;
					mirrorRay.origin.z = intersection.intersectionPoint.z - (mirrorRay.direction.z * 0.00001f)  ;
					// TODO Component by component maybe add or minus y opposite to what is done for x and z.. 
					// Also Maybe add Ray.origin

				}
				if (intersection.hitObjectIsSphere)
				{
					mirrorRay.origin = mirrorRay.origin + vec3(0, 0, 0);
					tempRay.origin = tempRay.origin + vec3(0, 0, 0);
				}
			}

			vec3 colour = FindColour(intersection, scene, camera, mirrorRay, ray);
			unsigned char col[3] = { 0, 0, 0 };
			// Reverse order of colours as FreeImage produces BGR image
			col[0] = unsigned char(min(colour[2] * 255, 255.0f));
			col[1] = unsigned char(min(colour[1] * 255, 255.0f));
			col[2] = unsigned char(min(colour[0] * 255, 255.0f));
			memcpy(&pixels[((j * width) + i)  * 3], &col, 3);
		}
	}

	FreeImage_Initialise();
	FIBITMAP* img = FreeImage_ConvertFromRawBits((BYTE*)pixels, width, height, 3 * width, 24, 0xFF0000, 0x00FF00, 0x0000FF, false);
	FreeImage_Save(FIF_PNG, img, outputFilename.c_str(), 0);
	FreeImage_DeInitialise();
}
