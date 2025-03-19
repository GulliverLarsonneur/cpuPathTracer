#define _CRT_SECURE_NO_WARNINGS 1

// Performance settings :

#define IMAGE_HEIGHT 128
#define IMAGE_WIDTH (int)(IMAGE_HEIGHT * 1.37)
#define RENDER_SAMPLES 1
#define GLOBAL_NUM_BOUNCES 3

// Aesthetics settings :
#define CAMERA_FOV 80.0
#define CAMERA_VERTICAL_ANGLE 0.0
#define CAMERA_LATERAL_ANGLE 10.0
#define CAMERA_X 15.0
#define CAMERA_Y 0.0
#define CAMERA_Z 45.0
#define LIGHT_INTENSITY 500 * 2e07
#define AIR_INDEX 1.0
#define AA_RADIUS 0.2
#define DEPTH_OF_FIELD_AMPLITUDE 0.1
#define DEPTH_OF_FIELD_DISTANCE 55.0

// Toggleable optimizations and techniques :
#define EXTENDED_SOURCE 0
#define ACTIVATE_ANTIALIASING 1
#define ACTIVATE_DEPTH_OF_FIELD 1
#define ACTIVATE_INDIRECT_LIGHTING 1
#define ACTIVATE_IMPORTANCE_SAMPLING 1
#define ACTIVATE_CUSTOM_MATERIALS 1
#define SPHERE_BBOX_OPTIMISATION 1
#define MESH_BVH_OPTIMIZATION 1
#define MESH_AABB_OPTIMIZATION 1

// Globals
#define M_PI 3.1415926535897932384626
#define EPSILON 0.00001
#define _OPENMP

#include <iostream>
#include <vector>
#include <chrono>
#include <ctime>

#include <random>
#ifdef _OPENMP
	#include <omp.h>
	static const int threadCount = omp_get_max_threads();
#else
	static const int threadCount = 1;
#endif
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "utilities.hpp"



class Scene
{
public:
	bool intersect(const Ray& ray, Vector3& intersectionPoint, Vector3& intersectionNormal, double& t, int& hitObjectIndex, Vector3& object_col) // objectIndex is a returned int to know th albedo
	{
		bool hasIntersected = false;

		double t_min = std::numeric_limits<double>::max();
#if EXTENDED_SOURCE
		int startIndex = 0;
#else
		int startIndex = 1; // We ignore the light sphere
#endif

		for (int objectIndex = startIndex; objectIndex < objects.size(); objectIndex++)
		{
			Vector3 localIntersectionPoint = { 0, 0, 0 };
			Vector3 localIntersectionNormal = { 0, 0, 0 };
			Vector3 localCol = { 0, 0, 0 };
			double local_t = 0;
			if (objects[objectIndex]->intersect(ray, localIntersectionPoint, localIntersectionNormal, local_t, localCol))
			{
				if (local_t <= t_min)
				{
					//std::cout << "local_t = " << local_t << "  -  old_t = " << old_t << "\n";
					intersectionPoint = localIntersectionPoint;
					intersectionNormal = localIntersectionNormal;
					t_min = local_t;
					t = t_min;
					hitObjectIndex = objectIndex;
					object_col = localCol;
				}
				hasIntersected = true;
			}
		}
		return hasIntersected;
	}

	Vector3 getColor(Ray ray, const int numBounces)
	{
		Vector3 intersectionPoint = { 0.0, 0.0 ,0.0 };
		Vector3 intersectionNormal = { 0.0, 0.0, 0.0 };

		Vector3 color = { 0, 0, 0 };
		int hitObjectIndex = 0;
		double t = 0;
		if (intersect(ray, intersectionPoint, intersectionNormal, t, hitObjectIndex, color))
		{
			if (numBounces < 0)
			{
				return { 0, 0, 0 };
				//return { 10000000, 0, 0 };
			}

			bool isEntering = true;

			if (objects[hitObjectIndex]->isTransparent || objects[hitObjectIndex]->isMirror)
			{
				if (objects[hitObjectIndex]->isMirror)
				{
					Ray bounceRay = { intersectionPoint + EPSILON * intersectionNormal, ray.direction - 2.0 * dot(ray.direction, intersectionNormal) * intersectionNormal };
					return getColor(bounceRay, numBounces - 1);
				}

				double indexRatio = 0;
				double n1 = AIR_INDEX;
				double n2 = objects[hitObjectIndex]->refractiveIndex;
				Ray fresnelRay = { { 0, 0, 0 }, { 0, 0, 0 } };
				if (dot(ray.direction, intersectionNormal) < 0) // We enter inside the object
				{
					indexRatio = (AIR_INDEX / objects[hitObjectIndex]->refractiveIndex);
					Vector3 Tt = indexRatio * (ray.direction - dot(ray.direction, intersectionNormal) * intersectionNormal);
					Vector3 Tn = -sqrt(1 - sqr(indexRatio) * (1 - sqr(dot(ray.direction, intersectionNormal)))) * intersectionNormal;
					fresnelRay = { intersectionPoint - EPSILON * intersectionNormal, Tn + Tt };
				}
				else // else, we exit the object
				{
					std::swap(n1, n2);
					indexRatio = (objects[hitObjectIndex]->refractiveIndex / AIR_INDEX);
					Vector3 Tt = indexRatio * (ray.direction - dot(ray.direction, intersectionNormal) * intersectionNormal);
					Vector3 Tn = sqrt(1 - sqr(indexRatio) * (1 - sqr(dot(ray.direction, intersectionNormal)))) * intersectionNormal;
					fresnelRay = { intersectionPoint + EPSILON * intersectionNormal, Tn + Tt };
					isEntering = false;
				}

				double k0 = sqr(n1 - n2) / sqr(n1 + n2);
				double alphaR = k0 + (1 - k0) * pow((1 - abs(dot(intersectionNormal, ray.direction))), 5); // alphaR = Reflection coefficient

				//int threadID = omp_get_thread_num(); // TO DO : fix this !
				//std::cout << hash13(intersectionPoint) << "\n";

				if (alphaR > uniform(engine))//uniform(engine))//uniform(engine[threadID])) // hash13(intersectionPoint)
				{
					Ray bounceRay = { intersectionPoint + EPSILON * intersectionNormal, ray.direction - 2.0 * dot(ray.direction, intersectionNormal) * intersectionNormal };
					return getColor(bounceRay, numBounces - 1);
				}
				else
				{
					return getColor(fresnelRay, numBounces - 1);
				}
			}

#if !EXTENDED_SOURCE
			const Sphere* light = dynamic_cast<const Sphere *>(objects[0]);

			Vector3 pointToLight = light->center - intersectionPoint;
			double d2 = pointToLight.norm2();
			pointToLight = pointToLight / sqrt(d2);

			Ray shadowRay = { intersectionPoint + EPSILON * intersectionNormal, pointToLight };
			double shadow_t = 0;
			bool isLightened = true;
			Vector3 a = { 0, 0, 0 };
			Vector3 b = { 0, 0, 0 };
			Vector3 c = { 0, 0, 0 };
			int hitObjectIndex2 = 0;

			if (intersect(shadowRay, a, b, shadow_t, hitObjectIndex2, c))
			{
				if (hitObjectIndex2 != 0)
				{
					//std::cout << shadow_t << "\n";
					if (shadow_t * shadow_t < d2)
					{
						isLightened = false;
					}
				}
			}

			// Direct lighting
			if (isLightened)
			{
				return color * lightIntensity * dot(intersectionNormal, pointToLight) / (4 * M_PI * d2 * M_PI);
			}
			else
			{
				return { 0, 0, 0 };
			}
#else // !EXTENDED_SOURCE

#if ACTIVATE_IMPORTANCE_SAMPLING
			Vector3 bounceOrigin = { 0, 0, 0 };

			if (isEntering)
			{
				bounceOrigin = EPSILON * intersectionNormal + intersectionPoint;
			}
			else
			{
				bounceOrigin = - EPSILON * intersectionNormal + intersectionPoint;
			}

			const Sphere* light = dynamic_cast<const Sphere *>(objects[0]);
			Vector3 shadowPoint = { 0, 0, 0 };
			Vector3 shadowNormal = { 0, 0, 0 };
			Vector3 shadowCol = { 0, 0, 0 };
			Vector3 directCol = { 0, 0, 0 };

			// We will shoot a ray to tell whether there is shadow on the current intersectionPoint
			const Vector3 pointToLight = normalize(light->center - intersectionPoint);
			const Vector3 lightNormal = normalize(random_cos( -1.0 * pointToLight, intersectionPoint));
			const Vector3 lightPointTarget = light->center + lightNormal * light->radius + EPSILON * lightNormal;
			const Vector3 directionToLight = normalize(lightPointTarget - bounceOrigin);
			double distanceToLight = sqrt((lightPointTarget - bounceOrigin).norm2());

			Ray shadow_ray = { bounceOrigin, directionToLight };
#if ACTIVATE_INDIRECT_LIGHTING
			// Indirect lighting
			Vector3 indirectCol = color * getColor({ bounceOrigin, random_cos(intersectionNormal, intersectionPoint) }, numBounces - 1);
#else // ACTIVATE_INDIRECT_LIGHTING
			Vector3 indirectCol = { 0, 0, 0 };
#endif // ACTIVATE_INDIRECT_LIGHTING
			for (int i = 0; i < objects.size(); ++i) // This loop MUST include the light itself, e.g. start with 0
			{
				bool shadowRayIntersects = objects[i]->intersect(shadow_ray, shadowPoint, shadowNormal, t, shadowCol);
				if (t < distanceToLight)
				{
					return indirectCol;
				}
			}

			directCol = color * lightIntensity * light->albedo / (4 * sqr(M_PI));
			directCol = directCol * dot(directionToLight, intersectionNormal) / dot(-1 * pointToLight, lightNormal);
			directCol = directCol * dot(-1 * directionToLight, lightNormal) / sqr(distanceToLight);

			return directCol + indirectCol;
#else // ACTIVATE_IMPORTANCE_SAMPLING
			// Extended source
			if (hitObjectIndex == 0)
			{
				return objects[hitObjectIndex]->albedo * lightIntensity / (4 * sqr(M_PI * dynamic_cast<const Sphere*>(objects[hitObjectIndex])->radius));
			}

			// Indirect lighting
			Vector3 randomRayDirection = random_cos(intersectionNormal, intersectionPoint); // randomly sample ray us ing random cos
			//Vector3 indirectLightingColor = objects[hitObjectIndex]->albedo * getColor({ intersectionPoint + EPSILON * intersectionNormal, randomRayDirection }, numBounces - 1);
			//return color + indirectLightingColor;

			return color * getColor({ intersectionPoint + EPSILON * intersectionNormal, randomRayDirection }, numBounces - 1);
#endif // ACTIVATE_IMPORTANCE_SAMPLING

#endif // !EXTENDED_SOURCE
		}
		else
		{
			return { 0, 0, 0 };
		}
	}

	void addObject(const Geometry* g)
	{
		objects.push_back(g);
	}

	std::vector<const Geometry*> objects;
	double lightIntensity = LIGHT_INTENSITY;
};


int main() 
{
	Timer timer;
	double cameraFOV = CAMERA_FOV * M_PI / 180.0;
	Vector3 camOrigin( CAMERA_X, CAMERA_Y, CAMERA_Z );
	double angleUp = CAMERA_VERTICAL_ANGLE * M_PI / 180.0;
	double lateralAngle = CAMERA_LATERAL_ANGLE * M_PI / 180.0;
	Vector3 cameraUp(0.0, cos(angleUp), sin(angleUp));
	Vector3 cameraDir(0.0, -sin(angleUp), cos(angleUp));
	Vector3 cameraRight = cross(cameraUp, cameraDir);
	Vector3 rotatedCameraDir = cameraDir * cos(lateralAngle) + cameraRight * sin(lateralAngle);
	cameraDir = rotatedCameraDir;
	cameraRight = cross(cameraUp, cameraDir);
	
	Scene scene;
	//                               Center,          albedo,    radius,  MateriaType,      refractiveIndex
	// The light MUST be the first object
	scene.addObject(new Sphere({ { 0, 31, 0 },  {1.0, 1.0, 1.0}, 4.0, MaterialType::ALBEDO,  1.3 })); // Extended light !!
	scene.addObject(new Sphere({ { 1000, 0, 0 },  {0.7, 0.3, 0.3}, 940.0,  MaterialType::ALBEDO, 1.0 })); // Right red wall
	scene.addObject(new Sphere({ { -1000, 0, 0 }, {0.2, 0.2, 0.7}, 940.0,  MaterialType::ALBEDO, 1.0 })); // Left blue wall
	scene.addObject(new Sphere({ { 0, 1000, 0 },  {1.0, 1.0, 1.0}, 940.0,  MaterialType::ALBEDO, 1.0 })); // Top white ceiling
	scene.addObject(new Sphere({ { 0, -950, 0 },  {0.3, 0.7, 0.3}, 940.0,  MaterialType::COLOR_CHECKERBOARD, 1.0 })); // Bottom floor
	scene.addObject(new Sphere({ { 0, 0, 1000 },  {1.0, 1.0, 1.0}, 940.0,  MaterialType::ALBEDO, 1.0 })); // Back (behind camera) wall
	scene.addObject(new Sphere({ { 0, 0, -1000 }, {0.2, 0.2, 0.7}, 940.0,  MaterialType::CHECKERBOARD, 1.0 })); // Front blue wall - CHECKERBOARD

	scene.addObject(new Sphere({ { 26, -1, -10 },  {0.6, 0.1, 0.8}, 8.0, MaterialType::ALBEDO, 1.0 }));  // Mat sphere
	scene.addObject(new Sphere({ { 8.5, 17, -10 },   {0.1, 0.1, 0.4}, 9.0,  MaterialType::MIRROR, 1.0 }));  // Mirror sphere
	scene.addObject(new Sphere({ { 7, -7, 19 },  {0.3, 0.7, 0.3}, 3.5,  MaterialType::TRANSPARENT,  2.0 }));  // Transparent sphere 1
	scene.addObject(new Sphere({ { 16, -8.5, 28 },  {0.3, 0.7, 0.3}, 2.0,  MaterialType::TRANSPARENT,  1.2 }));  // Transparent sphere 2

	//                                    Albedo,       MaterialType,       refractiveIndex,    modelFile,     Scale,   Position,     TextureFile  
	scene.addObject(new TriangleMesh({ {1.0, 1.0, 1.0}, MaterialType::ALBEDO, 1.0, "resources/dragon/scene.obj", 0.2, { 0, -10, 0 }, "resources/dragon/textures/DefaultMaterial_baseColor.jpeg" }));
	scene.addObject(new TriangleMesh({ {1.0, 1.0, 1.0}, MaterialType::MIRROR, 1.0, "resources/suzanne/suzanne.obj", 5.5, { -10, -4, 12.5 }}));
	scene.addObject(new TriangleMesh({ {1.0, 0.6, 0.1}, MaterialType::FIRE, 1.0, "resources/teapot/teapot.obj", 3.0, { 21, -10, 12 }}));
	//scene.addObject(new TriangleMesh({ {1.0, 1.0, 1.0}, MaterialType::ALBEDO, 1.0, "resources/tree/pine_tree.obj", 0.05, { 36, 6, 0 }, "resources/tree/10447_Pine_Tree_v1_Diffuse.jpg"}));
	//scene.addObject(new TriangleMesh({ {1.0, 1.0, 1.0}, MaterialType::ALBEDO, 1.0, "resources/cat/cat.obj", 0.2, { 7, -10, 23 }, "resources/cat/cat_diff.png"}));
	
	const double gamma = 2.2;

#ifdef _OPENMP
	std::cout << "[INFO] OpenMP is active.\n";
#endif
	char* image = new char[IMAGE_WIDTH * IMAGE_HEIGHT * 3];
	double deltaX, deltaY = 0;
	double deltaXDOF,deltaYDOF = 0;
	Vector3 intersectionPoint = { 0, 0, 0 };
	Vector3 intersectionNormal = {0, 0, 0};
	std::cout << "[INFO] Rendering started.\n";
	timer.start();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
	for (int x = 0; x < IMAGE_HEIGHT; ++x) 
	{
		for (int y = 0; y < IMAGE_WIDTH; ++y) 
		{
			Vector3 color = { 0, 0, 0 };

			for (int step = 0; step < RENDER_SAMPLES; ++step )
			{
#if ACTIVATE_ANTIALIASING
				boxMuller(AA_RADIUS, deltaX, deltaY);
				Vector3 camDirection = normalize(Vector3(y + deltaY - IMAGE_WIDTH / 2, -(x + deltaX - IMAGE_HEIGHT / 2), - IMAGE_WIDTH / ( 2 * tan(cameraFOV / 2) )));
#else
				Vector3 camDirection = normalize(Vector3(y - IMAGE_WIDTH / 2, -(x - IMAGE_HEIGHT / 2), - IMAGE_WIDTH / ( 2 * tan(cameraFOV / 2) )));
#endif
#if ACTIVATE_DEPTH_OF_FIELD
				boxMuller(DEPTH_OF_FIELD_AMPLITUDE, deltaXDOF, deltaYDOF);
				Vector3 locCamOrigin = camOrigin + Vector3(deltaXDOF, deltaYDOF, 0);
				Vector3 camDestination = camOrigin + camDirection * DEPTH_OF_FIELD_DISTANCE;
				camDirection = normalize(camDestination - locCamOrigin);
				Vector3 localCamDirection = camDirection[0] * cameraRight + camDirection[1] * cameraUp + camDirection[2] * cameraDir;

				Vector3 sceneColor = scene.getColor({ locCamOrigin, localCamDirection}, GLOBAL_NUM_BOUNCES); // TO DO : The diaphragm, in theory, should not be gaussian, but circular with rejection, or hexagonal
				color = color + sceneColor;
#else
				Vector3 localCamDirection = camDirection[0] * cameraRight + camDirection[1] * cameraUp + camDirection[2] * cameraDir;
				Vector3 sceneColor = scene.getColor({ camOrigin, localCamDirection }, GLOBAL_NUM_BOUNCES);
				color = color + sceneColor;
#endif
			}
			setImageColor(image, x * IMAGE_WIDTH + y, gammaCorrect(color / RENDER_SAMPLES, gamma));
		}
	}
	timer.stop();

	std::cout << "[INFO] Image generation done.\n";
	std::cout << "[INFO] Rendering time : " << timer.elapsedMilliseconds() / 1000.0 << "s.\n";


	stbi_write_png(("outputImage/final-" 
		+ std::to_string(IMAGE_WIDTH) 
		+ "x" 
		+ std::to_string(IMAGE_HEIGHT)
		+ "_pixels-" 
		+ std::to_string(RENDER_SAMPLES)
		+ "_samples-" 
		+ std::to_string(GLOBAL_NUM_BOUNCES)
		+ "_bounces-" 
#if EXTENDED_SOURCE
		+ "EXT_SOURCE-"
#else
		+ "PT_LIGHT-"
#endif
#if MESH_BVH_OPTIMIZATION
		+ "BVH-"
#endif
#if ACTIVATE_ANTIALIASING
		+ "AA-"
#endif
#if ACTIVATE_DEPTH_OF_FIELD
		+ "DOF-"
#endif
#if ACTIVATE_INDIRECT_LIGHTING
		+ "indirect-"
#endif
#if ACTIVATE_INDIRECT_LIGHTING
		+ "indirect-"
#endif
#if ACTIVATE_IMPORTANCE_SAMPLING
		+ "importance-"
#endif
		+ std::to_string(timer.elapsedMilliseconds() / 1000.0) + "_seconds.png").c_str(), IMAGE_WIDTH, IMAGE_HEIGHT, 3, &image[0], 0);
	delete[] image;

	return 0;
}