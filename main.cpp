#define _CRT_SECURE_NO_WARNINGS 1

#define IMAGE_WIDTH 256
#define CAMERA_FOV 70.0
#define M_PI 3.1415926535897932384626
#define EPSILON 0.00001
#define AIR_INDEX 1.0
#define RENDER_STEP_COUNT 300
#define AA_RADIUS 0.1
#define GLOBAL_NUM_BOUNCES 3
#define DEPTH_OF_FIELD_AMPLITUDE 0.0001
#define DEPTH_OF_FIELD_DISTANCE 1.0
#define _OPENMP


#define ACTIVATE_IMPORTANCE_SAMPLING 1
#define SPHERE_BBOX_OPTIMISATION 1

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


double clamp255(double x)
{
	if (x > 255.0) return 255.0;
	if (x < 0.0) return 0.0;
	return x;
}


void setImageColor(char* image, int index, Vector3 color)
{
	image[index * 3 + 0] = clamp255(color[0]);
	image[index * 3 + 1] = clamp255(color[1]);
	image[index * 3 + 2] = clamp255(color[2]);
}

Vector3 gammaCorrect(Vector3 v, double gamma)
{
	return Vector3({ std::pow(v[0], 1.0 / gamma), std::pow(v[1], 1.0 / gamma), std::pow(v[2], 1.0 / gamma) });
}

class Scene
{
public:
	bool intersect(const Ray& ray, Vector3& intersectionPoint, Vector3& intersectionNormal, double& t, int& hitObjectIndex, Vector3& object_col) // objectIndex is a returned int to know th albedo
	{
		bool hasIntersected = false;

		double t_min = std::numeric_limits<double>::max();
		for (int objectIndex = 0; objectIndex < objects.size(); objectIndex++)
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

			/*
			Vector3 pointToLight = lightPosition - intersectionPoint;
			double d2 = pointToLight.norm2();
			pointToLight = pointToLight / sqrt(d2);

			Ray shadowRay = { intersectionPoint + EPSILON * intersectionNormal, pointToLight };
			double shadow_t = 0;
			bool isLight = true;
			Vector3 a = { 0, 0, 0 };
			Vector3 b = { 0, 0, 0 };
			int _ = 0;


			if (intersect(shadowRay, a, b, shadow_t, _))
			{
				//std::cout << shadow_t << "\n";
				if (shadow_t * shadow_t < d2)
				{
					isLight = false;
				}
			}
			Vector3 color = { 0, 0, 0 };
			// Direct lighting
			if (isLight)
			{
				color = objects[hitObjectIndex].albedo * lightIntensity * dot(intersectionNormal, pointToLight) / (4 * M_PI * d2 * M_PI);
			}
			*/

			/*
			Principle of the importance sampling : On the LAST bounce :
			Camera -> many bounces -> Object intersection point P -> Light object origin L

			Then, construct a ray R coming from L and going towards P, with a random angle perturbation with random cos
			intersect R with the surface of the light object, get a L2 point

			PL2 is the ray that should then be used to compute the lighting 

			ALL vectors must be normalized !!
			Do not forget epsilon

			With the final normalized vector, we compute the scalar product with the normal distance. We divide by the square of the radius, and by 0 if we are in shadow
			*/

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

			const Sphere * light = dynamic_cast<const Sphere *>(objects[0]);
			const Vector3 pointToLight = normalize(light->center - intersectionPoint);
			const Vector3 lightNormal = normalize(random_cos( -1.0 * pointToLight, intersectionPoint));

			Vector3 shadowPoint = { 0, 0, 0 };
			Vector3 shadowNormal = { 0, 0, 0 };
			Vector3 shadowCol = { 0, 0, 0 };
			double best_t = std::numeric_limits<double>::max();

			// We will shoot a ray to tell whether there is shadow on the current intersectionPoint
			const Vector3 lightPointTarget = light->center + lightNormal * light->radius + EPSILON * lightNormal;
			const Vector3 directionToLight = normalize(lightPointTarget - bounceOrigin);
			double distanceToLight = sqrt((lightPointTarget - bounceOrigin).norm2());

			Ray shadow_ray = { bounceOrigin, directionToLight };

			for (int i = 1; i < objects.size(); ++i) // We can ignore index 0 which is the light itself
			{
				bool shadow_ray_intersected = objects[i]->intersect(shadow_ray, shadowPoint, shadowNormal, t, shadowCol);
				if (t <= best_t)
				{
					best_t = t;
				}
			}

			Vector3 lightColor = lightIntensity * light->albedo;
			Vector3 directCol = { 0, 0, 0 };

			if (!(best_t <= distanceToLight))
			{
				directCol = color * lightColor / (4 * sqr(M_PI));
				directCol = directCol * dot(directionToLight, intersectionNormal) / dot(-1 * pointToLight, lightNormal);
				directCol = directCol * dot(-1 * directionToLight, lightNormal) / sqr(distanceToLight);
			}

			// Indirect lighting
			Vector3 indirectCol = color * getColor({ bounceOrigin, random_cos(intersectionNormal, intersectionPoint) }, numBounces - 1);

			return directCol + indirectCol;
#else
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
#endif
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
	Vector3 lightPosition = { 10.0, 20.0, 40.0}; // x = -10 dans le pdf de cours
	double lightIntensity = 1000 * 2e07;
};

// TO DO : implement custom hash function

int main() 
{
	// Orientation canonique : Y est vers le haut
	/*
	for (int i = 0; i < threadCount; ++i)
	{
		engine[i].seed(i);
	}
	*/
	Timer timer;
	timer.start();
	const int W = IMAGE_WIDTH;
	const int H = IMAGE_WIDTH;
	double cameraFOV = CAMERA_FOV * M_PI / 180.0; // FOV angle
	Vector3 camOrigin( 0.0, 0.0, 40.0 );
	
	double angleUp = 0.0 * M_PI / 180.0;
	Vector3 cameraUp(0.0, cos(angleUp), sin(angleUp));
	Vector3 cameraDir(0.0, -sin(angleUp), cos(angleUp));
	Vector3 cameraRight = cross(cameraUp, cameraDir);
	Scene scene;
	scene.addObject(new Sphere({ { 0, 25, 0 },  {1.0, 1.0, 1.0}, 4.0, MaterialType::ALBEDO,  1.3 })); // Extended light !!
	/*
	//                                    Albedo,       MaterialType,       refractiveIndex,    modelFile,     Scale,   Position,     TextureFile  
	scene.addObject(new TriangleMesh({ {1.0, 1.0, 1.0}, MaterialType::ALBEDO, 1.0, "resources/dragon/scene.obj", 0.2, { 0, -10, 0 }, "resources/dragon/textures/DefaultMaterial_baseColor.jpeg" }));
	scene.addObject(new TriangleMesh({ {1.0, 1.0, 1.0}, MaterialType::MIRROR, 1.0, "resources/suzanne/suzanne.obj", 7.0, { -15, -3, 0 }}));
	scene.addObject(new TriangleMesh({ {1.0, 0.6, 0.1}, MaterialType::FIRE, 1.0, "resources/teapot/teapot.obj", 3.0, { 11, -10, 8 }}));
	//scene.addObject(new TriangleMesh({ {1.0, 1.0, 1.0}, MaterialType::ALBEDO, 1.0, "resources/cat/cat.obj", 0.15, { 0, -11, 17 }, "resources/cat/cat_diff.png"}));
	
	//                               Center,          albedo,    radius,  MateriaType,      refractiveIndex
	scene.addObject(new Sphere({ { 20, 0, -10 },  {0.6, 0.1, 0.8}, 10.0, MaterialType::ALBEDO, 1.0 }));  // Mat sphere 1
	scene.addObject(new Sphere({ { 12, 12, 0 },   {0.1, 0.1, 0.4}, 5.0,  MaterialType::MIRROR, 1.0 }));  // Mirror sphere
	scene.addObject(new Sphere({ { -12, 12, 0 },  {0.3, 0.7, 0.3}, 5.0,  MaterialType::TRANSPARENT,  2.0 }));  // Transparent sphere
	scene.addObject(new Sphere({ { -5, -3, 12 },  {0.1, 0.1, 0.4}, 3.0,  MaterialType::MIRROR, 1.3 }));   // Mirror sphere 2

	*/

	scene.addObject(new Sphere({ { 1000, 0, 0 },  {0.7, 0.3, 0.3}, 940.0,  MaterialType::ALBEDO, 1.0 })); // Right red wall
	scene.addObject(new Sphere({ { -1000, 0, 0 }, {0.2, 0.2, 0.7}, 940.0,  MaterialType::ALBEDO, 1.0 })); // Left blue wall
	scene.addObject(new Sphere({ { 0, 1000, 0 },  {1.0, 1.0, 1.0}, 940.0,  MaterialType::ALBEDO, 1.0 })); // Top white ceiling
	scene.addObject(new Sphere({ { 0, -950, 0 },  {0.3, 0.7, 0.3}, 940.0,  MaterialType::COLOR_CHECKERBOARD, 1.0 })); // Bottom floor
	scene.addObject(new Sphere({ { 0, 0, 1000 },  {1.0, 1.0, 1.0}, 940.0,  MaterialType::ALBEDO, 1.0 })); // Back (behind camera) wall
	scene.addObject(new Sphere({ { 0, 0, -1000 }, {0.2, 0.2, 0.7}, 940.0,  MaterialType::CHECKERBOARD, 1.0 })); // Front blue wall - CHECKERBOARD
	const double gamma = 2.2;

#ifdef _OPENMP
	std::cout << "[INFO] OpenMP is active.\n";
#endif

	char* image = new char[W * H * 3];
	std::cout << "[INFO] Rendering started.\n";
#pragma omp parallel for schedule(dynamic, 1)
	for (int x = 0; x < H; ++x) 
	{
		for (int y = 0; y < W; ++y) 
		{
			Vector3 color = { 0, 0, 0 };

			for (int step = 0; step < RENDER_STEP_COUNT; ++step )
			{
				double deltaX = 0;
				double deltaY = 0;
				boxMuller(AA_RADIUS, deltaX, deltaY);
				Vector3 u = normalize(Vector3(y + deltaY - W / 2, -(x + deltaX - H / 2), - W / ( 2 * tan(cameraFOV / 2) )));
				Vector3 intersectionPoint = {0, 0, 0}, intersectionNormal = {0, 0, 0};
				double deltaXDOF;
				double deltaYDOF;
				boxMuller(DEPTH_OF_FIELD_AMPLITUDE, deltaXDOF, deltaYDOF);
				Vector3 locCamOrigin = camOrigin + Vector3(deltaXDOF, deltaYDOF, 0);
				Vector3 camDestination = camOrigin + u * DEPTH_OF_FIELD_DISTANCE;
				Vector3 localCamDirection = normalize(camDestination - locCamOrigin);

				localCamDirection = localCamDirection[0] * cameraRight + localCamDirection[1] * cameraUp + localCamDirection[2] * cameraDir;

				color += scene.getColor({ locCamOrigin, localCamDirection}, GLOBAL_NUM_BOUNCES); // TO DO : The diaphragm, in theory, should not be gaussian, but circular with rejection, or hexagonal
			}
			setImageColor(image, x * W + y, gammaCorrect(color / RENDER_STEP_COUNT, gamma));
		}
	}
	timer.stop();
	stbi_write_png("outputImage/final.png", W, H, 3, &image[0], 0);
	delete[] image;
	std::cout << "[INFO] Image generation done.\n";
	std::cout << "[INFO] Rendering time : " << timer.elapsedMilliseconds() / 1000.0 << "s.\n";
	return 0;
}