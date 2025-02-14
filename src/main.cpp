#define _CRT_SECURE_NO_WARNINGS 1


// TO DO : put all the parameters in a static struct
#define IMAGE_WIDTH 128
#define M_PI 3.1415926535897932384626
#define EPSILON 0.00001
#define AIR_INDEX 1.0
#define RENDER_STEP_COUNT 120
#define AA_RADIUS 0.1 // TO DO : fix cross error when AA_RADIUS is set to 0 !!!
#define GLOBAL_NUM_BOUNCES 4
#define DEPTH_OF_FIELD_AMPLITUDE 0.000
#define DEPTH_OF_FIELD_DISTANCE 1.0
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

// TO DO : make ALL threads work on different parts of the image ! (otherwise, omp parallel for makes some threads work more than some other ones)
// TO DO : shift index with wavelengts -> Shoot rays of a certain wavelength
// TO DO : implement importance sampling (07/02/2025) for extended sources
// TO DO : Make macro to enable direct / extended sources lighting
// TO DO : make render progression percentage print
// TO DO : implement on GPU ?
// TO DO : allow other camera ratios
// TO DO : ALWAYS compare the std::uniform and custom hash output images !
// TO DO : Add a #define STOCHASTIC_SAMPLING 1
// TO DO : make multiple different random engines !
// TO DO : Implement Sobol and Halton sequences Quasi Monte-Carlo sampling
// TO DO : implement "ray killing" feature : the more bounces and the largest the distance is, the less bounces we compute
// TO DO : make it possible to input a transform matrix for the spheres
// TO DO : implement custom ratio
// 
// Methode de Sobol : polynômes irréductibles x^... + x^3 + x^2 + x + 1
// Use random sequences scrambling ! Owen scrambling 
// Sobol est assez imbattable, Halton est moins bien
// Autre technique qui marche très bien : choisir un point de départ et un vecteur d'offset bien choisi, que l'on itère, ça boucle sur les bords
// 	   -> Rank 1 sequencing
// 	   -> des ensembles de vecteurs qui marchent bien ont été trouvés dans toutes les dimensions
// Le bruit bleu, ça marche très bien en général : certains disent que Poisson c'est du blue noise, d'autre disent que pas
// 	   Plein de façons de le générer, avec du transport optimal, relaxation de points, GBN (gaussian blue noise), lent mais fait des points vraiment bien distribués
// L'analyse de ces points se fait avec le spectre moyen en amplitude de la transformée de Fourrier des coordonnées des points.
// Comme c'est à peu près isotrope selon la fréquence 0, 
// et quand on plot une tranche 1D de cette transformée de Fourier (elle est avec une symétrie de révolution autour de (0,0), 
// on a une caractéristique où on veut des valeurs les plus proches de zéro possible dans les basses fréquences, 
// avec une fréquence non nulle le plus à droite possible (excepté le pic précisément situé en (0,0) )
// Problème des grilles : marche très mal en dimension élevée, à cause du manque de contrôle sur les nombres de poins. 
// Les grilles c'est pas mal quand même, de base, mais Sobol c'est mieux.
// 	Pour implémenter le Owen Scrambling, on ne veut pas stocker tout l'arbre des permutations, 
// donc on a besoin d'une fonction de hashage qui dit si on permute en chaque point,
// 	et c'est assez compliqué.
// 	   C'est assez trivial de faire les méthodes de rang 1. Il y a des tables de vecteurs d'offset possibles.
// Autre nom des méthodes de rang 1 : random latices -> Voir publication de Pierre l'Ecuyer.
// Ces low discrepency sequences sont utilisées à chaque rebond
// 	Le Owen Scrambling étant assez compliqué, on peut aussi translater d'un vecteur aléatoire tous les points du pixel, avec un modulo sur les coordonnées
// static std::vector<std::default_random_engine> engine(threadCount); // random seed = 10
//static std::uniform_real_distribution<double> uniform( 0 , 1) ;
// Note : l'encodage des polynomes est en décimal pour représenter le binaire, et le coefficient de plus bas degré (1) et de plus haut degré (1 aussi) sont ommis





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
	bool intersect(const Ray& ray, Vector3& intersectionPoint, Vector3& intersectionNormal, double& t, int& hitObjectIndex) // objectIndex is a returned int to know th albedo
	{
		bool hasIntersected = false;

		double t_min = std::numeric_limits<double>::max();
		for (int objectIndex = 0; objectIndex < objects.size(); objectIndex++)
		{
			Vector3 localIntersectionPoint = { 0, 0, 0 };
			Vector3 localIntersectionNormal = { 0, 0, 0 };
			double local_t = 0;
			if (objects[objectIndex]->intersect(ray, localIntersectionPoint, localIntersectionNormal, local_t))
			{
				if (local_t <= t_min)
				{
					//std::cout << "local_t = " << local_t << "  -  old_t = " << old_t << "\n";
					intersectionPoint = localIntersectionPoint;
					intersectionNormal = localIntersectionNormal;
					t_min = local_t;
					t = t_min;
					hitObjectIndex = objectIndex;
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
		int hitObjectIndex = 0;
		double t = 0;
		if (intersect(ray, intersectionPoint, intersectionNormal, t, hitObjectIndex))
		{

			if (numBounces < 0)
			{
				return { 0, 0, 0 };
				//return { 10000000, 0, 0 };
			}

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
			Vector3 PL = lightPosition - intersectionPoint;
			double d2 = PL.norm2();
			PL = PL / sqrt(d2);

			Ray shadowRay = { intersectionPoint + EPSILON * intersectionNormal, PL };
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
				color = objects[hitObjectIndex].albedo * lightIntensity * dot(intersectionNormal, PL) / (4 * M_PI * d2 * M_PI);
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

			// Extended source
			if (hitObjectIndex == 0)
			{
				return objects[hitObjectIndex]->albedo * lightIntensity / (4 * sqr(M_PI * dynamic_cast<const Sphere*>(objects[hitObjectIndex])->radius));
			}

			Vector3 color = { 0, 0, 0 };
			// Indirect lighting
			Vector3 randomRayDirection = random_cos(intersectionNormal, intersectionPoint); // randomly sample ray us ing random cos
			Vector3 indirectLightingColor = objects[hitObjectIndex]->albedo * getColor({ intersectionPoint + EPSILON * intersectionNormal, randomRayDirection }, numBounces - 1);
			return color + indirectLightingColor;
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
	double cameraFOV = 60.0 * M_PI / 180.0; // FOV angle
	Vector3 camOrigin( 0.0, 0.0, 55.0 );
	Scene scene;
	scene.addObject(new Sphere({ { 0, 40, 0 },  {1.0, 1.0, 1.0}, 18.0,   false, false,  1.3 })); // Extended light !!
	TriangleMesh* mesh = new TriangleMesh({ {1.0, 1.0, 1.0}, false, false, 1.0 });
	mesh->readOBJ("resources/cat/cat.obj");
	mesh->scaleTranslate(0.6, { 0, -10, 0 });
	ComputeAABB(mesh->aabb, mesh->vertices);
	scene.addObject(mesh);
	//                Center,          albedo,          radius, isMirror, isTransparent, refractiveIndex
	//scene.addObject({ { 0, 0, 0 },     {0.6, 0.1, 0.8}, 10.0,  false, false, 1.0 }); // Mat sphere 1
	//scene.addObject(new Sphere({ { 12, 12, 0 },   {0.1, 0.1, 0.4}, 5.0,   true,  false, 1.0 })); // Mirror sphere
	//scene.addObject(new Sphere({ { -12, 12, 0 },  {0.3, 0.7, 0.3}, 5.0,   false, true,  2.0 })); // Transparent sphere
	//scene.addObject({ { 0, 12, 16 },   {0.1, 0.1, 0.4}, 2.0,   false, true,  1.3 }); // Transparent sphere 2
	//scene.addObject({ { 5, -3, 16 },  {0.1, 0.1, 0.4}, 3.0,    false, true,  1.3 }); // Transparent sphere 3
	//scene.addObject({ { -5, -3, 12 },   {0.1, 0.1, 0.4}, 3.0,  true,  false, 1.3 }); // Mirror sphere 2

	scene.addObject(new Sphere({ { 1000, 0, 0 },  {0.7, 0.3, 0.3}, 940.0, false, false, 1.0 })); // Right red wall
	scene.addObject(new Sphere({ { -1000, 0, 0 }, {0.2, 0.2, 0.7}, 940.0, false, false, 1.0 })); // Left blue wall
	scene.addObject(new Sphere({ { 0, 1000, 0 },  {0.3, 0.7, 0.3}, 940.0, false, false, 1.0 })); // Top green ceiling
	scene.addObject(new Sphere({ { 0, -950, 0 },  {0.3, 0.7, 0.3}, 940.0, false, false, 1.0 })); // Bottom green floor
	scene.addObject(new Sphere({ { 0, 0, 1000 },  {0.7, 0.7, 0.3}, 940.0, false, false, 1.0 })); // Back yellow wall
	scene.addObject(new Sphere({ { 0, 0, -1000 }, {0.2, 0.2, 0.7}, 940.0, false, false, 1.0 })); // Front blue wall
	const double gamma = 2.2;

#ifdef _OPENMP
	std::cout << "[INFO] OpenMP is active.\n";
#endif

	char* image = new char[W * H * 3];
	std::cout << "[INFO] Rendering started.\n";
#pragma omp parallel for
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
				color += scene.getColor({ locCamOrigin, localCamDirection}, GLOBAL_NUM_BOUNCES); // TO DO : The diaphragm, in theory, should not be gaussian, but circular with rejection, or hexagonal
			}
			setImageColor(image, x * W + y, gammaCorrect(color / RENDER_STEP_COUNT, gamma));
		}
	}
	timer.stop();
	stbi_write_png("outputImage/course5.png", W, H, 3, &image[0], 0);
	delete[] image;
	std::cout << "[INFO] Image generation done.\n";
	std::cout << "[INFO] Rendering time : " << timer.elapsedMilliseconds() << "ms.\n";
	return 0;
}