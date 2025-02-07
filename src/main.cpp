#define _CRT_SECURE_NO_WARNINGS 1


// TO DO : put all the parameters in a static struct
#define IMAGE_WIDTH 512
#define M_PI 3.1415926535897932384626
#define EPSILON 0.00001
#define AIR_INDEX 1.0
#define RENDER_STEP_COUNT 2048
#define AA_RADIUS 0.33
#define GLOBAL_NUM_BOUNCES 5
#define DEPTH_OF_FIELD_AMPLITUDE 10.0
#define DEPTH_OF_FIELD_DISTANCE 55.0
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
#include <thread>

// Each thread has its own random engine instance
thread_local std::default_random_engine engine{ std::random_device{}() };
static std::uniform_real_distribution<double> uniform(0, 1);

// generates a gaussian distribution with two uniform distributions as an input (x and y) and a standard deviation :
inline void boxMuller( double stdev , double &x , double &y )
{
	double r1 = uniform(engine) ;
	double r2 = uniform(engine) ;
	x = sqrt( -2 * log(r1)) * cos(2 * M_PI * r2) * stdev;
	y = sqrt( -2 * log(r1)) * sin(2 * M_PI * r2) * stdev;
}

static inline double sqr(double x) 
{ 
	return x * x; 
}

typedef struct Vector3 // TO DO : make this with simple x, y, z components
{
	Vector3(double x, double y, double z) 
	{
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}
	double& operator[](int i) { return coord[i]; }
	double operator[](int i) const { return coord[i]; }

	Vector3& operator+=(const Vector3& v) 
	{
		coord[0] += v[0];
		coord[1] += v[1];
		coord[2] += v[2];
		return *this;
	}

	double norm2() const 
	{
		return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
	}

	Vector3 getNormalized()
	{
		const double norm = sqrt(this->norm2());
		return { coord[0] / norm, coord[1] / norm, coord[2] / norm };
	}
	double coord[3];
};

// ostream operator
std::ostream& operator<<(std::ostream& os, const Vector3& v) {
	os << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
	return os;
}


Vector3 normalize(const Vector3& v)
{
	const double norm = sqrt(v.norm2());
	return { v[0] / norm, v[1] / norm, v[2] / norm };
}


Vector3 operator+(const Vector3& a, const Vector3& b) 
{
	return Vector3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector3 operator-(const Vector3& a, const Vector3& b) 
{
	return Vector3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector3 operator*(const Vector3& a, double b) 
{
	return Vector3(a[0]*b, a[1]*b, a[2]*b);
}
Vector3 operator*(double a, const Vector3& b) 
{
	return Vector3(a*b[0], a*b[1], a*b[2]);
}
Vector3 operator/(const Vector3& b, double a) 
{
	return Vector3(b[0]/a, b[1]/a, b[2]/a);
}

double dot(const Vector3& a, const Vector3& b) 
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector3 cross(const Vector3& a, const Vector3& b) 
{
	return { a[1] * b [2] - a[2] * b[1], a[2] * b[0] - b[2] * a[0], a[0] * b[1] - a[1] * b[0]};
}

Vector3 operator*(const Vector3& a, const Vector3& b) 
{
	return Vector3(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}

union FloatToInt {
	float f;
	int i;
};

inline int floatToInt(float value) {
	FloatToInt temp;
	temp.f = value;
	return temp.i;
}

inline float fract(float value) {
	return value - std::floor(value);
}


float hash11(float p) // From https://www.shadertoy.com/view/4djSRW
{
	p = fract(p * .1031);
	p *= p + 33.33;
	p *= p + p;
	return fract(p);
}

float hash13(Vector3 p3) // From https://www.shadertoy.com/view/4djSRW
{
	p3 = { fract((p3 * 0.1031)[0]), fract((p3 * 0.1031)[1]), fract((p3 * 0.1031)[2]) };
	p3 += Vector3({p3[0] * (p3[2] + 31.32), p3[1] * (p3[1] + 31.32), p3[2] * (p3[0] + 31.32)});
	return fract((p3[0] + p3[1]) * p3[2]);
}



Vector3 random_cos(const Vector3& intersectionNormal, const Vector3& intersectionPoint)
{
	double u1 = uniform(engine);
	double u2 = uniform(engine);
	//double u1 = hash13(165.265632 * hash13(2406.56216 * intersectionPoint) * intersectionNormal);
	//double u2 = hash13(144.2603585 * hash13(206.215616 * intersectionPoint) *intersectionNormal);
	
	//std::cout << u1 << "  " << u2 << "\n";
	/*
	if (u1 == 0)
	{
		std::cout << intersectionNormal << " " << intersectionPoint << " " << hash13(12635.265632 * intersectionNormal) << " " << hash13(245506.56216 * intersectionPoint) << "\n";
	}
	if (u2 == 0)
	{
		std::cout << intersectionNormal << " " << intersectionPoint << " " << hash13(12635.265632 * intersectionNormal) << " " << hash13(245506.56216 * intersectionPoint) << "\n";
	}
	*/
	//std::cout << hash13(12565.265632 * intersectionNormal) << "\n";
	double s = sqrt(1 - u2);
	double x = cos(2.0 * M_PI * u1) * s;
	double y = sin(2.0 * M_PI * u1) * s;
	double z = sqrt(u2);


	Vector3 tangentVector = { 0, 0, 0 };
	// L'idée est de garantir que <N|T> = 0 dans chacun des cas
	if (std::abs(intersectionNormal[0]) < std::abs(intersectionNormal[1]) && std::abs(intersectionNormal[0]) < std::abs(intersectionNormal[2]))
	{
		tangentVector = { 0, intersectionNormal[2], - intersectionNormal[1] };
	}
	else if (std::abs(intersectionNormal[1]) < std::abs(intersectionNormal[0]) && std::abs(intersectionNormal[1]) < std::abs(intersectionNormal[2]))
	{
		tangentVector = { intersectionNormal[2], 0, - intersectionNormal[0] };
	}
	else if (std::abs(intersectionNormal[2]) < std::abs(intersectionNormal[0]) && std::abs(intersectionNormal[2]) < std::abs(intersectionNormal[1]))
	{
		tangentVector = { intersectionNormal[1], - intersectionNormal[0], 0 };
	}
	tangentVector = normalize(tangentVector);
	
	Vector3 tangentVector2 = cross(tangentVector, intersectionNormal);
	return x * tangentVector + y * tangentVector2 + z * intersectionNormal;
}

typedef struct Ray
{
	Vector3 origin;
	Vector3 direction;
} Ray;



class Timer
{
public:
	void start()
	{
		m_StartTime = std::chrono::system_clock::now();
		m_bRunning = true;
	}

	void stop()
	{
		m_EndTime = std::chrono::system_clock::now();
		m_bRunning = false;
	}

	double elapsedMilliseconds()
	{
		std::chrono::time_point<std::chrono::system_clock> endTime;

		if(m_bRunning)
		{
			endTime = std::chrono::system_clock::now();
		}
		else
		{
			endTime = m_EndTime;
		}

		return std::chrono::duration_cast<std::chrono::milliseconds>(endTime - m_StartTime).count();
	}

	double elapsedSeconds()
	{
		return elapsedMilliseconds() / 1000.0;
	}

private:
	std::chrono::time_point<std::chrono::system_clock> m_StartTime;
	std::chrono::time_point<std::chrono::system_clock> m_EndTime;
	bool                                               m_bRunning = false;
};


typedef struct Sphere
{
	bool intersect(const Ray& ray, Vector3& intersectionPoint, Vector3& intersectionNormal, double& t)
	{
		// Resolve aX� + bX + c = 0 ; a = 1, b = 2 <ray.direction | ray.origin - center> , c = || origin - center || - radius * 2
		double b = 2 * dot(ray.direction, ray.origin - center);
		double delta = sqr(b) - 4 * ((ray.origin - center).norm2() - radius * radius);
		if (delta < 0)
			return false;

		double sqrtdelta = sqrt(delta);
		double t1 = (-b - sqrtdelta) / 2;
		double t2 = (-b + sqrtdelta) / 2;

		if (t2 < 0)
			return false;
		if (t1 >= 0)
			t = t1;
		else
			t = t2;

		//P = r.O + t * C
		intersectionPoint = ray.origin + t * ray.direction;
		intersectionNormal = normalize(intersectionPoint - center);
		return true;
	}
	Vector3 center;
	Vector3 albedo;
	double radius;
	bool isMirror;
	bool isTransparent;
	float refractiveIndex;
} Sphere;



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
			if (objects[objectIndex].intersect(ray, localIntersectionPoint, localIntersectionNormal, local_t))
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

			if (objects[hitObjectIndex].isTransparent || objects[hitObjectIndex].isMirror)
			{
				if (objects[hitObjectIndex].isMirror)
				{
					Ray bounceRay = { intersectionPoint + EPSILON * intersectionNormal, ray.direction - 2.0 * dot(ray.direction, intersectionNormal) * intersectionNormal };
					return getColor(bounceRay, numBounces - 1);
				}

				double indexRatio = 0;
				double n1 = AIR_INDEX;
				double n2 = objects[hitObjectIndex].refractiveIndex;
				Ray fresnelRay = { { 0, 0, 0 }, { 0, 0, 0 } };
				if (dot(ray.direction, intersectionNormal) < 0) // We enter inside the object
				{
					indexRatio = (AIR_INDEX / objects[hitObjectIndex].refractiveIndex);
					Vector3 Tt = indexRatio * (ray.direction - dot(ray.direction, intersectionNormal) * intersectionNormal);
					Vector3 Tn = -sqrt(1 - sqr(indexRatio) * (1 - sqr(dot(ray.direction, intersectionNormal)))) * intersectionNormal;
					fresnelRay = { intersectionPoint - EPSILON * intersectionNormal, Tn + Tt };
				}
				else // else, we exit the object
				{
					std::swap(n1, n2);
					indexRatio = (objects[hitObjectIndex].refractiveIndex / AIR_INDEX);
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

			// Indirect lighting
			Vector3 randomRayDirection = random_cos(intersectionNormal, intersectionPoint); // randomly sample ray us ing random cos
			Vector3 indirectLightingColor = objects[hitObjectIndex].albedo * getColor({ intersectionPoint + EPSILON * intersectionNormal, randomRayDirection }, numBounces - 1);
			return color + indirectLightingColor;
		}
		else
		{
			return { 0, 0, 0 };
		}
	}

	void addSphere(const Sphere& s)
	{
		objects.push_back(s);
	}

	std::vector<Sphere> objects;
	Vector3 lightPosition = { 10.0, 20.0, 40.0}; // x = -10 dans le pdf de cours
	double lightIntensity = 1000 * 2e07;
};

// TO DO : implement custom hash function

int main() 
{
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
	scene.addSphere({ { 0, 0, 0 },  {0.3, 0.7, 0.3}, 9.0,   false, false,  1.3 }); // Transparent sphere

	//                Center,          albedo,          radius, isMirror, isTransparent, refractiveIndex
	//scene.addSphere({ { 0, 0, 0 },     {0.6, 0.1, 0.8}, 10.0,  false, false, 1.0 }); // Mat sphere 1
	scene.addSphere({ { 12, 12, 0 },   {0.1, 0.1, 0.4}, 5.0,   true,  false, 1.0 }); // Mirror sphere
	scene.addSphere({ { -12, 12, 0 },  {0.3, 0.7, 0.3}, 5.0,   false, true,  2.0 }); // Transparent sphere
	//scene.addSphere({ { 0, 12, 16 },   {0.1, 0.1, 0.4}, 2.0,   false, true,  1.3 }); // Transparent sphere 2
	//scene.addSphere({ { 5, -3, 16 },  {0.1, 0.1, 0.4}, 3.0,    false, true,  1.3 }); // Transparent sphere 3
	//scene.addSphere({ { -5, -3, 12 },   {0.1, 0.1, 0.4}, 3.0,  true,  false, 1.3 }); // Mirror sphere 2

	scene.addSphere({ { 1000, 0, 0 },  {0.7, 0.3, 0.3}, 940.0, false, false, 1.0 }); // Right red wall
	scene.addSphere({ { -1000, 0, 0 }, {0.2, 0.2, 0.7}, 940.0, false, false, 1.0 }); // Left blue wall
	scene.addSphere({ { 0, 1000, 0 },  {0.3, 0.7, 0.3}, 940.0, false, false, 1.0 }); // Top green ceiling
	scene.addSphere({ { 0, -950, 0 },  {0.3, 0.7, 0.3}, 940.0, false, false, 1.0 }); // Bottom green floor
	scene.addSphere({ { 0, 0, 1000 },  {0.7, 0.7, 0.3}, 940.0, false, false, 1.0 }); // Back yellow wall
	scene.addSphere({ { 0, 0, -1000 }, {0.2, 0.2, 0.7}, 940.0, false, false, 1.0 }); // Front blue wall
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
				Vector3 camDestination = camOrigin + u * DEPTH_OF_FIELD_DISTANCE / u[2];
				Vector3 localCamDirection = normalize(camDestination - locCamOrigin);
				color += scene.getColor({ locCamOrigin, localCamDirection}, GLOBAL_NUM_BOUNCES); // The diaphragm, in theory, should not be gaussian, but circular with rejection, or hexagonal
			}
			setImageColor(image, x * W + y, gammaCorrect(color / RENDER_STEP_COUNT, gamma));
		}
	}
	timer.stop();
	stbi_write_png("outputImage/course3.png", W, H, 3, &image[0], 0);
	delete[] image;
	std::cout << "[INFO] Image generation done.\n";
	std::cout << "[INFO] Rendering time : " << timer.elapsedMilliseconds() << "ms.\n";
	return 0;
}