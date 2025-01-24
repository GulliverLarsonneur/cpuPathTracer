#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <iostream>
#include <chrono>
#include <ctime>
#include <cmath>
#include <limits>
#include <random>
#include <cmath>
#include <omp.h>

#define IMAGE_WIDTH 512
#define M_PI 3.1415926535897932384626
#define EPSILON 0.00001
#define AIR_INDEX 1.0
#define RENDER_STEP_COUNT 1024
#define AA_RADIUS 0.33

// TO DO : make multiple different random engines !
// TO DO : Use a custom hash
// TO DO : Implement Sobol and Halton sequences Quasi Monte-Carlo sampling



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

typedef struct Ray
{
	Vector3 origin;
	Vector3 direction;
} Ray;



static std::default_random_engine engine; // random s e ed = 10
static std::uniform_real_distribution<double> uniform( 0 , 1) ;

// generates a gaussian distribution with two uniform distributions as an input (x and y) and a standard deviation :
inline void boxMuller( double stdev , double &x , double &y )
{
	double r1 = uniform(engine) ;
	double r2 = uniform(engine) ;
	x = sqrt( -2 * log(r1)) * cos(2 * M_PI * r2) * stdev;
	y = sqrt( -2 * log(r1)) * sin(2 * M_PI * r2) * stdev;
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
				double alphaR = k0 + (1 - k0) * std::powf((1 - abs(dot(intersectionNormal, ray.direction))), 5); // alphaR = Reflection coefficient

				//int threadID = omp_get_thread_num(); // TO DO : fix this !
				if (alphaR > hash13(intersectionPoint))//uniform(engine))//uniform(engine[threadID]))
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

			Ray shadowRay = { intersectionPoint, PL };
			double shadow_t = 0;
			bool isLight = true;
			shadowRay.origin += EPSILON * intersectionNormal;
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

			if (isLight)
			{
				Vector3 color = objects[hitObjectIndex].albedo * lightIntensity * dot(intersectionNormal, PL) / (4 * M_PI * d2 * M_PI);
				return color;
			}
			else
			{
				return { 0, 0, 0 };
			}
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
	Timer timer;
	timer.start();
	const int W = IMAGE_WIDTH;
	const int H = IMAGE_WIDTH;
	double cameraFOV = 60.0 * M_PI / 180.0; // FOV angle
	Vector3 camOrigin( 0.0, 0.0, 55.0 );
	Scene scene;
	scene.addSphere({ { 0, 0, 0 },  {0.3, 0.7, 0.3}, 9.0,   false, true,  1.3 }); // Transparent sphere

	//                Center,          albedo,          radius, isMirror, isTransparent, refractiveIndex
	//scene.addSphere({ { 0, 0, 0 },     {0.6, 0.1, 0.8}, 10.0,  false, false, 1.0 }); // Mat sphere 1
	//scene.addSphere({ { 12, 12, 0 },   {0.1, 0.1, 0.4}, 9.0,   true,  false, 1.0 }); // Mirror sphere
	//scene.addSphere({ { -12, 12, 0 },  {0.3, 0.7, 0.3}, 9.0,   false, true,  1.1 }); // Transparent sphere
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

				color += scene.getColor({ camOrigin, u }, 10);
			}
			setImageColor(image, x * W + y, gammaCorrect(color / RENDER_STEP_COUNT, gamma));
		}
	}
	timer.stop();
	stbi_write_png("outputImage/course2.png", W, H, 3, &image[0], 0);
	delete[] image;
	std::cout << "[INFO] Image generation done.\n";
	std::cout << "[INFO] Rendering time : " << timer.elapsedMilliseconds() << "ms.\n";
	return 0;
}