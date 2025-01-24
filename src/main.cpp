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

#define M_PI 3.1415
#define EPSILON 0.00001

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


static inline double sqr(double x) 
{ 
	return x * x; 
}


typedef struct Vector3 
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
	Vector3 color;
	Vector3 albedo;
	double radius;
} Sphere;


#define IMAGE_WIDTH 512

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

	Vector3 getColor(Ray ray)
	{
		Vector3 intersectionPoint = { 0.0, 0.0 ,0.0 };
		Vector3 intersectionNormal = { 0.0, 0.0, 0.0 };
		int hitObjectIndex = 0;
		double t = 0;
		if (intersect(ray, intersectionPoint, intersectionNormal, t, hitObjectIndex))
		{
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

int main() 
{
	Timer timer;
	timer.start();
	const int W = IMAGE_WIDTH;
	const int H = IMAGE_WIDTH;
	double alpha = 60.0 * M_PI / 180.0; // FOV angle
	Vector3 camOrigin( 0.0, 0.0, 55.0 );
	Scene scene;
	scene.addSphere({ { 0, 0, 0 }, {0.5, 0.2, 0.9}, {0.6, 0.1, 0.8}, 10.0 }); // Center, color, albedo, radius
	scene.addSphere({ { 1000, 0, 0 }, {0.5, 0.2, 0.9}, {0.7, 0.3, 0.3}, 940.0 });
	scene.addSphere({ { -1000, 0, 0 }, {0.5, 0.2, 0.9}, {0.3, 0.3, 0.7}, 940.0 });
	scene.addSphere({ { 0, 1000, 0 }, {0.5, 0.2, 0.9}, {0.3, 0.7, 0.3}, 940.0 });
	scene.addSphere({ { 0, -950, 0 }, {0.5, 0.2, 0.9}, {0.3, 0.7, 0.3}, 940.0 });
	scene.addSphere({ { 0, 0, 1000 }, {0.5, 0.2, 0.9}, {0.3, 0.3, 0.7}, 940.0 });
	scene.addSphere({ { 0, 0, -1000 }, {0.5, 0.2, 0.9}, {0.3, 0.3, 0.7}, 940.0 });
	const double gamma = 2.2;
#ifdef _OPENMP
	std::cout << "OpenMP is active.\n";
#endif

	char* image = new char[W * H * 3];

#pragma omp parallel for
	for (int x = 0; x < H; x++) 
	{
		for (int y = 0; y < W; y++) 
		{
			Vector3 u = normalize(Vector3(y - W / 2, -(x - H / 2), - W / ( 2 * tan(alpha / 2) )));
			Vector3 intersectionPoint = {0, 0, 0}, intersectionNormal = {0, 0, 0};
			Vector3 color = scene.getColor({ camOrigin, u });
			setImageColor(image, x * W + y, gammaCorrect(color, gamma));
		}
	}
	stbi_write_png("outputImage/radial_gradient.png", W, H, 3, &image[0], 0);
	
	timer.stop();
	delete[] image;
	std::cout << "Milliseconds: " << timer.elapsedMilliseconds() << "\n";
	return 0;
}