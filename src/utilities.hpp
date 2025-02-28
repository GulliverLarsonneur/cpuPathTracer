#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <thread>
#define MAX_BVH_DEPTH 25 // TO DO : optimize this
#define BVH_LEAF_TRIANGLE_COUNT 4 // Comes from the cost of intersecting 4 triangles VS the cost of intersecting a AABB


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


class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};


class Geometry
{
public:
	Geometry(Vector3 _albedo,
	bool _isMirror,
	bool _isTransparent,
	float _refractiveIndex) :albedo(_albedo), isMirror(_isMirror), isTransparent(_isTransparent), refractiveIndex(_refractiveIndex){};
	virtual bool intersect(const Ray& ray, Vector3& intersectionPoint, Vector3& intersectionNormal, double& t, Vector3& col) const = 0;
public:
	Vector3 albedo;
	bool isMirror;
	bool isTransparent;
	float refractiveIndex;
};


typedef struct AABB
{
	bool intersect(const Ray& ray) const
	{
		double P1x = (minCoords[0] - ray.origin[0]) / ray.direction[0];
		double P2x = (maxCoords[0] - ray.origin[0]) / ray.direction[0];
		double xmin = std::min(P1x, P2x);
		double xmax = P1x + P2x - xmin;

		double P1y = (minCoords[1] - ray.origin[1]) / ray.direction[1];
		double P2y = (maxCoords[1] - ray.origin[1]) / ray.direction[1];
		double ymin = std::min(P1y, P2y);
		double ymax = P1y + P2y - ymin;

		double P1z = (minCoords[2] - ray.origin[2]) / ray.direction[2];
		double P2z = (maxCoords[2] - ray.origin[2]) / ray.direction[2];
		double zmin = std::min(P1z, P2z);
		double zmax = P1z + P2z - zmin;

		double minOfMax = std::min(xmax, std::min(ymax, zmax));
		double maxOfMin = std::max(xmin, std::max(ymin, zmin));
		
		if (minOfMax < 0) // There is a reason for this
			return false;
		return (minOfMax > maxOfMin);
	}
	AABB(): minCoords(0.0, 0.0, 0.0), maxCoords(0.0, 0.0, 0.0) {}
	AABB(Vector3 _minCoords, Vector3 _maxCoords): minCoords(_minCoords), maxCoords(_maxCoords) {}
	Vector3 minCoords;
	Vector3 maxCoords;
};


class BVH
{
public:
	BVH() : leftChild(nullptr), rightChild(nullptr), aabb() {};
	int startRange; // Triangle indices in the sorted array
	int endRange;   // Triangle indices in the sorted array
	BVH* leftChild;
	BVH* rightChild;
	AABB aabb;
};


class TriangleMesh : public Geometry
{
public:
  ~TriangleMesh() {}
	TriangleMesh(Vector3 _albedo,
		bool _isMirror,
		bool _isTransparent,
		float _refractiveIndex) : Geometry(_albedo, _isMirror, _isTransparent, _refractiveIndex), meshAABB(Vector3{0, 0, 0}, Vector3{0, 0, 0}) {};

	void scaleTranslate(double scale, const Vector3& translation)
	{
		for (int vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex)
		{
			vertices[vertexIndex] = scale * vertices[vertexIndex]; // TO do : define *= operator
			vertices[vertexIndex] += translation;
		}
	}
	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector3 vec = { 0, 0, 0 };

				Vector3 col = { 0, 0, 0 };
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector3 vec = { 0, 0, 0 };
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector3 vec = { 0, 0, 0 };
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}
	void loadTexture(const char* fileName)
	{
		int w;
		int h;
		int channels_in_file;
		unsigned char* image = stbi_load(fileName, &w, &h, &channels_in_file, 3);

		std::vector<float> catTexture;
		for (int i = 0; i < w * h * 3; i++)
		{
			catTexture.push_back(std::pow(double(image[i])/256.0, 2.2));
		}

		texW.push_back(w);
		texH.push_back(h);
		textures.push_back(catTexture);
	}

	void initBVH()
	{
		meshBVH = new BVH();
		buildBVH(meshBVH, 0, indices.size());
	}

	void buildBVH(BVH* bvh, int startRange, int endRange) // This is a quick sort algorithm
	{
		bvh->startRange = startRange;
		bvh->endRange = endRange;
		bvh->leftChild = nullptr;
		bvh->rightChild = nullptr;
		ComputeAABB(bvh->aabb, bvh->startRange, bvh->endRange);
		if (startRange - endRange < BVH_LEAF_TRIANGLE_COUNT) return;

		Vector3 diag = bvh->aabb.minCoords - bvh->aabb.maxCoords;
		int axis = 2;
		if (diag[0] >= diag[1] && diag[0] >= diag[2])
			axis = 0;
		if (diag[1] >= diag[0] && diag[1] >= diag[2])
			axis = 1;
		double milieuBoxAxis = (bvh->aabb.minCoords[axis] + bvh->aabb.maxCoords[axis]) / 2;
		int pivot = bvh->startRange;
		for (int triangleIndex = bvh->startRange; triangleIndex < bvh->endRange; triangleIndex++)
		{
			double barycentreTriangleAxis = (vertices[indices[triangleIndex].vtxi][axis] +
				vertices[indices[triangleIndex].vtxj][axis] + 
				vertices[indices[triangleIndex].vtxk][axis]) / 3;
			
			if (barycentreTriangleAxis < milieuBoxAxis)
			{
				if (bvh->leftChild)
				{
					std::swap(indices[triangleIndex], indices[pivot]);
					pivot++;
				}
			}
		}

		if (pivot - startRange == 0) return;
		if (endRange - pivot == 0) return; // We absolutely don't want all triangles to be on one side
		bvh->leftChild = new BVH();
		bvh->rightChild = new BVH();
		buildBVH(bvh->leftChild, bvh->startRange, pivot);
		buildBVH(bvh->rightChild, pivot, bvh->endRange);
	}

	bool intersect(const Ray& ray, Vector3& intersectionPoint, Vector3& intersectionNormal, double& best_t, Vector3& col) const
	{
		// optimization TO DO : add the best_t as a parameter to intersect, and compare the intersection point distance to best_t that we currently have.
		if (!meshAABB.intersect(ray))
		{
			return false;
		}

		bool hasIntersected = false;
		BVH* BVHstack[MAX_BVH_DEPTH];
		BVHstack[0] = meshBVH;
		int sizeBVHStack = 1;
		while (sizeBVHStack > 0)
		{
			BVH* currentBVH = BVHstack[sizeBVHStack - 1];
			sizeBVHStack--;
			if (currentBVH->leftChild) // A BVH with a right child ALWAYS has a left child
			{
				if (currentBVH->leftChild->aabb.intersect(ray))
				{
					sizeBVHStack++; // Test ??
					BVHstack[sizeBVHStack] = currentBVH->leftChild;
				}
				if (currentBVH->rightChild->aabb.intersect(ray))
				{
					sizeBVHStack++;
					BVHstack[sizeBVHStack] = currentBVH->rightChild;
				}
			}
			else
			{
				best_t = std::numeric_limits<double>::max();
				for (int triangleIndex = currentBVH->startRange; triangleIndex < currentBVH->endRange; ++triangleIndex)
				{
					const Vector3& A = vertices[indices[triangleIndex].vtxi];
					const Vector3& B = vertices[indices[triangleIndex].vtxj];
					const Vector3& C = vertices[indices[triangleIndex].vtxk];
					const Vector3 e1 = B - A;
					const Vector3 e2 = C - A;
					Vector3 N = cross(e1, e2);

					double invDet = 1.0 / dot(ray.direction, N);
					const Vector3 AO = ray.origin - A;
					const Vector3 AOcrossU = cross(AO, ray.direction);

					double beta = - dot(e2, AOcrossU) * invDet;
					if ((beta < 0) || (beta > 1)) continue;
					double gamma = dot(e1, AOcrossU) * invDet;
					if ((gamma < 0) || (gamma > 1)) continue;
					double alpha = 1 - beta - gamma;
					if (alpha < 0) continue; // No need to check alpha > 1
					double t = - dot(AO, N) * invDet;
					if (t < 0) continue; // Triangle behind the camera
					hasIntersected = true;
					if (t > best_t) continue;
					best_t = t;

					intersectionPoint = ray.origin + t * ray.direction;

					const Vector3& NA = normals[indices[triangleIndex].ni];
					const Vector3& NB = normals[indices[triangleIndex].nj];
					const Vector3& NC = normals[indices[triangleIndex].nk];
					intersectionNormal = normalize(NA * alpha + NB * beta + NC * gamma);
					if (textures.size() != 0)
					{
						// This should be a Vector2
						Vector3 uv = alpha * uvs[indices[triangleIndex].uvi] + beta * uvs[indices[triangleIndex].uvj] + gamma * uvs[indices[triangleIndex].uvk];// group sert à donner l'ID de la texture dans le cas où il y a plusieurs textures
						int w = texW[indices[triangleIndex].group];
						int h = texH[indices[triangleIndex].group];
						//std::cout << "(u0,v0) = (" << uv[0] << "," << uv[1] << ")" << std::endl;
						uv[0] = fmod(uv[0] + 10000, 1.) * w; // +10000 => Trick to prevent the error with fmod of negative values
						uv[1] = (1.0 - fmod(uv[1] + 10000, 1.)) * h;
						int uvx = uv[0];
						int uvy = uv[1];

						//std::cout << "(h, w) = (" << h << ", " << w << ")  (u1,v1) = (" << uvx << "," << uvy << "), ID = " << 3 * (uvx + uvy * w) << " group = " << indices[triangleIndex].group << " texture size : " << textures[indices[triangleIndex].group].size() << std::endl;

						col = Vector3(textures[indices[triangleIndex].group][3 * (uvx + uvy * w)], textures[indices[triangleIndex].group][3 * (uvx + uvy * w) + 1], textures[indices[triangleIndex].group][3 * (uvx + uvy * w) + 2]); // TO DO : finish this !!
						//std::cout << "col = " << col << "\n" << std::endl;
					}
				}
			}
		}
		return hasIntersected;
	}

	void ComputeAABB(AABB& aabb, int startRange = 0, int endRange = -1)
	{
		if (endRange == -1)
		{
			endRange = vertices.size();
		}
		Vector3 minCoords = std::numeric_limits<double>::max() * Vector3{ 1.0, 1.0, 1.0 };
		Vector3 maxCoords = std::numeric_limits<double>::max() * Vector3{ -1.0, -1.0, -1.0 };

		for (int triangleIndex = startRange; triangleIndex < endRange; ++triangleIndex)
		{
			for (int j = 0; j < 3; ++j)
			{
				minCoords[j] = std::min(vertices[indices[triangleIndex].vtxi][j], minCoords[j]);
				minCoords[j] = std::min(vertices[indices[triangleIndex].vtxj][j], minCoords[j]);
				minCoords[j] = std::min(vertices[indices[triangleIndex].vtxk][j], minCoords[j]);
				maxCoords[j] = std::max(vertices[indices[triangleIndex].vtxi][j], maxCoords[j]);
				maxCoords[j] = std::max(vertices[indices[triangleIndex].vtxj][j], maxCoords[j]);
				maxCoords[j] = std::max(vertices[indices[triangleIndex].vtxk][j], maxCoords[j]);
			}
		}
		aabb.minCoords = minCoords;
		aabb.maxCoords = maxCoords;
		//std::cout << minCoords << " " << maxCoords << "\n";
	}
public:
	std::vector<int> texW;
	std::vector<int> texH;
	std::vector<std::vector<float>> textures;
	std::vector<TriangleIndices> indices;
	std::vector<Vector3> vertices;
	std::vector<Vector3> normals;
	std::vector<Vector3> uvs;
	std::vector<Vector3> vertexcolors;
	AABB meshAABB;
	BVH* meshBVH;
};

typedef struct Sphere : public Geometry
{
public:
	Sphere(Vector3 _center, Vector3 _albedo,
		float _radius,
		bool _isMirror,
		bool _isTransparent,
		float _refractiveIndex) : center(_center), radius(_radius), Geometry(_albedo, _isMirror, _isTransparent, _refractiveIndex){}
	bool intersect(const Ray& ray, Vector3& intersectionPoint, Vector3& intersectionNormal, double& t, Vector3& col) const
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
		//col = albedo;
		double frequency = 0.2;
		//double lightness = abs(fmod(frequency * (intersectionPoint[0] + intersectionPoint[1] + intersectionPoint[2]), 1.0));
		//double lightness = abs(fmod(frequency * (intersectionPoint[0]), 1.0) * fmod(intersectionPoint[1], 1.0) * fmod(intersectionPoint[2], 1.0));
		
		int a0 = (int)(frequency * intersectionPoint[0] + 100000000);
		int a1 = (int)(frequency * intersectionPoint[1] + 100000000);
		int a2 = (int)(frequency * intersectionPoint[2] + 100000000);
		double lightness = (double)((a0 + a1 + a2) % 2);


		//col = Vector3(lightness + hash13({(double)a0 - 100000000, (double)a1 - 100000000, (double)a2 - 100000000}), lightness + hash13({(double)a0 - 100000000, (double)a2 - 100000000, (double)a1 - 100000000}), lightness + hash13({(double)a2 - 100000000, (double)a1 - 100000000, (double)a0 - 100000000}));
		double colorFactor = 0.6;
		col = Vector3(lightness + colorFactor * hash13({(double)a0 - 100000000, (double)a1 - 100000000, (double)a2 - 100000000}), lightness + colorFactor * hash13({(double)a0 - 100000000, (double)a2 - 100000000, (double)a1 - 100000000}), lightness + colorFactor * hash13({(double)a2 - 100000000, (double)a1 - 100000000, (double)a0 - 100000000}));

		return true;
	}
	Vector3 center;
	double radius;
} Sphere;


