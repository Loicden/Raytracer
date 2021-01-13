#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		coords[0] = x;
		coords[1] = y;
		coords[2] = z;
	};
	double operator[](int i) const { return coords[i]; };
	double& operator[](int i) { return coords[i]; };
	double sqrNorm() const {
		return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
	}
	Vector get_normalized() {
		double n = sqrt(sqrNorm());
		return Vector(coords[0] / n, coords[1] / n, coords[2] / n);
	}
private:
	double coords[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(double a, const Vector& b) {
	return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector& a, double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(double a, const Vector& b) {
	return Vector(a / b[0], a / b[1], a / b[2]);
}
Vector operator/(const Vector& a, double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
double sqr(double x) {
	return x * x;
}

class Ray {
public:
	Ray(const Vector& C, const Vector& u) : C(C), u(u) {};
	Vector C, u;
};

class Sphere {
public:
	Sphere(const Vector& O, double R, const Vector &albedo, bool isMirror = false, bool isTransparent = false) : O(O), R(R), albedo(albedo), isMirror(isMirror), isTransparent(isTransparent) {}
	bool intersect(const Ray& r, Vector& P, Vector& N, double& t) {
		// solves a*t^2 + b*t + c = 0
		double a = 1;
		double b = 2 * dot(r.u, r.C - O);
		double c = (r.C - O).sqrNorm() - R * R;
		double delta = b * b - 4 * a*c;
		if (delta < 0) return false;

		double sqDelta = sqrt(delta);
		double t2 = (-b + sqDelta) / (2 * a);
		if (t2 < 0) return false;
		double t1 = (-b - sqDelta) / (2 * a);
		if (t1 > 0)
			t = t1;
		else
			t = t2;

		P = r.C + t * r.u; // la direction P est définié par l'origine du rayon et t* la direction du rayon
		N = (P - O).get_normalized();

		return true;
	};
	Vector O;
	double R;
	Vector albedo;
	bool isMirror, isTransparent;
};

class Scene {
public:
	Scene() {};
	bool intersect(const Ray& r, Vector& P, Vector& N, Vector &albedo, bool &mirror, bool &transp, double &t) {

		t = 1E10;
		bool has_inter = false;
		for (int i = 0; i < objects.size(); i++) {
			Vector localP, localN;
			double localt;
			bool local_has_inter = objects[i].intersect(r, localP, localN, localt);
			if (local_has_inter && localt < t) {
				has_inter = true;
				t = localt;
				albedo = objects[i].albedo;
				P = localP;
				N = localN;
				mirror = objects[i].isMirror;
				transp = objects[i].isTransparent;
			}
		}
		return has_inter;
	}

	Vector getColor(const Ray& r, int rebond) {
		if (rebond > 5) {
			return Vector(0., 0., 0.);
		}
		Vector P, N, albedo;
		double t;
		double M_PI = 3.14159;
		bool transp, mirror;
		bool inter = intersect(r, P, N, albedo, mirror, transp, t);
		Vector color(0, 0, 0);

		if (inter) {

			if (mirror) {
				Vector reflectedDir = r.u - 2 * dot(r.u, N) * N;
				Ray reflectedRay(P + 0.001*N, reflectedDir);
				return getColor(reflectedRay, rebond + 1);
			}
			else {
				if (transp) {
					Vector reflectedDir = r.u - 2 * dot(r.u, N) * N;
					Ray reflectedRay(P + 0.001 * N, reflectedDir);
					double n1 = 1, n2 = 1.4;
					Vector N2 = N; 
					if (dot(r.u, N) > 0) { //on sort de la sphère
						std::swap(n1, n2);
						N2 = -N;
					}
					Vector Tt = n1 / n2 * (r.u - dot(r.u, N2) * N2);
					double rad = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.u, N2)));
					if (rad < 0) {
						Vector reflectedDir = r.u - 2 * dot(r.u, N) * N;
						Ray reflectedRay(P + 0.0001 * N, reflectedDir);
						return getColor(reflectedRay, rebond + 1);
					}
					Vector Tn = -sqrt(rad)*N2;

					Vector refractedDir = Tt + Tn;

					return getColor(Ray(P - 0.0001*N2, refractedDir), rebond + 1);
				}
				else {
					Vector PL = Lum - P;
					double d = sqrt(PL.sqrNorm());
					double ddot = dot(N, PL / d);
					Vector shadowP, shadowN, shadowAlbedo;
					double shadowt;
					bool shadowMirror, shadowTransp;
					Ray shadowRay(P + 0.001 * N, PL / d);
					bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransp, shadowt);
					if (shadowInter && shadowt < d) {
						color = Vector(0., 0., 0.);
					}
					else {
						color = I / (4 * M_PI * d * d) * albedo / M_PI * std::max(0., dot(N, PL / d));
					}
				}
			}
		}
		return color;
	}
	std::vector<Sphere> objects;
	Vector Lum;
	double I;

};

int main() {
	int W = 1024;
	int H = 1024;
	
	Sphere s1(Vector(0, 0, 0), 10, Vector(1, 1, 1),false, true);
	//Sphere s2(Vector(5, 5, 0), 7, Vector(1, 1, 1), false, false);
	//Sphere s3(Vector(-5, 5, 0), 7, Vector(1, 1, 1), false, false);

	Sphere s4(Vector(0, 0, 60), 4, Vector(1, 0, 1), false);

	Sphere ssol(Vector(0, -1000, 0), 990, Vector(1, 1, 1)); //sol
	Sphere splafond(Vector(0, 1000, 0), 940, Vector(1, 1, 1)); //plafond
	Sphere sgauche(Vector(-1000, 0, 0), 940, Vector(1, 0, 0)); //mur gauche
	Sphere sdroite(Vector(1000, 0, 0), 940, Vector(0, 1, 0)); //mur droite
	Sphere sfond(Vector(0, 0, -1000), 940, Vector(0, 0, 1)); //mur du fond
	Sphere sarriere(Vector(0, 0, 1000), 940, Vector(1, 1, 0)); //mur arrière

	Scene scene;

	scene.objects.push_back(s1);
	//scene.objects.push_back(s2);
	//scene.objects.push_back(s3);
	scene.objects.push_back(s4);
	scene.objects.push_back(ssol);
	scene.objects.push_back(splafond);
	scene.objects.push_back(sgauche);
	scene.objects.push_back(sdroite);
	scene.objects.push_back(sfond);
	scene.objects.push_back(sarriere);

	double M_PI = 3.14159;
	double fov = 60 * M_PI / 180;
	Vector C(0, 0, 55);
	scene.I = 5E9;
	scene.Lum = Vector (-10, 20, 40);

	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			Vector u(j - W/2, i - H/2, -W/(2*tan(fov/2)));
			u = u.get_normalized();
			Ray r(C, u);
			Vector color = scene.getColor(r, 0);
			
		
			image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));
			image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));
			image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));

		}
	}
	stbi_write_png("titouan.png", W, H, 3, &image[0], 0);

	return 0;
}