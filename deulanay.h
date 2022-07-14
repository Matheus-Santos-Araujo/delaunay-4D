#pragma once

#include <iostream>
#include <vector>
#include<GL/freeglut.h>
#include "vectorops.h"
#include "convexhull.h"

double angle(const Point& p1, const Point& p2) {

	double norm1 = sqrtf(p1.x * p1.x + p1.y * p1.y + p1.z * p1.z);

	double norm2 = sqrtf(p2.x * p2.x + p2.y * p2.y + p2.z * p2.z);

	double cos = dot(p1, p2) / (norm1 * norm2);

	if (cos < -1.0) {
		cos = -1.0;
	}
	else if (cos > 1.0) {
		cos = 1.0;
	}

	return acos(cos);
}

struct Tetra {
    union {
        struct { Point* p1, * p2, * p3, * p4; };
        Point* p[4];
    };
    vector<Face> faces;

    Point center() {
        return Point{ (p1->x + p2->x + p3->x + p4->x) / 4.0,  (p1->y + p2->y + p3->y + +p4->y) / 4.0, (p1->z + p2->z + p3->z + +p4->z) / 4.0 };
    }

    Tetra(Point *v1, Point *v2, Point *v3, Point *v4) {
        p1 = v1;
        p2 = v2;
        p3 = v3;
        p4 = v4;

		//face 1
		faces.push_back(Face{ v1, v2, v3 }); 

		//face 2
		faces.push_back(Face{ v2, v4, v3 }); 

		//face 3
		faces.push_back(Face{ v1, v4, v2 }); 

		//face 4
		faces.push_back(Face{ v1, v3, v4 }); 
    }
};

struct PolygonT {
	vector<Tetra> tetras;
};

	std::vector< std::vector<Point*> > vectorObj;
	std::vector<Point*> pointvector;
	std::vector<Face*> global_triangleList;
	std::vector<Face*> triangleList;
	std::vector<Face*> ftriangleList;
	std::vector<Tetra*> tetraList;

	bool addTriangle(Face* t) {
		for (std::vector<Face*>::iterator it = triangleList.begin(); it != triangleList.end(); it++) {
			if (t == (*it)) {
				t = (Face*)(*it);
				ftriangleList.erase(remove(ftriangleList.begin(), ftriangleList.end(), (*it)), ftriangleList.end());
				return false;
			}
		}

		global_triangleList.push_back(t);
		triangleList.push_back(t);
		ftriangleList.push_back(t);

		return true;
	}

	bool interceptTriangle(Face* t) {
		for (std::vector<Face*>::iterator it = triangleList.begin(); it != triangleList.end(); it++) {
			if (t == (*it) ||
				((t->p1 == (*it)->p1 || t->p2 == (*it)->p1 || t->p3 == (*it)->p1) && (t->p1 == (*it)->p2 || t->p2 == (*it)->p2 || t->p3 == (*it)->p2)) ||
				((t->p1 == (*it)->p1 || t->p2 == (*it)->p1 || t->p3 == (*it)->p1) && (t->p1 == (*it)->p3 || t->p2 == (*it)->p3 || t->p3 == (*it)->p3)) ||
				((t->p1 == (*it)->p2 || t->p2 == (*it)->p2 || t->p3 == (*it)->p2) && (t->p1 == (*it)->p3 || t->p2 == (*it)->p3 || t->p3 == (*it)->p3))) {
				continue;
			}

			/*if (intercept(t, *it)) {
				cout << "intercept \n";
				return true;
			}*/

		}

		return false;
	}

	void Delaunay3D(vector<Hull>& hulls) {
		for (Hull& h : hulls) {
			std::vector<PolygonT> polygon;
			PolygonT t;
			pointvector = h.vertices;
			
			// First Triangle - Convex Hull Face ----------------------------
			double minY = 1000000000000;
			Point* p0 = pointvector[0];
			for (int i = 0; i < pointvector.size(); i++) {

				if (pointvector[i]->y < minY) {
					minY = pointvector[i]->y;
					p0 = pointvector[i];
				}
			}

			Point p0_1 = Point{ p0->x + 1.0, p0->y, p0->z };
			Point p0_2 = Point{ p0->x, p0->y, p0->z - 1.0 };

			Point n1 = cross(Point{ p0->x - p0_1.x, p0->y - p0_1.y, p0->z - p0_1.z }, Point{ p0->x - p0_2.x, p0->y - p0_2.y, p0->z - p0_2.z });

			n1.norm();

			double ang = 360.0;

			for (std::vector<Point*>::iterator it = pointvector.begin(); it != pointvector.end(); it++) {
				if ((*it)->x == p0->x && (*it)->y == p0->y && (*it)->z == p0->z)
					continue;

				Point n2 = cross(Point{ p0->x - (*it)->x, p0->y - (*it)->y, p0->z - (*it)->z }, Point{ p0_2.x - (*it)->x, p0_2.y - (*it)->y, p0_2.z - (*it)->z });

				n2.norm();

				double ang_temp = angle(n2, n1);

				if (ang >= ang_temp) {
					ang = ang_temp;
					p0_1 = Point{ (*it)->x, (*it)->y, (*it)->z };
				}
			}

			n1 = cross(Point{ p0_1.x - p0->x, p0_1.y - p0->y, p0_1.z - p0->z }, Point{ p0_2.x - p0->x, p0_2.y - p0->y, p0_2.z - p0->z });

			n1.norm();

			ang = 360.0;

			for (std::vector<Point*>::iterator it = pointvector.begin(); it != pointvector.end(); it++) {
				if ((*it)->x == p0->x && (*it)->y == p0->y && (*it)->z == p0->z ||
					(*it)->x == p0_1.x && (*it)->y == p0_1.y && (*it)->z == p0_1.z)
					continue;

				Point n2 = cross(Point{ p0->x - (*it)->x, p0->y - (*it)->y, p0->z - (*it)->z }, Point{ p0_1.x - (*it)->x, p0_1.y - (*it)->y, p0_1.z - (*it)->z });

				n2.norm();

				double ang_temp = angle(n2, n1);

				if (ang >= ang_temp) {
					ang = ang_temp;
					p0_2 = Point{ (*it)->x, (*it)->y, (*it)->z };
				}
			}

			Face triangle_aux = Face{ &p0_2, &p0_1, p0 };

			global_triangleList.push_back(&triangle_aux);
			ftriangleList.push_back(&triangle_aux);
			triangleList.push_back(&triangle_aux);

			// First Tetra ----------------------------
			Face* t0 = ftriangleList.front();
			Face f0 = *ftriangleList.front();

			double max_ang = 0.0;
			Point* p = NULL;

			for (std::vector<Point*>::iterator it = pointvector.begin(); it != pointvector.end(); it++) {
				if (t0->p1->x == (*it)->x && t0->p1->y == (*it)->y && t0->p1->z == (*it)->z ||
					t0->p2->x == (*it)->x && t0->p2->y == (*it)->y && t0->p2->z == (*it)->z ||
					t0->p3->x == (*it)->x && t0->p3->y == (*it)->y && t0->p3->z == (*it)->z) {
					continue;
				}
				double ang = f0.solidAngle(&f0, (*it), pointvector, true);
				if (max_ang <= ang) {
					max_ang = ang;
					p = (*it);
				}
			}

			Face t0_1 = Face{ t0->p3, p, t0->p1 };
			Face t0_2 = Face{ t0->p2, t0->p1, p };
			Face t0_3 = Face{ t0->p3, t0->p2, p };

			addTriangle(&t0_1);
			addTriangle(&t0_2);
			addTriangle(&t0_3);

			ftriangleList.erase(remove(ftriangleList.begin(), ftriangleList.end(), t0), ftriangleList.end());

			Tetra tetra(t0->p1, t0->p2, t0->p3, p);
			tetraList.push_back(&tetra);
			t.tetras.push_back(tetra);

			// Main Loop ----------------------------
			int k = 0;
			while (ftriangleList.size() > 0) {
				k++;
				Face* n0 = ftriangleList.front();
				Face s0 = *ftriangleList.front();

				double max_ang2 = 0.0;
				Point* pt = NULL;

				for (std::vector<Point*>::iterator it = pointvector.begin(); it != pointvector.end(); it++) {
					if (n0->p1->x == (*it)->x && n0->p1->y == (*it)->y && n0->p1->z == (*it)->z ||
						n0->p2->x == (*it)->x && n0->p2->y == (*it)->y && n0->p2->z == (*it)->z ||
						n0->p3->x == (*it)->x && n0->p3->y == (*it)->y && n0->p3->z == (*it)->z) {
						continue;
					}

					double ang2 = n0->OsolidAngle(&s0, (*it), pointvector, true);

					if (max_ang2 > ang2) {
						max_ang2 = ang2;
						pt = (*it);
					}
				}

				if (pt == NULL) {
					ftriangleList.erase(remove(ftriangleList.begin(), ftriangleList.end(), n0), ftriangleList.end());
					continue;
				}

				Face t0_1n = Face{ n0->p1, pt, n0->p3 };
				Face t0_2n = Face{ n0->p1, n0->p2, pt };
				Face t0_3n = Face{ n0->p2, n0->p3, pt };

				if (!interceptTriangle(&t0_1n) && !interceptTriangle(&t0_2n) && !interceptTriangle(&t0_3n)) {
					addTriangle(new Face { n0->p1, pt, n0->p3 });
					addTriangle(new Face{ n0->p1, n0->p2, pt });
					addTriangle(new Face{ n0->p2, n0->p3, pt });

					Tetra tetran(n0->p1, n0->p2, n0->p3, pt);
					tetraList.push_back(&tetran);

					t.tetras.push_back(tetran);
				}

				ftriangleList.erase(remove(ftriangleList.begin(), ftriangleList.end(), n0), ftriangleList.end());
			}
			// Main Loop ----------------------------
			polygon.push_back(t);

			pointvector.clear();
			triangleList.clear();
			ftriangleList.clear();
			tetraList.clear();

			glLineWidth(1.0f);
			int j = 0;
			for (Tetra& tt : t.tetras) {
				int k = 0;
				j++;
				for (const Face& t : tt.faces) {

					k++;
					srand(95);
					glColor3f(((rand() % 1000) * 0.001f) * fabs(dot(t.normal(), Point{ 1, 1, 1 }) * 0.7f + 0.3f), ((rand() % 1000) * 0.001f) * fabs(dot(t.normal(), Point{ 1, 1, 1 })), ((rand() % 1000) * 0.001f) * fabs(dot(t.normal(), Point{ 1, 1, 1 })));
					glBegin(GL_TRIANGLES);

					glVertex3f(t.p1->x + 0.1 * tt.center().x, t.p1->y + 0.1 * tt.center().y, t.p1->z + 0.1 * tt.center().z);
					glVertex3f(t.p2->x + 0.1 * tt.center().x, t.p2->y + 0.1 * tt.center().y, t.p2->z + 0.1 * tt.center().z);
					glVertex3f(t.p3->x + 0.1 * tt.center().x, t.p3->y + 0.1 * tt.center().y, t.p3->z + 0.1 * tt.center().z);
					glEnd();
				}
			}
		}
	}