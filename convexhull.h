#pragma once
#include <algorithm>
#include <vector>
#include "vectorops.h"
#define M_PI           3.14159265358979323846  /* pi */

bool Place(Point& A, Point& B, Point& C, Point P) {

    double det = (B.y - C.y) * (A.x - C.x) + (C.x - B.x) * (A.y - C.y);
    double factor_alpha = (B.y - C.y) * (P.x - C.x) + (C.x - B.x) * (P.y - C.y);
    double factor_beta = (C.y - A.y) * (P.x - C.x) + (A.x - C.x) * (P.y - C.y);
    double alpha = factor_alpha / det;
    double beta = factor_beta / det;
    double gamma = 1.0 - alpha - beta;

    bool In = false;

    if (((A.x == P.x) && (A.y == P.y)) || ((B.x == P.x) && (B.y == P.y)) || ((C.x == P.x) && (C.y == P.y)))
        In = true; 
    if ((alpha == 0) || (beta == 0) || (gamma == 0))
        In = true; 
    if (((0 < alpha) && (alpha < 1)) && ((0 < beta) && (beta < 1)) && ((0 < gamma) && (gamma < 1)))
        In = true; 

    return In;
}

struct Edge {
    Point* p1, * p2;

    bool equal(Point* v) const {
        bool h = p1->x == v->x && p1->y == v->y && p1->z == v->z;
        return h;
    }

    inline bool operator==(const Edge& e) const { return equal(e.p1) && equal(e.p2); }

    bool left(Point* other) {

        std::vector<std::vector<double>> m{
       {p1->x, p2->x, other->x},
       {p1->y, p2->y, other->y},
       {1, 1, 1} };

        if (get3x3Determinant(m) < 0) {
            return false;
        }
        else {
            return true;
        }
    }
};

struct EdgeStack {
    vector<Edge> data;

    Edge get() {
        Edge b = data.back();
        data.pop_back();

        return b;
    }
    void put(Edge e) {
        data.push_back(e);
    }

    int findIndex(const vector<Edge>& arr, Edge item) {

        for (int i = 0; i < arr.size(); i++) {
            if ((arr[i].p1->x == item.p1->x && arr[i].p1->y == item.p1->y && arr[i].p1->z == item.p1->z)
                && (arr[i].p2->x == item.p2->x && arr[i].p2->y == item.p2->y && arr[i].p2->z == item.p2->z)
                || (arr[i].p2->x == item.p1->x && arr[i].p2->y == item.p1->y && arr[i].p2->z == item.p1->z)
                && (arr[i].p1->x == item.p2->x && arr[i].p1->y == item.p2->y && arr[i].p1->z == item.p2->z)) {
                return i;
            }
        }

        return -1;
    }

    void puts(Edge e) {

        int index = findIndex(data, e);
        if (index == -1) {
            data.push_back(e);
        }
        else {
            data.erase(data.begin() + index);
        }
    }

    bool isEmpty() {
        return data.empty();
    }
};

struct Face {
    union {
        struct { Point* p1, * p2, * p3; };
        Point* p[3];
    };

    inline Point normal() const {
        return cross(Point{ p2->x - p1->x, p2->y - p1->y, p2->z - p1->z }, Point{ p3->x - p1->x, p3->y - p1->y, p3->z - p1->z }).unit();
    }

    inline Point center() const {
        return Point{ (p1->x + p2->x + p3->x) / 3.0,  (p1->y + p2->y + p3->y) / 3.0, (p1->z + p2->z + p3->z) / 3.0 };
    }


    bool outside(const Face& f, const Point& p) {
        return dot(f.normal(), (Point{ p.x - f.center().x, p.y - f.center().y, p.z - f.center().z }).unit()) <= 0;
    }

    bool accordingToNormal(const Face* f, const Point* p, bool insideTest) const {
        Point n = cross(Point{ f->p2->x - f->p1->x, f->p2->y - f->p1->y, f->p2->z - f->p1->z }, Point{ f->p3->x - f->p1->x, f->p3->y - f->p1->y, f->p3->z - f->p1->z });
        Point v = Point{ p->x - f->p1->x,  p->y - f->p1->y, p->z - f->p1->z };

        n.norm();
        v.norm();

        return (insideTest) ? dot(n, v) > -0.0001 : dot(n, v) > 0.0001;
    }

    double solidAngle(const Face* f, const Point* p, vector<Point*> pointvector, bool insideTest) const {

        Face test = Face{ f->p1, f->p2, f->p3 };

        if (Place(*f->p1, *f->p2, *f->p3, *p)) {
            return fabs(2.0 * M_PI);
        }

        Point* p1 = f->p1;
        Point* p2 = f->p2;
        Point* p3 = f->p3;

        Point a{ p1->x - p->x, p1->y - p->y, p1->z - p->z };
        Point b{ p2->x - p->x, p2->y - p->y, p2->z - p->z };
        Point c{ p3->x - p->x, p3->y - p->y, p3->z - p->z };

        double an = a.normd();
        double bn = b.normd();
        double cn = c.normd();

        double den = an * bn * cn + dot(a, b) * cn + dot(b, c) * an + dot(c, a) * bn;

        if (fabs(den) < 0.0001) {
            return 0.0;
        }

        double num = dot(a, cross(c, b));

        double ang = 2.0 * atan2(num, den);

        if (fabs(ang) < 0.0001) {
            ang = 0.0;
        }

        ang = test.accordingToNormal(f, p, insideTest) ?
            ((ang < 0.0) ? 2.0 * M_PI + ang : ang) :
            ((ang > 0.0) ? -2.0 * M_PI + ang : ang);

        return fabs(ang);
    }

    double OsolidAngle(const Face* f, const Point* p, vector<Point*> pointvector, bool insideTest) const {

        Face test = Face{ f->p1, f->p2, f->p3 };

        Point* p1 = f->p1;
        Point* p2 = f->p2;
        Point* p3 = f->p3;

        Point a{ p1->x - p->x, p1->y - p->y, p1->z - p->z };
        Point b{ p2->x - p->x, p2->y - p->y, p2->z - p->z };
        Point c{ p3->x - p->x, p3->y - p->y, p3->z - p->z };

        double an = a.normd();
        double bn = b.normd();
        double cn = c.normd();

        double den = an * bn * cn + dot(a, b) * cn + dot(b, c) * an + dot(c, a) * bn;

        if (fabs(den) < 0.0001) {
            return 0.0;
        }

        double num = dot(a, cross(c, b));

        double ang = 2.0 * atan2(num, den);

        if (fabs(ang) < 0.0001) {
            ang = 0.0;
        }

        ang = test.accordingToNormal(f, p, insideTest) ?
            ((ang < 0.0) ? 2.0 * M_PI + ang : ang) :
            ((ang > 0.0) ? -2.0 * M_PI + ang : ang);

        return ang;
    }


    bool operator==(const Face& t) {
        return this->has(t.p1) && this->has(t.p2) && this->has(t.p3);
    }

    bool has(Point* v) {
        return (p1 == v || p2 == v || p3 == v);
    }
};

struct Hull {
    vector<Point*> vertices;
    vector<Face> faces;

    void clear() { vertices.clear(); faces.clear(); }
    void addv(Point* v) { vertices.push_back(v); }
};

vector<Hull> convexhull(vector<Hull>& objects) {

    vector<Hull> hulls;
    for (Hull& h : objects) {
        Hull hull;
        hull.vertices = h.vertices;

        EdgeStack es = EdgeStack();

        srand(1);
        // Random noise
        for (int k = 0; k < hull.vertices.size(); k++) {
            hull.vertices[k]->x = hull.vertices[k]->x + ((float)rand() / RAND_MAX) / 10000;
            hull.vertices[k]->y = hull.vertices[k]->y + ((float)rand() / RAND_MAX) / 10000;
            hull.vertices[k]->z = hull.vertices[k]->z + ((float)rand() / RAND_MAX) / 10000;
        }

        // Search 2D
        double minY = 1000000000000;
        Point* firstpoint = hull.vertices[0];
        Point* secondpoint = hull.vertices[0];
        int p1 = 0;
        for (int i = 0; i < hull.vertices.size(); i++) {

            if (hull.vertices[i]->y < minY) {
                minY = hull.vertices[i]->y;
                firstpoint = hull.vertices[i];
                p1 = i;
            }
        }

        int p2 = 0;
        for (int j = 0; j < hull.vertices.size(); j++) {
            Edge testEdge = Edge{ firstpoint, secondpoint };
            if (testEdge.left(hull.vertices[j])) {
                secondpoint = hull.vertices[j];
            }
        }

        Edge e = Edge{ firstpoint, secondpoint };

        es.data.push_back(e);
        es.data.push_back(Edge{ e.p2, e.p1 });

        vector<Face> hullFaces;
        while (!es.isEmpty()) {
            e = es.get();
            // Search
            int i;

            for (i = 0; ((hull.vertices[i]->x == e.p1->x && hull.vertices[i]->y == e.p1->y && hull.vertices[i]->z == e.p1->z))
                || ((hull.vertices[i]->x == e.p2->x && hull.vertices[i]->y == e.p2->y && hull.vertices[i]->z == e.p2->z)); i++) {
            }
            Point* candidate = hull.vertices[i];
            int c = 0;
            Face candh = Face{ e.p1, e.p2, candidate };
            for (i++; i < hull.vertices.size(); i++) {
                bool b = !candh.outside(candh, *hull.vertices[i]);
                if ((hull.vertices[i]->x != e.p1->x || hull.vertices[i]->y != e.p1->y || hull.vertices[i]->z != e.p1->z)
                    && (hull.vertices[i]->x != e.p2->x || hull.vertices[i]->y != e.p2->y || hull.vertices[i]->z != e.p2->z)
                    && !candh.outside(candh, *hull.vertices[i])) {
                    candidate = hull.vertices[i];
                    candh = Face{ e.p1, e.p2, candidate };
                    c = i;
                }
            }
            Point* cand = candidate;
            hullFaces.push_back(Face{ e.p1, e.p2, cand });
            es.puts(Edge{ e.p1, cand });
            es.puts(Edge{ cand, e.p2 });
        }

        for (Face& f : hullFaces) {
            hull.faces.push_back(f);
        }
        hulls.push_back(hull);
    }

    return hulls;
}