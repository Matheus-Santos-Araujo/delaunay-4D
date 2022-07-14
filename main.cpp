#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include <algorithm>
#include<GL/freeglut.h>
#include<chrono>
#include<thread>
#include "convexhull.h"
#include "deulanay.h"

int width = 800, height = 800;
Camera camera;

float resolutionx = width * 1.0f;
float resolutiony = height * 1.0f;
float fovx = 2.0f / width;
float fovy = 2.0f / height;
bool tint, fourd, tri = false;
int t;
int T = 8;

Hull obj;
vector<Hull> objects, hulls;

void printhull(const Hull& o) {

	// Point cloud
	glPointSize(4.0f);
	glBegin(GL_POINTS);
	for (const Point* v : o.vertices) {
		if (v->x > 0.9 && v->x < 1.1 && v->y > -1.1 && v->y < -0.9 && v->z > -1.1 && v->z < -0.9) { glColor3f(0.1f, 0.1f, 0.8f); }
		else { glColor3f(0.8f, 0.1f, 0.1f); }
		glVertex3f(v->x, v->y, v->z);
	}
	glEnd();

	// Faces
	glLineWidth(1.0f);
	for (const Face& t : o.faces) {

		if (tint) {
			glColor3f(0.8, 0.1, 0.1);
			glBegin(GL_TRIANGLES);
			for (int i = 0; i < 3; i++) {
				glVertex3f(t.p[i]->x, t.p[i]->y, t.p[i]->z);
			}
			glEnd();
		}
		Point n = t.normal();
		n.norm();
		glColor3f(0.1f, 0.1f, 0.1f);
		glBegin(GL_LINES);
		for (int i = 0; i < 3; i++) {
			glVertex3f(t.p[i]->x, t.p[i]->y, t.p[i]->z);
			glVertex3f(t.p[(i + 1) % 3]->x, t.p[(i + 1) % 3]->y, t.p[(i + 1) % 3]->z);

			n.x *= -1; n.y *= -1; n.z *= -1;
			glVertex3f(t.p[i]->x, t.p[i]->y, t.p[i]->z);
			glVertex3f(t.p[(i + 1) % 3]->x, t.p[(i + 1) % 3]->y, t.p[(i + 1) % 3]->z);
		}
		glEnd();
	}
}

void loader(string name) {

	string tag;
	int num_objs = -1;
	int a, b, c;
	int A, B, C;
	int i, j, k;
	double x, y, z;

	std::ifstream in(name, std::ios::in);
	if (!in)
	{
		std::cerr << "Cannot open " << name << std::endl;
		exit(1);

	}
	std::string line;
	while (std::getline(in, line))
	{

		if (line.substr(0, 2) == "o ") {
			num_objs++;
			objects.push_back(Hull());
		}

		else if (line.substr(0, 2) == "v ") {
			std::istringstream v(line.substr(2));
			v >> x; v >> y; v >> z;
			Point* p = new Point{ x, y, z };

			obj.addv(p);
			objects[num_objs].addv(p);
		}
	}
}

void disp(void) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	{

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glViewport(0, 0, resolutionx, resolutiony);
		int largura = fovx * resolutionx, altura = fovy * resolutiony;
		glFrustum(-largura / 2, largura / 2, -altura / 2, altura / 2, 1, 100);
		glMatrixMode(GL_MODELVIEW);

		glLoadIdentity();

		gluLookAt(camera.position.x, camera.position.y, camera.position.z,
			camera.view.x, camera.view.y, camera.view.z,
			camera.up.x, camera.up.y, camera.up.z);

		glPointSize(10.0f);
		glColor3f(0, 0, 1);

		glPointSize(4.0f);
		glLineWidth(2.0f);
		glColor3f(1, 1, 1);
		//look at
		glBegin(GL_POINTS);
		glVertex3f(0, 0, 0); // view
		glEnd();
	}

	if (tri == false) {
		if (fourd == true && t < T) {
			t += 1;
			objects.clear();
			obj.clear();
			hulls.clear();
			std::stringstream ss;
			ss << "mywind" << t << ".obj";
			string frame = ss.str();
			cout << frame << "\n";
			loader(frame);
			hulls = convexhull(objects);

			this_thread::sleep_for(chrono::milliseconds(2000));
		}

		for (const Hull& o : hulls) {
			printhull(o);
		}

	}
	else if (tri == true && fourd == true) {
		if (fourd == true && t < T) {
			t += 1;
			objects.clear();
			obj.clear();
			hulls.clear();
			std::stringstream ss;
			ss << "mywind" << t << ".obj";
			string frame = ss.str();
			cout << frame << "\n";
			loader(frame);
			hulls = convexhull(objects);

			Delaunay3D(hulls);

			this_thread::sleep_for(chrono::milliseconds(2000));
		}
	} else {
		Delaunay3D(hulls);
	}

	glutSwapBuffers();
}

void keyboard(unsigned char key, int /*x*/, int /*y*/) {

	switch (key) {

	case(27): exit(0);
	case(' '):
		break;
	case('o'):
		tri = false;
		fourd = false;
		objects.clear();
		// Descomentar modelo .obj 3D desejado
		//loader("mywind1.obj");
		loader("cube.obj");
		hulls = convexhull(objects);
		break;
	case('p'):
		tint = true;
		break;
	case('a'):
		tint = false;
		break;
		// Tema 4D
	case('t'):
		tri = false;
		fourd = true;
		t = 0;
		break;
	case('y'):
		tri = true;
		fourd = true;
		t = 0;
		break;
	case('x'):
		fourd = false;
		// Descomentar modelo .obj 3D desejado
		//loader("mywind1.obj");
		loader("cube.obj");
		hulls = convexhull(objects);
		tri = !tri;
		break;
	}
	glutPostRedisplay();
}

void tick(int m) {

	//Redraw the display
	glutPostRedisplay();

	/*recall timer func in 33 milliseconds
	redraw every 33 milliseconds which is
	approximatly=1000/33=30 frames per seconds
	*/
	glutTimerFunc(33, tick, 0);

}

int lx = 0, ly = 0;
int theButtonState = 0;
int theModifierState = 0;

void motion(int x, int y) {
	int cX = lx - x;
	int cY = ly - y;

	if (cX != 0 || cY != 0) {

		if (theButtonState == GLUT_LEFT_BUTTON) {
			Point v = cross(camera.view - camera.position, camera.up);
			camera.Transform(rotateArbitrary(cY, v));
			camera.Transform(rotateY(cX));
		}

		if (theButtonState == GLUT_RIGHT_BUTTON) { camera.position += camera.position * (cY * 0.01f); }

		lx = x;
		ly = y;
		glutPostRedisplay();
	}
}

void mouse(int b, int state, int x, int y) {
	theButtonState = b;
	theModifierState = glutGetModifiers();
	lx = x;
	ly = y;

	motion(x, y);
}


int main(int argc, char** argv) {

	camera = Camera{ Point{6,2,-2}, Point{0,2,0}, Point{0,1,0} };
	camera.resolution = f2 { width * 1.0f, height * 1.0f };
	camera.fov = f2 { 2.0f / width, 2.0f / height };

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(width, height);
	glutCreateWindow("Fecho Convexo e Tetraedralização");

	glClearColor(0.8f, 0.8f, 0.8f, 0.0f);

	glMatrixMode(GL_PROJECTION);
	glViewport(0, 0, camera.resolution.x, camera.resolution.y);
	int largura = camera.fov.x * camera.resolution.x, altura = camera.fov.y * camera.resolution.y;
	glFrustum(-largura / 2, largura / 2, altura / 2, -altura / 2, 1, 100);
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	glutDisplayFunc(disp);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutKeyboardFunc(keyboard);
	glutTimerFunc(33, tick, 0);
	glutMainLoop();

	cin.get();
	cin.get();
	system("pause");
	return 0;
}