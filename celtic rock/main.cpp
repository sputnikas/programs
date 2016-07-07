#include <cstdio>
#include <cstdlib>

#include "Function.h"
#include "Vector3.h"
#include "Matrix3.h"

const double g = 9.81;

double ea, eb, ec;
double m, M;
double xa, ya;
double k;
double Mc;
Matrix J;
double dt;

struct QP {
    double psi;
    double theta;
    double phi;
    Vector omega;
};

int N;
QP *allQP;

QP operator + (QP a, QP b) {
    return {a.psi + b.psi, a.theta + b.theta, a.phi + b.phi, a.omega + b.omega};
}

QP operator * (double a, QP b) {
    return {a*b.psi, a*b.theta, a*b.phi, a*b.omega};
}

QP operator * (QP b, double a) {
    return {a*b.psi, a*b.theta, a*b.phi, a*b.omega};
}

QP right(double t, QP qp) {
    QP r;
    double sth = sin(qp.theta), cth = cos(qp.theta);
    double sph = sin(qp.phi), cph = cos(qp.phi);
    r.psi   = (sth != 0.) ? (qp.omega.x*sph + qp.omega.y*cph) / sth : 0.;
    r.theta = qp.omega.x*cph - qp.omega.y*sph;
    r.phi   = (cth != 0.) ? qp.omega.z - (qp.omega.x*sph + qp.omega.y*cph) / cth : qp.omega.z;

    double NA = - 1/sqrt(ea*ea*sth*sth*sph*sph + eb*eb*sth*sth*cph*cph + ec*ec*cth*cth);
    Vector R = Vector(ea*ea*sth*sph, - eb*eb*sth*cph, ec*ec*cth)*fabs(NA);
    Vector A = Vector(sth*sph, - sth*cph, cth)*NA;
    Vector T = A*R;
    Matrix Jd = Matrix(Mc/NA/NA*T.x*T.x, Mc/NA/NA*T.x*T.y, Mc/NA/NA*T.x*T.z,
                       Mc/NA/NA*T.x*T.y, Mc/NA/NA*T.y*T.y, Mc/NA/NA*T.y*T.z,
                       Mc/NA/NA*T.z*T.x, Mc/NA/NA*T.y*T.z, Mc/NA/NA*T.z*T.z);
    Matrix Jf = J - Jd;
//    printf("Jf = "); to_console(Jf); printf("\n");
//    printf("T = "); to_console(T); printf("\n");
//    printf("R = "); to_console(R); printf("\n");
//    printf("psi = %e\n", qp.psi);
//    printf("theta = %e\n", qp.theta);
//    printf("phi = %e\n", qp.phi);
    Vector Mfr = - Mc*g*k*ea*eb*qp.omega;
    //r.omega = Jf.inversed()*(Mfr);
    r.omega = Jf.inversed()*(Mc*g*T/NA + Mfr - qp.omega*(J*qp.omega) + Mc*T*dot(qp.omega, T)/NA/NA);
    return r;
}

double energy(QP qp) {
    return 0.5*dot(J*qp.omega, qp.omega);
}

//////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////

#include <GL/gl.h>
#include <GL/glut.h>

enum
{
	MOUSE_LEFT_BUTTON = GLUT_LEFT_BUTTON,
	MOUSE_MIDDLE_BUTTON = GLUT_MIDDLE_BUTTON,
	MOUSE_RIGHT_BUTTON = GLUT_RIGHT_BUTTON,
	MOUSE_SCROLL_UP = 3,
	MOUSE_SCROLL_DOWN = 4
};

int nt_opengl;
double angle_x = 30.0;
double angle_y = 70.0;
double delta_angle = 1;
double center_scene = -8.0;
double x_origin;
double y_origin;
int gl_width = 800;
int gl_height = 600;
double scale = 1;//./pow(ea*eb*ec, 1./3);
int ngrid = 20;
double gscale;

const GLfloat light_ambient[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
const GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position[] = { 2.0f, 5.0f, 5.0f, 0.0f };

const GLfloat mat_ambient[]    = { 0.7f, 0.7f, 0.7f, 1.0f };
const GLfloat mat_diffuse[]    = { 0.8f, 0.8f, 0.8f, 1.0f };
const GLfloat mat_specular[]   = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat high_shininess[] = { 100.0f };

void glVertex3f(Vector r) {
    glVertex3f(r.x, r.y, r.z);
};

void display() {
    glMatrixMode(GL_PROJECTION);            // Действия будут производиться с матрицей проекции
    glLoadIdentity();                       // Текущая матрица (проекции) сбрасывается на единичную
    glFrustum(-scale, scale, -scale, scale, 1, 10);         // Устанавливается перспективная проекция

    glMatrixMode(GL_MODELVIEW);             // Действия будут производиться с матрицей модели
    glLoadIdentity();                       // Текущая матрица (модели) сбрасывается на единичную
    glTranslatef(0.0, 0.0, center_scene);           // Система координат переносится на 8 единиц вглубь сцены
    glRotatef(angle_x, 1.0, 0.0, 0.0);         // и поворачивается на 30 градусов вокруг оси x,
    glRotatef(angle_y, 0.0, 1.0, 0.0);         // а затем на 70 градусов вокруг оси y

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  // Очищается буфер кадра и буфер глубины
    glPushMatrix();                                      // Запоминается матрица модели

    glColor4f(1.0f, 0.0f, 0.0f, 0.0f);                         // Задается текущий цвет примитивов

    double sth = sin(allQP[nt_opengl].theta), cth = cos(allQP[nt_opengl].theta);
    double sph = sin(allQP[nt_opengl].phi), cph = cos(allQP[nt_opengl].phi);

    glTranslated(0, sqrt(ea*ea*sth*sth*sph*sph + eb*eb*sth*sth*cph*cph + ec*ec*cth*cth), 0);
    glRotated(allQP[nt_opengl].psi*180/M_PI,   0, 1, 0);
    glRotated(allQP[nt_opengl].theta*180/M_PI, 1, 0, 0);
    glRotated(allQP[nt_opengl].phi*180/M_PI,   0, 1, 0);
    glScaled(ea, eb, ec);

    glutSolidSphere(1., 64, 64);

    glPopMatrix();

    glPushMatrix();

    glColor3f(1.0f, 1.0f, 1.0f);
    gscale = 0.1*sqrt(ea*eb)*1/scale;
    glBegin(GL_LINES);
    for (int i = 0; i<ngrid; i++) {
        glVertex3f(-gscale + i*2*gscale/(ngrid - 1), 0, -gscale);
        glVertex3f(-gscale + i*2*gscale/(ngrid - 1), 0, gscale);
        glVertex3f(-gscale, 0, -gscale + i*2*gscale/(ngrid - 1));
        glVertex3f(gscale, 0, -gscale + i*2*gscale/(ngrid - 1));
    }
    glEnd();

    glPopMatrix();

    glutSwapBuffers();
}

void timer(int value) {
    nt_opengl = (nt_opengl + 20 < N) ? nt_opengl + 20 : 0;
    printf("\r %6d %2.4e %2.4e %2.4e %2.4e",
           nt_opengl,
           allQP[nt_opengl].omega.x,
           allQP[nt_opengl].omega.y,
           allQP[nt_opengl].omega.z,
           energy(allQP[nt_opengl]));
    glutPostRedisplay();
    glutTimerFunc(100, timer, 0);
}

void reshape(int width, int height) {
    if (height == 0) {
        height = 1;
    }
    gl_width = width;
    gl_height = height;

    glViewport(0, 0, width, height);                     // Устанавливается область просмотра

    glMatrixMode(GL_PROJECTION);            // Действия будут производиться с матрицей проекции
    glLoadIdentity();                       // Текущая матрица (проекции) сбрасывается на единичную
    glFrustum(-scale, scale, -scale, scale, -10, 10);         // Устанавливается перспективная проекция

    glMatrixMode(GL_MODELVIEW);             // Действия будут производиться с матрицей модели
    glLoadIdentity();                       // Текущая матрица (модели) сбрасывается на единичную
    glTranslatef(0.0, 0.0, center_scene);           // Система координат переносится на 8 единиц вглубь сцены
    glRotatef(angle_x, 1.0, 0.0, 0.0);         // и поворачивается на 30 градусов вокруг оси x,
    glRotatef(angle_y, 0.0, 1.0, 0.0);         // а затем на 70 градусов вокруг оси y
}

void mouse_pressed_move(int x, int y) {
    angle_y = (x - x_origin) / gl_width * 360;
    angle_x = (y - y_origin) / gl_width * 360;
}

void mouse_button(int button, int state, int x, int y) {
    switch (button) {
    case MOUSE_LEFT_BUTTON:
        if (state == GLUT_DOWN) {
            x_origin = x;
            y_origin = y;
		}
        break;
    case MOUSE_MIDDLE_BUTTON:
        break;
    case MOUSE_RIGHT_BUTTON:
        break;
    case MOUSE_SCROLL_UP:
        scale += 0.1*scale;
        break;
    case MOUSE_SCROLL_DOWN:
        scale -= 0.1*scale;
        break;
    }
}

int main(int argc, char** argv) {
    ea = 0.1, eb = 0.01, ec = 0.01;
    m = 0.1, M = 0.;
    xa = 0.05, ya = 0.002;
    k = 1;
    dt = 0.0001;
    N = 10000;
    allQP = new QP[N];
    QP qp0 = {0., 0.2, 0., Vector(0., 0., 4)};
    Mc = m + M;
    J = Matrix( 0.2*m*(eb*eb + ec*ec) + M*ya*ya, M*xa*ya, 0,
                M*xa*ya, 0.2*m*(ea*ea + ec*ec) + M*xa*xa, 0,
                0,          0,        0.2*m*(ea*ea + eb*eb));
    QP k1, k2, k3, k4;
    allQP[0] = qp0;
    for (int i = 1; i<N; i++) {
        //dt = 2*M_PI/norm(allQP[i - 1].omega)/10000;
        k1 = right(i*dt, allQP[i - 1])*dt;
        k2 = right(i*dt, allQP[i - 1] + 0.5*k1)*dt;
        k3 = right(i*dt, allQP[i - 1] + 0.5*k2)*dt;
        k4 = right(i*dt, allQP[i - 1] + k3)*dt;
        allQP[i] = allQP[i - 1] + 1./6*(k1 + 2*k2 + 2*k3 + k4);
    }

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(gl_width, gl_height);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Celtic rock");

    glClearColor(0.2f, 0.5f, 0.75f, 1.0f);

    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);

    glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);

    glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);

    glutTimerFunc(100, timer, 0);
    glutMouseFunc(mouse_button);
	glutMotionFunc(mouse_pressed_move);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMainLoop();

    delete [] allQP;
    return 0;
}
