#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#if defined(_WIN64)||defined(_WIN32)
    #include <windows.h>
#endif // _WIN32, _WIN64
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include "main.h"
#include "math_algorithms.h"

#include <chrono>

Params h;
Particle *particles[N_COUNTS];
double xx = 0.;

void displayCall() {
    int nt = (int) trunc(glutGet(GLUT_ELAPSED_TIME) / 1000.0);
    while (nt >= N_COUNTS)
        nt -= N_COUNTS;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    gluLookAt(0., 0., xx - h.a, 0., 0., 0., 0., 1., 0.);
    glMatrixMode(GL_MODELVIEW);

    glColor3f(1,0,0);

    glPushMatrix();

    printf("\r%d", nt);
    //glTranslated(xx, 0., 0.1);
    //xx += 0.1;
    //glutSolidSphere(0.05, 32, 16);
    for (int i = 0; i<h.n_particle; i++){
        glTranslated(particles[nt][i].r.x, particles[nt][i].r.y, particles[nt][i].r.z);
        glutSolidSphere(particles[nt][i].R, 32, 16);
    }
    glPopMatrix();

    glutSwapBuffers();
}

void key(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 'q':
            exit(0);
            break;
        case 'w':
            xx += 0.2;
            break;
        case 's':
            xx -= 0.2;
            break;
    }
}

int main(int argc, char *argv[]){
    #if ((defined(_WIN64)||defined(_WIN32))&&(defined(RUS)))
        system("chcp 65001");
    #endif // __WIN__

    using namespace std::chrono;
    steady_clock::time_point t1 = steady_clock::now();

    init("test.ini");
    //for (int i = 0; i<N_COUNTS; i++){
    //    calculate(i);
    //    printf("\r%d", i);
    //}
    //printf("\n");

    //to_bfile("test.particles");
    from_bfile("test.particles");

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(300, 200);
    glutCreateWindow("Hello World!");
    glOrtho(-h.a, h.a, -h.a, h.a, h.a, -h.a);
    //glFrustum(-640./480, 640./480, -1.0, 1.0, 2.0, 100.0);
    glClearColor(1., 1., 1., 1.);
    glutDisplayFunc(displayCall);
    glutIdleFunc(displayCall);
    glutMainLoop();

    free();

    steady_clock::time_point t2 = steady_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    printf(SECS, time_span.count());

    return 0;
}
