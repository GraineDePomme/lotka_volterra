#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A 1
#define B 0.5
#define C 0.5
#define D 2

#define T0 0.0
#define TF 30.0
#define STEPS_COUNT 1000
#define H (TF - T0) / STEPS_COUNT


double dx_dt(double x, double y) {
    return A * x - B * x * y;
}


double dy_dt(double x, double y) {
    return C * x * y - D * y;
}


void next_step_rk4(double* x, double* y) {
    double k1_x = H * dx_dt(*x, *y);
    double k1_y = H * dy_dt(*x, *y);

    double k2_x = H * dx_dt(*x + k1_x / 2, *y + k1_y / 2);
    double k2_y = H * dy_dt(*x + k1_x / 2, *y + k1_y / 2);

    double k3_x = H * dx_dt(*x + k2_x / 2, *y + k2_y / 2);
    double k3_y = H * dy_dt(*x + k2_x / 2, *y + k2_y / 2);

    double k4_x = H * dx_dt(*x + k3_x, *y + k3_y);
    double k4_y = H * dy_dt(*x + k3_x, *y + k3_y);

    *x += (k1_x + 2*k2_x + 2*k3_x + k4_x) / 6;
    *y += (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6;
}


int main() {
    double t = T0;
    double x = 2;
    double y = 2;

    printf("t,x,y\n");

    for (t = T0 ; t <= TF+H ; t += H) {
        printf("%f,%f,%f\n", t, x, y);
        next_step_rk4(&x, &y);
    }

    return 0;
}