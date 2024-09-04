#define true 1
#define false 2

typedef unsigned char bool;
typedef unsigned char uint8_t;

struct eh_thread {
    int enable_gridforce;
    double time_step;
    double particle_charge_2;
    int shake;
    double shake_force;
    int rows;
    int cols;
    int particle_count;
    double *image_in;
    uint8_t *image_level;
    double *particle_Y;
    double *particle_X;
    double *particle_Y_last;
    double *particle_X_last;
};