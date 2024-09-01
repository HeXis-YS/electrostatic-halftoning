#include "cv.hpp"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define possibility(x, y) __builtin_expect_with_probability(x, 1, y)
#define likely(x) __builtin_expect(x, 1)
#define unlikely(x) __builtin_expect(x, 0)

static double rand_double() {
	return (double)rand() / ((double)RAND_MAX + 1.0);
}

int ElectrostaticHalftoning2010(struct CMat src,
								struct CMat *dst,
								int max_iterations,
								int enable_initial_charge,
								int enable_gridforce,
								int enable_shake,
								int enable_early_stop,
								int enable_debug) {
	//////////////////////////////////////////////////////////////////////////
	///// exceptions
	max_iterations = (max_iterations > 0) ? max_iterations : 8;
	enable_initial_charge = enable_initial_charge ? 1 : 0;
	enable_gridforce = enable_gridforce ? 1 : 0;
	enable_shake = enable_shake ? 1 : 0;
	enable_early_stop = enable_early_stop ? enable_early_stop : 0;
	enable_debug = enable_debug ? 1 : 0;
	printf("Max iterations = %d\n", max_iterations);
	printf("Initial charge = %s\n", enable_initial_charge ? "Enabled" : "Disabled");
	printf("Grid Force = %s\n", enable_gridforce ? "Enabled" : "Disabled");
	printf("Shake = %s\n", enable_shake ? "Enabled" : "Disabled");
	printf("Early stop = ");
	enable_early_stop ? printf("%d iterations.\n", enable_early_stop) : printf("Disabled.\n");
	if (enable_shake && max_iterations <= 64) {
		printf("Error: max_iterations > 64, when enable_shake = 1\n");
		return 1;
	}
	int rows = src.rows;
	int cols = src.cols;
	int pixel_count = rows * cols;
	double *image_in = (double *)malloc(sizeof(double) * pixel_count);
	unsigned char *image_dst = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);
	dst->rows = rows;
	dst->cols = cols;
	dst->data = image_dst;

	//////////////////////////////////////////////////////////////////////////
	///// Initialization
	const int color_depth = 2;
	const int pix_level_max = (1 << color_depth) - 1;
	const double particle_charge = 1.0 / (double)pix_level_max;
	const double particle_charge_2 = 1.0 / (double)(pix_level_max * pix_level_max);
	printf("%f\n", particle_charge);
	double pix_level_float[128] = {0.0};
	unsigned char pix_level_array[128] = {0};
	for (int i = 1; i <= pix_level_max; i++) {
		pix_level_float[i] = (double)i / (double)pix_level_max;
		pix_level_array[i] = (unsigned char)(pix_level_float[i] * 255.0 + 0.5);
	}
	int particle_count = 0;
	memset(image_dst, pix_level_max, sizeof(unsigned char) * pixel_count);
	for (int p = 0; p < pixel_count; p++) {
		int tmp = 255 - src.data[p];
		image_in[p] = (double)tmp / 255.0;
		particle_count += tmp;
	}
	particle_count = particle_count * pix_level_max / 255;
	printf("The number of black pixel(charge) = %d\n", particle_count);

	//////////////////////////////////////////////////////////////////////////
	///// Initialize the Particle's position
	double *particle_Y = (double *)malloc(sizeof(double) * particle_count);
	double *particle_X = (double *)malloc(sizeof(double) * particle_count);
	for (int particle = 0; particle < particle_count;) {
		int rand_Y = rand() % rows;
		int rand_X = rand() % cols;
		int p = rand_Y * cols + rand_X;
		if (enable_initial_charge && rand() % 256 <= src.data[p]) {
			continue;
		}
		image_dst[p] = (image_dst[p] > 0) ? (image_dst[p] - 1) : 0;
		particle_Y[particle] = (double)rand_Y + rand_double();
		particle_X[particle] = (double)rand_X + rand_double();
		particle++;
	}
	if (enable_debug) {
		for (int p = 0; p < pixel_count; p++) {
			image_dst[p] = pix_level_array[image_dst[p]];
		}
		cv_imwrite(".\\output\\0.bmp", *dst);
	}

	//////////////////////////////////////////////////////////////////////////
	///// process
	double *particle_Y_last = (double *)malloc(sizeof(double) * particle_count);
	double *particle_X_last = (double *)malloc(sizeof(double) * particle_count);
	double *distance_X = (double *)malloc(sizeof(double) * cols);
	double *distance_X_2 = (double *)malloc(sizeof(double) * cols);
	double *force_Y_array = (double *)malloc(sizeof(double) * particle_count);
	double *force_X_array = (double *)malloc(sizeof(double) * particle_count);
	// Shake
	int shake = 0;
	double shake_tmp;
	double shake_tmp1;
	if (enable_shake) {
		shake_tmp = log10((double)max_iterations) / log10(1024.0) - 0.6;
		shake_tmp1 = 0.0;
	}
	// Time step
	double time_step = 0.1;
	double mean;
	// Early stop
	unsigned char *image_last = NULL;
	int early_stop_counter = 0;
	if (enable_early_stop) {
		image_last = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);
	}
	if (enable_debug) {
		printf("Time step = %f\n", time_step);
	}
	for (int current_iteration = 1; current_iteration <= max_iterations; current_iteration++) {
		printf("Iteration %d\n", current_iteration);
		memcpy(particle_Y_last, particle_Y, sizeof(double) * particle_count);
		memcpy(particle_X_last, particle_X, sizeof(double) * particle_count);
		memset(force_Y_array, 0, sizeof(double) * particle_count);
		memset(force_X_array, 0, sizeof(double) * particle_count);
		if (enable_shake) {
			if (possibility(current_iteration % 10 == 0, 0.1)) {
				shake_tmp1 = shake_tmp * exp(current_iteration / 1000.0);
				shake = 1;
			} else {
				shake = 0;
			}
		}
		if (enable_debug) {
			mean = 0.0;
		}
		if (enable_early_stop) {
			memcpy(image_last, image_dst, sizeof(unsigned char) * pixel_count);
		}
		memset(image_dst, pix_level_max, sizeof(unsigned char) * pixel_count);
		for (int current_particle = 0; current_particle < particle_count; current_particle++) {
			double force_Y = 0.0;
			double force_X = 0.0;
			double particle_Y_current = particle_Y_last[current_particle];
			double particle_X_current = particle_X_last[current_particle];
			double tmp;

			// Attraction
			for (int x = 0; x < cols; x++) {
				tmp = x + 0.5 - particle_X_current;
				distance_X[x] = tmp;
				distance_X_2[x] = tmp * tmp;
			}
			for (int y = 0, p = 0; y < rows; y++) {
				double distance_Y = y + 0.5 - particle_Y_current;
				double distance_Y_2 = distance_Y * distance_Y;
				for (int x = 0; x < cols; x++, p++) {
					if (unlikely(image_in[p] == 0.0)) {
						continue;
					}
					tmp = distance_Y_2 + distance_X_2[x];
					tmp = (tmp <= 0.1) ? (20.0 - 100.0 * tmp) : (1.0 / tmp);
					tmp *= image_in[p] * particle_charge;
					force_Y += distance_Y * tmp;
					force_X += distance_X[x] * tmp;
				}
			}

			// Repulsion
			for (int particle = current_particle + 1; particle < particle_count; particle++) {
				double distance_Y = particle_Y_last[particle] - particle_Y_current;
				double distance_X = particle_X_last[particle] - particle_X_current;
				tmp = distance_Y * distance_Y + distance_X * distance_X;
				tmp = (tmp <= 0.1) ? (20.0 - 100.0 * tmp) : (1.0 / tmp);
				tmp *= particle_charge_2;
				double repulsion_force = distance_Y * tmp;
				force_Y -= repulsion_force;
				force_Y_array[particle] += repulsion_force;
				repulsion_force = distance_X * tmp;
				force_X -= repulsion_force;
				force_X_array[particle] += repulsion_force;
			}

			// Add GridForce to find discrete particle locations
			if (enable_gridforce) {
				double grid_distance_Y = particle_Y_current - (int)particle_Y_current;
				double grid_distance_X = particle_X_current - (int)particle_X_current;
				grid_distance_Y = possibility(grid_distance_Y <= 0.5, 0.5) ? -grid_distance_Y : 1 - grid_distance_Y;
				grid_distance_X = possibility(grid_distance_X <= 0.5, 0.5) ? -grid_distance_X : 1 - grid_distance_X;
				tmp = grid_distance_Y * grid_distance_Y + grid_distance_X * grid_distance_X;
				if (likely(tmp != 0.0)) {
					tmp = sqrt(tmp);
					tmp = 3.5 / (tmp + 10000.0 * pow(tmp, 9.0));
					force_Y += grid_distance_Y * tmp;
					force_X += grid_distance_X * tmp;
				}
			}

			force_Y += force_Y_array[current_particle];
			force_X += force_X_array[current_particle];
			force_Y_array[current_particle] = force_Y;
			force_X_array[current_particle] = force_X;

			if (enable_debug) {
				mean += sqrt(force_Y * force_Y + force_X * force_X);
			}
		}
		if (enable_debug) {
			mean /= (double)particle_count;
			printf("Mean force = %f\n", mean);
		}
		for (int current_particle = 0; current_particle < particle_count; current_particle++) {
			double particle_Y_current = particle_Y_last[current_particle] + force_Y_array[current_particle] * time_step;
			double particle_X_current = particle_X_last[current_particle] + force_X_array[current_particle] * time_step;

			// Shake
			if (shake) {
				particle_Y_current += shake_tmp1;
				particle_X_current += shake_tmp1;
			}

			// Result (new position of particles)
			particle_Y_current -= floor(particle_Y_current / (double)rows) * (double)rows;
			particle_X_current -= floor(particle_X_current / (double)cols) * (double)cols;
			particle_Y[current_particle] = particle_Y_current;
			particle_X[current_particle] = particle_X_current;

			// Output
			int p = (int)particle_Y_current * cols + (int)particle_X_current;
			image_dst[p] = (image_dst[p] > 0) ? (image_dst[p] - 1) : 0;
		}
		if (enable_debug && shake) {
			printf("Shake performed (%f)\n", shake_tmp1);
		}
		for (int p = 0; p < pixel_count; p++) {
			image_dst[p] = pix_level_array[image_dst[p]];
		}
		if (enable_early_stop) {
			if (memcmp(image_dst, image_last, sizeof(unsigned char) * pixel_count) == 0) {
				early_stop_counter++;
				printf("Result unchanged for %d iterations.\n", early_stop_counter);
				if (early_stop_counter >= enable_early_stop) {
					printf("Early stop.\n");
					break;
				}
			} else {
				early_stop_counter = 0;
			}
		}
		if (enable_debug) {
			char out_file[50];
			sprintf(out_file, ".\\output\\%d.bmp", current_iteration);
			cv_imwrite(out_file, *dst);
		}
	}

	free(image_in);
	free(image_last);
	free(particle_Y);
	free(particle_X);
	free(force_Y_array);
	free(force_X_array);
	free(particle_Y_last);
	free(particle_X_last);
	free(distance_X);
	free(distance_X_2);

	return 0;
}
