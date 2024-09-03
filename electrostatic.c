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

int electrostatic_halftoning(const char *src_path,
							 const char *dst_path,
							 int color_depth,
							 int max_iterations,
							 int enable_initial_charge,
							 int enable_gridforce,
							 int enable_shake,
							 int enable_early_stop,
							 int enable_debug) {
	/* Exceptions */
	max_iterations = (max_iterations > 0) ? max_iterations : 8;
	enable_initial_charge = enable_initial_charge ? 1 : 0;
	enable_gridforce = enable_gridforce ? 1 : 0;
	enable_shake = enable_shake ? 1 : 0;
	enable_early_stop = enable_early_stop ? enable_early_stop : 0;
	enable_debug = enable_debug ? 1 : 0;
	if (color_depth < 1 || color_depth > 7) {
		printf("Error: Color depth = [1, 7]\n");
	}
	printf("Color depth = %d\n", color_depth);
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

	// Load image
	struct CMat src;
	if (cv_imread(src_path, &src) == 0) {
		return 2;
	}
	const int rows = src.rows;
	const int cols = src.cols;
	const int pixel_count = rows * cols;
	double *image_in = (double *)malloc(sizeof(double) * pixel_count);
	unsigned char *image_dst = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);
	unsigned char *image_level = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);
	struct CMat dst;
	dst.rows = rows;
	dst.cols = cols;
	dst.data = image_dst;

	/* Color Depth */
	const int pixel_level_max = (1 << color_depth) - 1;
	const double particle_charge = 1.0 / (double)pixel_level_max;
	const double particle_charge_2 = 1.0 / (double)(pixel_level_max * pixel_level_max);
	unsigned char pixel_level[128];
	for (int i = 0; i <= pixel_level_max; i++) {
		pixel_level[i] = (255 * i + pixel_level_max / 2) / pixel_level_max;
	}

	/* Particle Count */
	int particle_count = 0;
	for (int p = 0; p < pixel_count; p++) {
		int tmp = 255 - src.data[p];
		image_in[p] = (double)tmp * particle_charge / 255.0;
		particle_count += tmp;
	}
	particle_count = (particle_count * pixel_level_max + 127) / 255;
	printf("Particle count = %d\n", particle_count);

	/* Particle Initialization */
	double *particle_Y = (double *)malloc(sizeof(double) * particle_count);
	double *particle_X = (double *)malloc(sizeof(double) * particle_count);
	memset(image_level, pixel_level_max, sizeof(unsigned char) * pixel_count);
	for (int particle = 0; particle < particle_count;) {
		int rand_Y = rand() % rows;
		int rand_X = rand() % cols;
		int p = rand_Y * cols + rand_X;
		if (enable_initial_charge && (rand() % 257 < src.data[p])) {
			continue;
		}
		image_level[p] -= (image_level[p] > 0) ? 1 : 0;
		particle_Y[particle] = (double)rand_Y + rand_double();
		particle_X[particle] = (double)rand_X + rand_double();
		particle++;
	}

	if (enable_debug) {
		for (int p = 0; p < pixel_count; p++) {
			image_dst[p] = pixel_level[image_level[p]];
		}
		cv_imwrite(".\\output\\0.bmp", dst);
	}

	/* Process */
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
	printf("\n");
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
		if (enable_early_stop) {
			memcpy(image_last, image_level, sizeof(unsigned char) * pixel_count);
		}
		memset(image_level, pixel_level_max, sizeof(unsigned char) * pixel_count);
		for (int current_particle = 0; current_particle < particle_count; current_particle++) {
			double force_Y = 0.0;
			double force_X = 0.0;
			double particle_Y_current = particle_Y_last[current_particle];
			double particle_X_current = particle_X_last[current_particle];
			double tmp;

			// Attraction
			double distance_Y = 0.5 - particle_Y_current;
			tmp = 0.5 - particle_X_current;
			for (int x = 0; x < cols; x++) {
				distance_X[x] = tmp;
				distance_X_2[x] = tmp * tmp;
				tmp += 1.0;
			}
			for (int y = 0, p = 0; y < rows; y++) {
				double distance_Y_2 = distance_Y * distance_Y;
				for (int x = 0; x < cols; x++, p++) {
					if (unlikely(image_in[p] == 0.0)) {
						continue;
					}
					tmp = distance_Y_2 + distance_X_2[x];
					tmp = (tmp <= 0.1) ? (20.0 - 100.0 * tmp) : (1.0 / tmp);
					// image_in[p] * particle_charge is already done in the initialization step
					tmp *= image_in[p];
					force_Y += distance_Y * tmp;
					force_X += distance_X[x] * tmp;
				}
				distance_Y += 1.0;
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

			// For debug only
			force_Y_array[current_particle] = force_Y;
			force_X_array[current_particle] = force_X;

			particle_Y_current = particle_Y_last[current_particle] + force_Y * time_step;
			particle_X_current = particle_X_last[current_particle] + force_X * time_step;

			// Shake
			if (shake) {
				particle_Y_current += shake_tmp1;
				particle_X_current += shake_tmp1;
			}

			// Result (new position of particles)
			particle_Y[current_particle] = particle_Y_current - floor(particle_Y_current / (double)rows) * (double)rows;
			particle_X[current_particle] = particle_X_current - floor(particle_X_current / (double)cols) * (double)cols;

			// Output
			int p = (int)particle_Y_current * cols + (int)particle_X_current;
			image_level[p] -= (image_level[p] > 0) ? 1 : 0;
		}
		if (enable_debug) {
			mean = 0.0;
			for (int current_particle = 0; current_particle < particle_count; current_particle++) {
				double force_Y = force_Y_array[current_particle];
				double force_X = force_X_array[current_particle];
				mean += sqrt(force_Y * force_Y + force_X * force_X);
			}
			mean /= (double)particle_count;
			printf("Mean force = %f\n", mean);
			if (shake) {
				printf("Shake performed (%f)\n", shake_tmp1);
			}
			for (int p = 0; p < pixel_count; p++) {
				image_dst[p] = pixel_level[image_level[p]];
			}
			char out_file[50];
			sprintf(out_file, ".\\output\\%d.bmp", current_iteration);
			cv_imwrite(out_file, dst);
		}
		if (enable_early_stop) {
			if (memcmp(image_level, image_last, sizeof(unsigned char) * pixel_count) == 0) {
				early_stop_counter++;
				printf("Result unchanged for %d iterations.\n", early_stop_counter);
				if (early_stop_counter >= enable_early_stop) {
					printf("Early stop.\n");
					for (int p = 0; p < pixel_count; p++) {
						image_dst[p] = pixel_level[image_level[p]];
					}
					break;
				}
			} else {
				early_stop_counter = 0;
			}
		}
	}

	cv_imwrite(dst_path, dst);

	free(image_in);
	free(image_dst);
	free(image_level);
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
