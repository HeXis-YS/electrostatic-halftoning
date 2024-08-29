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
								int enable_adaptive_learning_rate,
								int enable_early_stop,
								int enable_debug) {
	//////////////////////////////////////////////////////////////////////////
	///// exceptions
	max_iterations = (max_iterations > 0) ? max_iterations : 8;
	enable_initial_charge = enable_initial_charge ? 1 : 0;
	enable_gridforce = enable_gridforce ? 1 : 0;
	enable_shake = enable_shake ? 1 : 0;
	enable_adaptive_learning_rate = enable_adaptive_learning_rate ? 1 : 0;
	enable_early_stop = enable_early_stop ? enable_early_stop : 0;
	enable_debug = enable_debug ? 1 : 0;
	printf("Max iterations = %d\n", max_iterations);
	printf("Initial charge = %s\n", enable_initial_charge ? "Enabled" : "Disabled");
	printf("Grid Force = %s\n", enable_gridforce ? "Enabled" : "Disabled");
	printf("Shake = %s\n", enable_shake ? "Enabled" : "Disabled");
	printf("Adaptive learning rate = %s\n", enable_adaptive_learning_rate ? "Enabled" : "Disabled");
	printf("Early stop = ");
	enable_early_stop ? printf("%d iterations.\n", enable_early_stop) : printf("Disabled.\n");
	if (enable_shake) {
		if (max_iterations <= 64) {
			printf("Error: max_iterations > 64, when enable_shake = 1\n");
			return 1;
		} else if (enable_adaptive_learning_rate) {
			printf("Error: max_iterations != 1, when enable_shake = 1\n");
			return 2;
		}
	}
	if (enable_adaptive_learning_rate && enable_gridforce) {
		printf("Warning: Not recommended for use with grid force with adaptive learning rate.\n");
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
	int particle_count = 0;
	memset(image_dst, 255, sizeof(unsigned char) * pixel_count);
	for (int p = 0; p < pixel_count; p++) {
		int tmp = 255 - src.data[p];
		image_in[p] = (double)tmp / 255.0;
		particle_count += tmp;
	}
	particle_count = (particle_count + 127) / 255;
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
		image_dst[p] = 0;
		particle_Y[particle] = (double)rand_Y + rand_double();
		particle_X[particle] = (double)rand_X + rand_double();
		particle++;
	}
	if (enable_debug) {
		cv_imwrite(".\\output\\0.bmp", *dst);
	}

	//////////////////////////////////////////////////////////////////////////
	///// process
	// Early stop
	unsigned char *image_last = NULL;
	int early_stop_counter = 0;
	if (enable_early_stop) {
		image_last = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);
	}
	double *particle_Y_last = (double *)malloc(sizeof(double) * particle_count);
	double *particle_X_last = (double *)malloc(sizeof(double) * particle_count);
	double *distance_X = (double *)malloc(sizeof(double) * cols);
	double *distance_X_2 = (double *)malloc(sizeof(double) * cols);
	double *particle_Y_offset_array = (double *)malloc(sizeof(double) * particle_count);
	double *particle_X_offset_array = (double *)malloc(sizeof(double) * particle_count);
	int shake = 0;
	double shake_tmp = log10((double)max_iterations) / log10(1024.0) - 0.6;
	double shake_tmp1 = 0.0;
	double learning_rate = 0.1;
	double mean = 0.0;
	double mean_last = DBL_MAX;
	if (enable_debug) {
		printf("Initial learning rate = %f\n", learning_rate);
	}
	for (int current_iteration = 1; current_iteration <= max_iterations; current_iteration++) {
		printf("Iteration %d\n", current_iteration);
		if (enable_early_stop) {
			memcpy(image_last, image_dst, sizeof(unsigned char) * pixel_count);
		}
		memset(image_dst, 255, sizeof(unsigned char) * pixel_count);
		memcpy(particle_Y_last, particle_Y, sizeof(double) * particle_count);
		memcpy(particle_X_last, particle_X, sizeof(double) * particle_count);
		if (enable_shake) {
			if (possibility(current_iteration % 10 == 0, 0.1)) {
				shake_tmp1 = shake_tmp * exp(current_iteration / 1000.0);
				shake = 1;
			} else {
				shake = 0;
			}
		}
		if (enable_adaptive_learning_rate) {
			mean = 0.0;
		}
		for (int current_particle = 0; current_particle < particle_count; current_particle++) {
			double particle_Y_offset = 0.0;
			double particle_X_offset = 0.0;
			double particle_Y_current = particle_Y_last[current_particle];
			double particle_X_current = particle_X_last[current_particle];

			// Attraction
			for (int x = 0; x < cols; x++) {
				double tmp = x + 0.5 - particle_X_current;
				distance_X[x] = tmp;
				distance_X_2[x] = tmp * tmp;
			}
			for (int y = 0, p = 0; y < rows; y++) {
				double distance_Y = y + 0.5 - particle_Y_current;
				double distance_Y_2 = distance_Y * distance_Y;
				for (int x = 0; x < cols; x++, p++) {
					double image_in_tmp = image_in[p];
					if (unlikely(image_in[p] == 0.0)) {
						continue;
					}
					double tmp = distance_Y_2 + distance_X_2[x];
					if (likely(tmp != 0.0)) {
						tmp = image_in[p] / tmp;
						particle_Y_offset += distance_Y * tmp;
						particle_X_offset += distance_X[x] * tmp;
					}
				}
			}

			// Repulsion
			for (int particle = 0; particle < particle_count; particle++) {
				if (unlikely(current_particle == particle)) {
					continue;
				}
				double distance_Y = particle_Y_last[particle] - particle_Y_current;
				double distance_X = particle_X_last[particle] - particle_X_current;
				double tmp = 0.0;
				if (likely(distance_Y != 0.0)) {
					tmp += distance_Y * distance_Y;
				}
				if (likely(distance_X != 0.0)) {
					tmp += distance_X * distance_X;
				}
				if (likely(tmp != 0.0)) {
					tmp = 1.0 / tmp;
					if (likely(distance_Y != 0.0)) {
						particle_Y_offset -= distance_Y * tmp;
					}
					if (likely(distance_X != 0.0)) {
						particle_X_offset -= distance_X * tmp;
					}
				}
			}

			// Add GridForce to find discrete particle locations
			if (enable_gridforce) {
				double grid_distance_Y = particle_Y_current - (int)particle_Y_current;
				double grid_distance_X = particle_X_current - (int)particle_X_current;
				double tmp = 0.0;
				if (likely(grid_distance_Y != 0.0)) {
					grid_distance_Y = possibility(grid_distance_Y <= 0.5, 0.5) ? -grid_distance_Y : 1 - grid_distance_Y;
					tmp += grid_distance_Y * grid_distance_Y;
				}
				if (likely(grid_distance_X != 0.0)) {
					grid_distance_X = possibility(grid_distance_X <= 0.5, 0.5) ? -grid_distance_X : 1 - grid_distance_X;
					tmp += grid_distance_X * grid_distance_X;
				}
				if (likely(tmp != 0.0)) {
					tmp = sqrt(tmp);
					tmp = 3.5 / (tmp + 10000.0 * pow(tmp, 9.0));
					if (likely(grid_distance_Y != 0.0)) {
						particle_Y_offset += grid_distance_Y * tmp;
					}
					if (likely(grid_distance_X != 0.0)) {
						particle_X_offset += grid_distance_X * tmp;
					}
				}
			}

			particle_Y_offset_array[current_particle] = particle_Y_offset;
			particle_X_offset_array[current_particle] = particle_X_offset;

			if (enable_adaptive_learning_rate) {
				mean += sqrt(particle_Y_offset * particle_Y_offset + particle_X_offset * particle_X_offset);
			}
		}
		if (enable_adaptive_learning_rate) {
			mean /= (double)pixel_count;
			while (mean * learning_rate >= mean_last) {
				learning_rate *= 0.5;
				printf("Learning rate reduced to %f\n", learning_rate);
			}
			printf("Mean force/offset = %f/%f\n", mean, mean * learning_rate);
			mean_last = mean * learning_rate;
		}
		for (int current_particle = 0; current_particle < particle_count; current_particle++) {
			double particle_Y_current = particle_Y_last[current_particle] + particle_Y_offset_array[current_particle] * learning_rate;
			double particle_X_current = particle_X_last[current_particle] + particle_X_offset_array[current_particle] * learning_rate;

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
			image_dst[(int)particle_Y_current * cols + (int)particle_X_current] = 0;
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
	free(particle_Y_offset_array);
	free(particle_X_offset_array);
	free(particle_Y_last);
	free(particle_X_last);
	free(distance_X);
	free(distance_X_2);

	return 0;
}
