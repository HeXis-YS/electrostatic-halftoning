#include "cv.hpp"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int ElectrostaticHalftoning2010(struct CMat src, struct CMat *dst, int enable_initial_charge, int max_iterations, int enable_gridforce, int enable_shake, int enable_debug) {

	//////////////////////////////////////////////////////////////////////////
	///// exceptions
	// For backward compatibility
	int error = 0;
	// if (src.type() != CV_8U) {
	// 	CV_Error(CV_BadNumChannels, "[pixkit::halftoning::ElectrostaticHalftoning] image should be grayscale");
	// }
	if (enable_initial_charge != 0 && enable_initial_charge != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] enable_initial_charge should be 0 or 1");
		error = 1;
	} else if (max_iterations < 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] max_iterations should be bigger than 1");
		error = 2;
	} else if (enable_gridforce != 0 && enable_gridforce != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] enable_gridforce should be 0 or 1");
		error = 3;
	} else if (enable_shake != 0 && enable_shake != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] enable_shake should be 0 or 1");
		error = 4;
	} else if (enable_shake == 1 && max_iterations <= 64) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] max_iterations should be bigger than 64");
		error = 5;
	} else if (enable_debug != 0 && enable_debug != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] enable_debug should be 0 or 1");
		error = 6;
	}
	if (error) {
		return error;
	}

	int rows = src.rows;
	int cols = src.cols;
	int pixel_count = rows * cols;
	double *image_in = (double *)malloc(sizeof(double) * pixel_count);
	dst->rows = rows;
	dst->cols = cols;
	dst->data = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);

	//////////////////////////////////////////////////////////////////////////
	///// Initialization
	int particle_count = 0;
	memset(dst->data, 255, sizeof(unsigned char) * pixel_count);
	for (int p = 0; p < pixel_count; p++) {
		int tmp = 255 - src.data[p];
		image_in[p] = (double)tmp / 255.0;
		particle_count += tmp;
	}
	particle_count = (particle_count + 127) / 255;
	printf("The number of black pixel(charge) = %d\n", particle_count);

	//////////////////////////////////////////////////////////////////////////
	///// Initialize the Particle's position
	double *position_Y = (double *)malloc(sizeof(double) * particle_count);
	double *position_X = (double *)malloc(sizeof(double) * particle_count);
	for (int particle = 0; particle < particle_count;) {
		int rand_Y = rand() % rows;
		int rand_X = rand() % cols;
		int p = rand_Y * cols + rand_X;
		if (dst->data[p] == 0 || (enable_initial_charge && rand() % 256 <= src.data[p])) {
			continue;
		}
		dst->data[p] = 0;
		position_Y[particle] = (double)rand_Y + 0.5;
		position_X[particle] = (double)rand_X + 0.5;
		particle++;
	}
	if (enable_debug) {
		cv_imwrite(".\\output\\0.bmp", *dst);
	}

	//////////////////////////////////////////////////////////////////////////
	///// process
	double *position_Y_tmp = (double *)malloc(sizeof(double) * particle_count);
	double *position_X_tmp = (double *)malloc(sizeof(double) * particle_count);
	double *distance_X_array = (double *)malloc(sizeof(double) * cols);
	double *distance_X_2_array = (double *)malloc(sizeof(double) * cols);
	double shake_tmp = log10((double)max_iterations) / log10(1024.0) - 0.6;
	double shake_tmp1 = 0.0;
	for (int current_iteration = 1; current_iteration <= max_iterations; current_iteration++) {
		printf("Iteration %d\n", current_iteration);
		memset(dst->data, 255, sizeof(unsigned char) * pixel_count);
		memcpy(position_Y_tmp, position_Y, sizeof(double) * particle_count);
		memcpy(position_X_tmp, position_X, sizeof(double) * particle_count);
		if (enable_shake == 1 && max_iterations > 64 && current_iteration % 10 == 0) {
			shake_tmp1 = shake_tmp * exp(current_iteration / 1000.0);
		}
		for (int current_particle = 0; current_particle < particle_count; current_particle++) {
			double position_X_offset = 0.0;
			double position_Y_offset = 0.0;
			double position_Y_current = position_Y_tmp[current_particle];
			double position_X_current = position_X_tmp[current_particle];

			// Attraction
			for (int x = 0; x < cols; x++) {
				double distance_X = x + 0.5 - position_X_current;
				distance_X_array[x] = distance_X;
				distance_X_2_array[x] = distance_X * distance_X;
			}
			for (int y = 0, p = 0; y < rows; y++) {
				double distance_Y = y + 0.5 - position_Y_current;
				double distance_Y_2 = distance_Y * distance_Y;
				for (int x = 0; x < cols; x++, p++) {
					double image_in_tmp = image_in[p];
					if (image_in[p] == 0.0) {
						continue;
					}
					double tmp = distance_Y_2 + distance_X_2_array[x];
					if (tmp != 0.0) {
						tmp = image_in[p] / tmp;
						position_X_offset += distance_Y * tmp;
						position_Y_offset += distance_X_array[x] * tmp;
					}
				}
			}

			// Repulsion
			for (int particle = 0; particle < particle_count; particle++) {
				if (current_particle != particle) {
					double distance_Y = position_Y_tmp[particle] - position_Y_current;
					double distance_X = position_X_tmp[particle] - position_X_current;
					double tmp = 0.0;
					if (distance_Y != 0.0) {
						tmp += distance_Y * distance_Y;
					}
					if (distance_X != 0.0) {
						tmp += distance_X * distance_X;
					}
					if (tmp != 0.0) {
						tmp = 1.0 / tmp;
						if (distance_Y != 0.0) {
							position_X_offset -= distance_Y * tmp;
						}
						if (distance_X != 0.0) {
							position_Y_offset -= distance_X * tmp;
						}
					}
				}
			}

			// Add GridForce to find discrete particle locations
			if (enable_gridforce) {
				double grid_distance_Y = position_Y_current - (int)position_Y_current;
				double grid_distance_X = position_X_current - (int)position_X_current;
				double tmp = 0.0;
				if (grid_distance_Y != 0.0) {
					grid_distance_Y = (grid_distance_Y < 0.5) ? -grid_distance_Y : 1 - grid_distance_Y;
					tmp += grid_distance_Y * grid_distance_Y;
				}
				if (grid_distance_X != 0.0) {
					grid_distance_X = (grid_distance_X < 0.5) ? -grid_distance_X : 1 - grid_distance_X;
					tmp += grid_distance_X * grid_distance_X;
				}
				if (tmp != 0.0) {
					tmp = sqrt(tmp);
					tmp = 3.5 / (tmp + 10000.0 * pow(tmp, 9.0));
					if (grid_distance_Y != 0.0) {
						position_X_offset += grid_distance_Y * tmp;
					}
					if (grid_distance_X != 0.0) {
						position_X_offset += grid_distance_X * tmp;
					}
				}
			}

			position_X_offset *= 0.1;
			position_Y_offset *= 0.1;

			// Shake
			if (enable_shake == 1 && max_iterations > 64 && current_iteration % 10 == 0) {
				position_X_offset += shake_tmp1;
				position_Y_offset += shake_tmp1;
			}

			// Result (new position of particles)
			position_Y_current += position_X_offset;
			position_X_current += position_Y_offset;
			position_Y_current -= floor(position_Y_current / (double)rows) * (double)rows;
			position_X_current -= floor(position_X_current / (double)cols) * (double)cols;

			// Output
			dst->data[(int)position_Y_current * cols + (int)position_X_current] = 0;

			position_Y[current_particle] = position_Y_current;
			position_X[current_particle] = position_X_current;
		}

		if (enable_debug) {
			char out_file[50];
			sprintf(out_file, ".\\output\\%d.bmp", current_iteration);
			cv_imwrite(out_file, *dst);
		}
	}

	free(image_in);
	free(position_Y);
	free(position_X);
	free(position_Y_tmp);
	free(position_X_tmp);
	free(distance_X_array);
	free(distance_X_2_array);

	return 0;
}
