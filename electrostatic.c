#include "cv.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int ElectrostaticHalftoning2010(struct CMat src, struct CMat *dst, int InitialCharge, int Iterations, int GridForce, int Shake, int Debug) {

	//////////////////////////////////////////////////////////////////////////
	///// exceptions
	// For backward compatibility
	int error = 0;
	// if (src.type() != CV_8U) {
	// 	CV_Error(CV_BadNumChannels, "[pixkit::halftoning::ElectrostaticHalftoning] image should be grayscale");
	// }
	if (InitialCharge != 0 && InitialCharge != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] InitialCharge should be 0 or 1");
		error = 1;
	} else if (Iterations < 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] Iterations should be bigger than 1");
		error = 2;
	} else if (GridForce != 0 && GridForce != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] GridForce should be 0 or 1");
		error = 3;
	} else if (Shake != 0 && Shake != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] Shake should be 0 or 1");
		error = 4;
	} else if (Shake == 1 && Iterations <= 64) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] Iterations should be bigger than 64");
		error = 5;
	} else if (Debug != 0 && Debug != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] Debug should be 0 or 1");
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
	double *Particle_Y = (double *)malloc(sizeof(double) * particle_count);
	double *Particle_X = (double *)malloc(sizeof(double) * particle_count);
	for (int particle = 0; particle < particle_count;) {
		int RandY = rand() % rows;
		int RandX = rand() % cols;
		int p = RandY * cols + RandX;
		if (dst->data[p] == 0 || (InitialCharge && rand() % 256 <= src.data[p])) {
			continue;
		}
		dst->data[p] = 0;
		Particle_Y[particle] = (double)RandY + 0.5;
		Particle_X[particle] = (double)RandX + 0.5;
		particle++;
	}
	if (Debug) {
		cv_imwrite(".\\output\\0.bmp", *dst);
	}

	//////////////////////////////////////////////////////////////////////////
	///// process
	int Particle = particle_count;
	double *distance_X_array = (double *)malloc(sizeof(double) * cols);
	double *distance_X_2_array = (double *)malloc(sizeof(double) * cols);
	double shake_tmp = log10((double)Iterations) / log10(1024.0) - 0.6;
	double shake_tmp1 = 0.0;
	for (int iterations = 1; iterations <= Iterations; iterations++) {
		printf("Iterations %d\n", iterations);
		memset(dst->data, 255, sizeof(unsigned char) * pixel_count);
		if (Shake == 1 && Iterations > 64 && iterations % 10 == 0) {
			shake_tmp1 = shake_tmp * exp(iterations / 1000.0);
		}
		for (int NowCharge = 0; NowCharge < Particle; NowCharge++) {
			double NewPosition_Y = 0, NewPosition_X = 0;

			// Attraction
			for (int x = 0; x < cols; x++) {
				double distance_X = x + 0.5 - Particle_X[NowCharge];
				distance_X_array[x] = distance_X;
				distance_X_2_array[x] = distance_X * distance_X;
			}
			for (int y = 0, p = 0; y < rows; y++) {
				double distance_Y = y + 0.5 - Particle_Y[NowCharge];
				double distance_Y_2 = distance_Y * distance_Y;
				for (int x = 0; x < cols; x++, p++) {
					if (image_in[p] == 0.0) {
						continue;
					}
					double tmp = distance_Y_2 + distance_X_2_array[x];
					if (tmp != 0) {
						tmp = image_in[p] / tmp;
						NewPosition_Y += distance_Y * tmp;
						NewPosition_X += distance_X_array[x] * tmp;
					}
				}
			}

			// Repulsion
			for (int OtherCharge = 0; OtherCharge < Particle; OtherCharge++) {
				if (NowCharge != OtherCharge) {
					double distance_Y = Particle_Y[OtherCharge] - Particle_Y[NowCharge];
					double distance_X = Particle_X[OtherCharge] - Particle_X[NowCharge];
					double tmp = 0.0;
					if (distance_Y != 0.0) {
						tmp += distance_Y * distance_Y;
					}
					if (distance_X != 0.0) {
						tmp += distance_X * distance_X;
					}
					if (tmp != 0) {
						tmp = 1.0 / tmp;
						if (distance_Y != 0.0) {
							NewPosition_Y -= distance_Y * tmp;
						}
						if (distance_X != 0.0) {
							NewPosition_X -= distance_X * tmp;
						}
					}
				}
			}

			// Add GridForce to find discrete particle locations
			if (GridForce) {
				double grid_distance_Y = Particle_Y[NowCharge] - (int)Particle_Y[NowCharge];
				double grid_distance_X = Particle_X[NowCharge] - (int)Particle_X[NowCharge];
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
					tmp = 3.5 / (tmp + 10000 * pow(tmp, 9));
					if (grid_distance_Y != 0) {
						NewPosition_Y += grid_distance_Y * tmp;
					}
					if (grid_distance_X != 0) {
						NewPosition_Y += grid_distance_X * tmp;
					}
				}
			}

			// Result (new position of particles)
			Particle_Y[NowCharge] += 0.1 * NewPosition_Y;
			Particle_X[NowCharge] += 0.1 * NewPosition_X;

			// Shake
			if (Shake == 1 && Iterations > 64 && iterations % 10 == 0) {
				Particle_Y[NowCharge] += shake_tmp1;
				Particle_X[NowCharge] += shake_tmp1;
			}

			Particle_Y[NowCharge] = Particle_Y[NowCharge] - floor(Particle_Y[NowCharge] / (double)rows) * (double)rows;
			Particle_X[NowCharge] = Particle_X[NowCharge] - floor(Particle_X[NowCharge] / (double)cols) * (double)cols;

			// Output
			dst->data[(int)Particle_Y[NowCharge] * cols + (int)Particle_X[NowCharge]] = 0;
		}

		if (Debug) {
			char out_file[50];
			sprintf(out_file, ".\\output\\%d.bmp", iterations);
			cv_imwrite(out_file, *dst);
		}
	}

	// dst = dst->clone();

	free(image_in);

	return 0;
}
