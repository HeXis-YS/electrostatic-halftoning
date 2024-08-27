#include "electrostatic.h"
#include "cv.hpp"
#include <stdio.h>
#include <stdlib.h>

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
	for (int p = 0; p < pixel_count; p++) {
		image_in[p] = (double)src.data[p] / 255;
		dst->data[p] = 255;
	}

	//////////////////////////////////////////////////////////////////////////
	///// Find the number of Particle
	double CountParticle = 0;
	for (int p = 0; p < pixel_count; p++) {
		CountParticle = CountParticle + (1 - image_in[p]);
	}
	printf("The number of black pixel(charge) = %d\n", (int)CountParticle);

	//////////////////////////////////////////////////////////////////////////
	///// Initialize the Particle's position
	double *Particle_Y = (double *)malloc(sizeof(double) * (int)CountParticle);
	double *Particle_X = (double *)malloc(sizeof(double) * (int)CountParticle);
	int Particle = CountParticle;
	while (Particle > 0) {
		int RandY = rand() % rows;
		int RandX = rand() % cols;
		int p = RandY * cols + RandX;
		if (dst->data[p] != 0) {
			if (InitialCharge == 1) {
				int RandNumber = rand() % 256;
				if (RandNumber > src.data[p]) {
					dst->data[p] = 0;
					Particle--;
				}
			} else if (InitialCharge == 0) {
				dst->data[p] = 0;
				Particle--;
			}
		}
	}
	if (Debug) {
		cv_imwrite(".\\output\\0.bmp", *dst);
	}

	//////////////////////////////////////////////////////////////////////////
	///// Record the Particle's position
	int ParticleNumber = 0;
	for (int i = 0, p = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++, p++) {
			if (dst->data[p] == 0) {
				Particle_Y[ParticleNumber] = (double)i;
				Particle_X[ParticleNumber] = (double)j;
				ParticleNumber++;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	///// process
	double instead_y, instead_x;
	Particle = CountParticle;
	for (int iterations = 1; iterations <= Iterations; iterations++) {
		printf("Iterations %d\n", iterations);

		for (int NowCharge = 0; NowCharge < Particle; NowCharge++) {
			double NewPosition_Y = 0, NewPosition_X = 0;
			double GridForce_Y = 0, GridForce_X = 0;

			// Attraction
			double i = Particle_Y[NowCharge] - 0.5;
			double j = Particle_X[NowCharge] - 0.5;
			for (int y = 0, p = 0; y < rows; y++) {
				for (int x = 0; x < cols; x++, p++) {
					NewPosition_Y += (1 - image_in[p]) * (y - i) / ((y - i) * (y - i) + (x - j) * (x - j));
					NewPosition_X += (1 - image_in[p]) * (x - j) / ((y - i) * (y - i) + (x - j) * (x - j));
				}
			}

			// Repulsion
			for (int OtherCharge = 0; OtherCharge < Particle; OtherCharge++) {
				if (NowCharge != OtherCharge) {
					instead_y = Particle_Y[OtherCharge] - Particle_Y[NowCharge];
					instead_x = Particle_X[OtherCharge] - Particle_X[NowCharge];
					if (!(instead_y == 0 && instead_x == 0)) {
						NewPosition_Y -= instead_y / (instead_y * instead_y + instead_x * instead_x);
						NewPosition_X -= instead_x / (instead_y * instead_y + instead_x * instead_x);
					}
				}
			}

			// Add GridForce to find discrete particle locations
			double real_y = Particle_Y[NowCharge] - (int)Particle_Y[NowCharge];
			double real_x = Particle_X[NowCharge] - (int)Particle_X[NowCharge];
			if (real_y == 0 && real_x == 0) {
				GridForce_Y = 0;
				GridForce_X = 0;
			} else {
				if (real_y < 0.5) {
					if (real_x < 0.5) {
						real_y = (0 - real_y);
						real_x = (0 - real_x);
					} else {
						real_y = (0 - real_y);
						real_x = (1 - real_x);
					}
				} else {
					if (real_x < 0.5) {
						real_y = (1 - real_y);
						real_x = (0 - real_x);
					} else {
						real_y = (1 - real_y);
						real_x = (1 - real_x);
					}
				}
				double vector3 = sqrt(real_y * real_y + real_x * real_x);
				if (real_y == 0) {
					GridForce_Y = 0;
				} else {
					GridForce_Y = 3.5 * real_y / (vector3 * (1 + pow(vector3, 8) * 10000));
				}
				if (real_x == 0) {
					GridForce_X = 0;
				} else {
					GridForce_X = 3.5 * real_x / (vector3 * (1 + pow(vector3, 8) * 10000));
				}
			}

			// resault (new position of particles)
			if (GridForce == 0) {
				Particle_Y[NowCharge] = Particle_Y[NowCharge] + 0.1 * NewPosition_Y;
				Particle_X[NowCharge] = Particle_X[NowCharge] + 0.1 * NewPosition_X;
			} else if (GridForce == 1) {
				Particle_Y[NowCharge] = Particle_Y[NowCharge] + 0.1 * (NewPosition_Y + GridForce_Y);
				Particle_X[NowCharge] = Particle_X[NowCharge] + 0.1 * (NewPosition_X + GridForce_X);
			}

			// Shake
			if (Shake == 1 && iterations % 10 == 0 && Iterations > 64) {
				Particle_Y[NowCharge] += (log10((double)Iterations) / log10(2.0) - 6) * exp(iterations / 1000.0) / 10;
				Particle_X[NowCharge] += (log10((double)Iterations) / log10(2.0) - 6) * exp(iterations / 1000.0) / 10;
			}

			if (Particle_Y[NowCharge] < 0) {
				Particle_Y[NowCharge] = 0;
			}
			if (Particle_Y[NowCharge] >= rows) {
				Particle_Y[NowCharge] = rows - 1;
			}
			if (Particle_X[NowCharge] < 0) {
				Particle_X[NowCharge] = 0;
			}
			if (Particle_X[NowCharge] >= cols) {
				Particle_X[NowCharge] = cols - 1;
			}
		}

		// Output
		for (int p = 0; p < pixel_count; p++) {
			dst->data[p] = 255;
		}
		int output_position;
		int out_Y, out_X;
		double count_errorY = 0, count_errorX = 0;
		for (int NowCharge = 0; NowCharge < Particle; NowCharge++) {
			out_Y = Particle_Y[NowCharge] + 0.5;
			out_X = Particle_X[NowCharge] + 0.5;
			if (out_Y >= rows) {
				out_Y = rows - 1;
			}
			if (out_X >= cols) {
				out_X = cols - 1;
			}
			dst->data[out_Y * cols + out_X] = 0;
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
