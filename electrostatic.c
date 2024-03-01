#include "electrostatic.h"
#include "cv.hpp"
#include <stdio.h>
#include <stdlib.h>

int ElectrostaticHalftoning2010(struct CMat src, struct CMat *dst, int InitialCharge, int Iterations, int GridForce, int Shake, int Debug) {

	//////////////////////////////////////////////////////////////////////////
	///// exceptions
	// if(src.type()!=CV_8U){
	// 	CV_Error(CV_BadNumChannels,"[pixkit::halftoning::ElectrostaticHalftoning] image should be grayscale");
	// }
	if (InitialCharge != 0 && InitialCharge != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] InitialCharge should be 0 or 1");
		system("pause");
		exit(0);
	}
	if (Iterations < 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] Iterations should be bigger than 1");
		system("pause");
		exit(0);
	}
	if (GridForce != 0 && GridForce != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] GridForce should be 0 or 1");
		system("pause");
		exit(0);
	}
	if (Shake != 0 && Shake != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] Shake should be 0 or 1");
		system("pause");
		exit(0);
	}
	if (Shake == 1 && Iterations <= 64) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] Iterations should be bigger than 64");
		system("pause");
		exit(0);
	}
	if (Debug != 0 && Debug != 1 && Debug != 2) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] Debug should be 0, 1 or 2");
		system("pause");
		exit(0);
	}

	char out_file[50];
	double **image_in = (double **)malloc(sizeof(double *) * src.rows);
	unsigned char **image_tmp = (unsigned char **)malloc(sizeof(unsigned char *) * src.rows);
	for (int i = 0; i < src.rows; i++) {
		image_in[i] = (double *)malloc(sizeof(double) * src.cols);
		image_tmp[i] = (unsigned char *)malloc(sizeof(unsigned char) * src.cols);
	}
	dst->rows = src.rows;
	dst->cols = src.cols;
	dst->data = (unsigned char *)malloc(sizeof(unsigned char) * src.rows * src.cols);

	//////////////////////////////////////////////////////////////////////////
	///// Initialization
	for (int i = 0; i < src.rows; i++) {
		for (int j = 0; j < src.cols; j++) {
			image_in[i][j] = (double)src.data[i * src.cols + j] / 255;
			image_tmp[i][j] = 255;
			dst->data[i * src.cols + j] = 255;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	///// Find the number of Particle
	double CountParticle = 0;
	for (int i = 0; i < src.rows; i++) {
		for (int j = 0; j < src.cols; j++) {
			CountParticle = CountParticle + (1 - image_in[i][j]);
		}
	}
	printf("The number of black pixel(charge) = %d\n", (int)CountParticle);

	//////////////////////////////////////////////////////////////////////////
	///// Initialize the Particle's position
	double *Particle_Y = (double *)malloc(sizeof(double) * (int)CountParticle);
	double *Particle_X = (double *)malloc(sizeof(double) * (int)CountParticle);
	int Particle = CountParticle;
	while (Particle > 0) {
		int RandY = rand() % src.rows;
		int RandX = rand() % src.cols;
		if (image_tmp[RandY][RandX] == 0) {
			continue;
		}
		// int RandNumber = rand() % 256;
		if (InitialCharge && rand() % 256 <= src.data[RandY * src.cols + RandX]) {
			continue;
		}
		image_tmp[RandY][RandX] = 0;
		if (Debug == 1) {
			dst->data[RandY * src.cols + RandX] = 0;
		}
		Particle--;
	}
	if (Debug == 1) {
		cv_imwrite("output.bmp", *dst);
	} else if (Debug == 2) {
		sprintf(out_file, ".\\output\\0.bmp");
		cv_imwrite(out_file, *dst);
	}

	//////////////////////////////////////////////////////////////////////////
	///// Record the Particle's position
	int ParticleNumber = 0;
	for (int i = 0; i < src.rows; i++) {
		for (int j = 0; j < src.cols; j++) {
			if (image_tmp[i][j] == 0) {
				Particle_Y[ParticleNumber] = (double)i;
				Particle_X[ParticleNumber] = (double)j;
				ParticleNumber++;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	///// Create Forcefield Table
	printf("Create Forcefield Table, \n");
	double **forcefield_y = (double **)malloc(sizeof(double *) * src.rows);
	double **forcefield_x = (double **)malloc(sizeof(double *) * src.rows);
	for (int i = 0; i < src.rows; i++) {
		forcefield_y[i] = (double *)malloc(sizeof(double *) * src.cols);
		forcefield_x[i] = (double *)malloc(sizeof(double *) * src.cols);
	}
	for (int i = 0; i < src.rows; i++) {
		for (int j = 0; j < src.cols; j++) {
			forcefield_y[i][j] = 0;
			forcefield_x[i][j] = 0;
			int a = -i;
			for (int y = 0; y < src.rows; y++) {
				int aa = a * a;
				int b = -j;
				for (int x = 0; x < src.cols; x++) {
					if (!(i == y && j == x)) {
						double t = (1 - image_in[y][x]) / (aa + b * b);
						// forcefield_y[i][j] += (1 - image_in[y][x]) * (y - i) / ((y - i) * (y - i) + (x - j) * (x - j));
						// forcefield_x[i][j] += (1 - image_in[y][x]) * (x - j) / ((y - i) * (y - i) + (x - j) * (x - j));
						forcefield_y[i][j] += a * t;
						forcefield_x[i][j] += b * t;
					}
					b++;
				}
				a++;
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

			// Attraction(by using bilinear interpolation)
			if (Particle_Y[NowCharge] - (int)Particle_Y[NowCharge] == 0 && Particle_X[NowCharge] - (int)Particle_X[NowCharge] == 0) {
				NewPosition_Y = forcefield_y[(int)Particle_Y[NowCharge]][(int)Particle_X[NowCharge]];
				NewPosition_X = forcefield_x[(int)Particle_Y[NowCharge]][(int)Particle_X[NowCharge]];
			} else {
				int Bilinear_y1 = Particle_Y[NowCharge];
				int Bilinear_x1 = Particle_X[NowCharge];
				int Bilinear_y2 = Bilinear_y1 + 1;
				int Bilinear_x2 = Bilinear_x1 + 1;
				if (Bilinear_y2 < src.rows && Bilinear_x2 < src.cols) {
					double a = (double)Bilinear_x2 - Particle_X[NowCharge];
					double b = (double)Bilinear_y2 - Particle_Y[NowCharge];
					double c = Particle_X[NowCharge] - (double)Bilinear_x1;
					double d = Particle_Y[NowCharge] - (double)Bilinear_y1;
					NewPosition_Y = forcefield_y[Bilinear_y1][Bilinear_x1] * a * b + forcefield_y[Bilinear_y1][Bilinear_x2] * c * b + forcefield_y[Bilinear_y2][Bilinear_x1] * a * d + forcefield_y[Bilinear_y2][Bilinear_x2] * c * d;
					NewPosition_X = forcefield_x[Bilinear_y1][Bilinear_x1] * a * b + forcefield_x[Bilinear_y1][Bilinear_x2] * c * b + forcefield_x[Bilinear_y2][Bilinear_x1] * a * d + forcefield_x[Bilinear_y2][Bilinear_x2] * c * d;
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
					real_y = -real_y;
				} else {
					real_y = 1 - real_y;
				}
				if (real_x < 0.5) {
					real_x = -real_x;
				} else {
					real_x = 1 - real_x;
				}
				double vector3 = sqrt(real_y * real_y + real_x * real_x);
				if (real_y != 0) {
					GridForce_Y = 3.5 * real_y / (vector3 * (1 + pow(vector3, 8) * 10000));
				}
				if (real_x != 0) {
					GridForce_X = 3.5 * real_x / (vector3 * (1 + pow(vector3, 8) * 10000));
				}
			}

			// result (new position of particles)
			if (GridForce) {
				Particle_Y[NowCharge] = Particle_Y[NowCharge] + 0.1 * (NewPosition_Y + GridForce_Y);
				Particle_X[NowCharge] = Particle_X[NowCharge] + 0.1 * (NewPosition_X + GridForce_X);
			} else {
				Particle_Y[NowCharge] = Particle_Y[NowCharge] + 0.1 * NewPosition_Y;
				Particle_X[NowCharge] = Particle_X[NowCharge] + 0.1 * NewPosition_X;
			}

			// Shake
			if (Shake == 1 && iterations % 10 == 0 && Iterations > 64) {
				Particle_Y[NowCharge] += (log10((double)Iterations) / log10(2.0) - 6) * exp(iterations / 1000.0) / 10;
				Particle_X[NowCharge] += (log10((double)Iterations) / log10(2.0) - 6) * exp(iterations / 1000.0) / 10;
			}

			if (Particle_Y[NowCharge] < 0) {
				Particle_Y[NowCharge] = 0;
			} else if (Particle_Y[NowCharge] >= src.rows) {
				Particle_Y[NowCharge] = src.rows - 1;
			}
			if (Particle_X[NowCharge] < 0) {
				Particle_X[NowCharge] = 0;
			} else if (Particle_X[NowCharge] >= src.cols) {
				Particle_X[NowCharge] = src.cols - 1;
			}
		}

		// Output
		for (int y = 0; y < src.rows; y++) {
			for (int x = 0; x < src.cols; x++) {
				dst->data[y * src.cols + x] = 255;
				image_tmp[y][x] = 255;
			}
		}
		int output_position;
		int out_Y, out_X;
		double count_errorY = 0, count_errorX = 0;
		for (int NowCharge = 0; NowCharge < Particle; NowCharge++) {
			out_Y = Particle_Y[NowCharge] + 0.5;
			out_X = Particle_X[NowCharge] + 0.5;
			if (out_Y >= src.rows) {
				out_Y = src.rows - 1;
			}
			if (out_X >= src.cols) {
				out_X = src.cols - 1;
			}
			image_tmp[out_Y][out_X] = 0;
		}

		for (int y = 0; y < src.rows; y++) {
			for (int x = 0; x < src.cols; x++) {
				dst->data[y * src.cols + x] = image_tmp[y][x];
			}
		}

		if (Debug == 1) {
			cv_imwrite("output.bmp", *dst);
		} else if (Debug == 2) {
			sprintf(out_file, ".\\output\\%d.bmp", iterations);
			cv_imwrite(out_file, *dst);
		}
	}

	// dst = dst->clone();

	free(image_in);
	free(image_tmp);

	return 1;
}
