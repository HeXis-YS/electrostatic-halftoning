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
	int pixel_count = src.cols * src.rows;
	double *image_in = (double *)malloc(sizeof(double) * pixel_count);
	unsigned char *image_tmp = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);
	dst->cols = src.cols;
	dst->rows = src.rows;
	dst->data = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);

	//////////////////////////////////////////////////////////////////////////
	///// Initialization
	for (int p = 0; p < pixel_count; p++) {
		image_in[p] = (double)src.data[p] / 255;
		image_tmp[p] = 255;
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
		int RandY = rand() % src.rows;
		int RandX = rand() % src.cols;
		int p = RandY * src.cols + RandX;
		if (image_tmp[p] == 0) {
			continue;
		}
		// int RandNumber = rand() % 256;
		if (InitialCharge && rand() % 256 <= src.data[p]) {
			continue;
		}
		image_tmp[p] = 0;
		if (Debug == 1) {
			dst->data[p] = 0;
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
	for (int i = 0, p = 0; i < src.rows; i++) {
		for (int j = 0; j < src.cols; j++) {
			if (image_tmp[p] == 0) {
				Particle_Y[ParticleNumber] = (double)i;
				Particle_X[ParticleNumber] = (double)j;
				ParticleNumber++;
			}
			p++;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	///// Create Forcefield Table
	printf("Create Forcefield Table, \n");
	double *forcefield_y = (double *)malloc(sizeof(double) * pixel_count);
	double *forcefield_x = (double *)malloc(sizeof(double) * pixel_count);
	for (int i = 0, p1 = 0; i < src.rows; i++) {
		for (int j = 0; j < src.cols; j++) {
			forcefield_y[p1] = 0;
			forcefield_x[p1] = 0;
			int a = -i;
			for (int y = 0, p2 = 0; y < src.rows; y++) {
				int aa = a * a;
				int b = -j;
				for (int x = 0; x < src.cols; x++) {
					if (i != y || j != x) {
						double t = (1 - image_in[p2]) / (aa + b * b);
						// forcefield_y[i][j] += (1 - image_in[y][x]) * (y - i) / ((y - i) * (y - i) + (x - j) * (x - j));
						// forcefield_x[i][j] += (1 - image_in[y][x]) * (x - j) / ((y - i) * (y - i) + (x - j) * (x - j));
						forcefield_y[p1] += a * t;
						forcefield_x[p1] += b * t;
					}
					b++;
					p2++;
				}
				a++;
			}
			p1++;
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
			double real_y = Particle_Y[NowCharge] - (int)Particle_Y[NowCharge];
			double real_x = Particle_X[NowCharge] - (int)Particle_X[NowCharge];

			// Attraction(by using bilinear interpolation)
			if (real_y == 0 && real_x == 0) {
				int p = (int)Particle_Y[NowCharge] * src.cols + (int)Particle_X[NowCharge];
				NewPosition_Y = forcefield_y[p];
				NewPosition_X = forcefield_x[p];
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
					int p11 = Bilinear_y1 * src.cols + Bilinear_x1;
					int p12 = Bilinear_y1 * src.cols + Bilinear_x2;
					int p21 = Bilinear_y2 * src.cols + Bilinear_x1;
					int p22 = Bilinear_y2 * src.cols + Bilinear_x2;
					NewPosition_Y = forcefield_y[p11] * a * b + forcefield_y[p12] * c * b + forcefield_y[p21] * a * d + forcefield_y[p22] * c * d;
					NewPosition_X = forcefield_x[p11] * a * b + forcefield_x[p12] * c * b + forcefield_x[p21] * a * d + forcefield_x[p22] * c * d;
				}
			}

			// Repulsion
			for (int OtherCharge = 0; OtherCharge < Particle; OtherCharge++) {
				if (NowCharge != OtherCharge) {
					instead_y = Particle_Y[OtherCharge] - Particle_Y[NowCharge];
					instead_x = Particle_X[OtherCharge] - Particle_X[NowCharge];
					if (instead_y != 0 || instead_x != 0) {
						double t = instead_y * instead_y + instead_x * instead_x;
						NewPosition_Y -= instead_y / t;
						NewPosition_X -= instead_x / t;
					}
				}
			}

			// result (new position of particles)
			Particle_Y[NowCharge] += 0.1 * NewPosition_Y;
			Particle_X[NowCharge] += 0.1 * NewPosition_X;
			if (GridForce) {
				// Add GridForce to find discrete particle locations
				// double real_y = Particle_Y[NowCharge] - (int)Particle_Y[NowCharge];
				// double real_x = Particle_X[NowCharge] - (int)Particle_X[NowCharge];
				if (real_y != 0 || real_x != 0) {
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
					double t = 0.35 / (vector3 + pow(vector3, 9) * 10000);
					if (real_y != 0) {
						// GridForce_Y = 3.5 * real_y / (vector3 * (1 + pow(vector3, 8) * 10000));
						Particle_Y[NowCharge] += real_y * t;
					}
					if (real_x != 0) {
						// GridForce_X = 3.5 * real_x / (vector3 * (1 + pow(vector3, 8) * 10000));
						Particle_X[NowCharge] += real_x * t;
					}
				}
			}

			// Shake
			if (Shake && iterations % 10 == 0 && Iterations > 64) {
				double t = (log10((double)Iterations) / log10(2.0) - 6) * exp(iterations / 1000.0) / 10;
				Particle_Y[NowCharge] += t;
				Particle_X[NowCharge] += t;
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
		for (int p = 0; p < pixel_count; p++) {
			dst->data[p] = 255;
			image_tmp[p] = 255;
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
			image_tmp[out_Y * src.cols + out_X] = 0;
		}

		for (int p = 0; p < pixel_count; p++) {
			dst->data[p] = image_tmp[p];
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
