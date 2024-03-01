#include "electrostatic.h"
#include "cv.hpp"
#include <stdio.h>
#include <stdlib.h>

static const double Pixel_LUT[] = {
	255.0 / 255.0,
	254.0 / 255.0,
	253.0 / 255.0,
	252.0 / 255.0,
	251.0 / 255.0,
	250.0 / 255.0,
	249.0 / 255.0,
	248.0 / 255.0,
	247.0 / 255.0,
	246.0 / 255.0,
	245.0 / 255.0,
	244.0 / 255.0,
	243.0 / 255.0,
	242.0 / 255.0,
	241.0 / 255.0,
	240.0 / 255.0,
	239.0 / 255.0,
	238.0 / 255.0,
	237.0 / 255.0,
	236.0 / 255.0,
	235.0 / 255.0,
	234.0 / 255.0,
	233.0 / 255.0,
	232.0 / 255.0,
	231.0 / 255.0,
	230.0 / 255.0,
	229.0 / 255.0,
	228.0 / 255.0,
	227.0 / 255.0,
	226.0 / 255.0,
	225.0 / 255.0,
	224.0 / 255.0,
	223.0 / 255.0,
	222.0 / 255.0,
	221.0 / 255.0,
	220.0 / 255.0,
	219.0 / 255.0,
	218.0 / 255.0,
	217.0 / 255.0,
	216.0 / 255.0,
	215.0 / 255.0,
	214.0 / 255.0,
	213.0 / 255.0,
	212.0 / 255.0,
	211.0 / 255.0,
	210.0 / 255.0,
	209.0 / 255.0,
	208.0 / 255.0,
	207.0 / 255.0,
	206.0 / 255.0,
	205.0 / 255.0,
	204.0 / 255.0,
	203.0 / 255.0,
	202.0 / 255.0,
	201.0 / 255.0,
	200.0 / 255.0,
	199.0 / 255.0,
	198.0 / 255.0,
	197.0 / 255.0,
	196.0 / 255.0,
	195.0 / 255.0,
	194.0 / 255.0,
	193.0 / 255.0,
	192.0 / 255.0,
	191.0 / 255.0,
	190.0 / 255.0,
	189.0 / 255.0,
	188.0 / 255.0,
	187.0 / 255.0,
	186.0 / 255.0,
	185.0 / 255.0,
	184.0 / 255.0,
	183.0 / 255.0,
	182.0 / 255.0,
	181.0 / 255.0,
	180.0 / 255.0,
	179.0 / 255.0,
	178.0 / 255.0,
	177.0 / 255.0,
	176.0 / 255.0,
	175.0 / 255.0,
	174.0 / 255.0,
	173.0 / 255.0,
	172.0 / 255.0,
	171.0 / 255.0,
	170.0 / 255.0,
	169.0 / 255.0,
	168.0 / 255.0,
	167.0 / 255.0,
	166.0 / 255.0,
	165.0 / 255.0,
	164.0 / 255.0,
	163.0 / 255.0,
	162.0 / 255.0,
	161.0 / 255.0,
	160.0 / 255.0,
	159.0 / 255.0,
	158.0 / 255.0,
	157.0 / 255.0,
	156.0 / 255.0,
	155.0 / 255.0,
	154.0 / 255.0,
	153.0 / 255.0,
	152.0 / 255.0,
	151.0 / 255.0,
	150.0 / 255.0,
	149.0 / 255.0,
	148.0 / 255.0,
	147.0 / 255.0,
	146.0 / 255.0,
	145.0 / 255.0,
	144.0 / 255.0,
	143.0 / 255.0,
	142.0 / 255.0,
	141.0 / 255.0,
	140.0 / 255.0,
	139.0 / 255.0,
	138.0 / 255.0,
	137.0 / 255.0,
	136.0 / 255.0,
	135.0 / 255.0,
	134.0 / 255.0,
	133.0 / 255.0,
	132.0 / 255.0,
	131.0 / 255.0,
	130.0 / 255.0,
	129.0 / 255.0,
	128.0 / 255.0,
	127.0 / 255.0,
	126.0 / 255.0,
	125.0 / 255.0,
	124.0 / 255.0,
	123.0 / 255.0,
	122.0 / 255.0,
	121.0 / 255.0,
	120.0 / 255.0,
	119.0 / 255.0,
	118.0 / 255.0,
	117.0 / 255.0,
	116.0 / 255.0,
	115.0 / 255.0,
	114.0 / 255.0,
	113.0 / 255.0,
	112.0 / 255.0,
	111.0 / 255.0,
	110.0 / 255.0,
	109.0 / 255.0,
	108.0 / 255.0,
	107.0 / 255.0,
	106.0 / 255.0,
	105.0 / 255.0,
	104.0 / 255.0,
	103.0 / 255.0,
	102.0 / 255.0,
	101.0 / 255.0,
	100.0 / 255.0,
	99.0 / 255.0,
	98.0 / 255.0,
	97.0 / 255.0,
	96.0 / 255.0,
	95.0 / 255.0,
	94.0 / 255.0,
	93.0 / 255.0,
	92.0 / 255.0,
	91.0 / 255.0,
	90.0 / 255.0,
	89.0 / 255.0,
	88.0 / 255.0,
	87.0 / 255.0,
	86.0 / 255.0,
	85.0 / 255.0,
	84.0 / 255.0,
	83.0 / 255.0,
	82.0 / 255.0,
	81.0 / 255.0,
	80.0 / 255.0,
	79.0 / 255.0,
	78.0 / 255.0,
	77.0 / 255.0,
	76.0 / 255.0,
	75.0 / 255.0,
	74.0 / 255.0,
	73.0 / 255.0,
	72.0 / 255.0,
	71.0 / 255.0,
	70.0 / 255.0,
	69.0 / 255.0,
	68.0 / 255.0,
	67.0 / 255.0,
	66.0 / 255.0,
	65.0 / 255.0,
	64.0 / 255.0,
	63.0 / 255.0,
	62.0 / 255.0,
	61.0 / 255.0,
	60.0 / 255.0,
	59.0 / 255.0,
	58.0 / 255.0,
	57.0 / 255.0,
	56.0 / 255.0,
	55.0 / 255.0,
	54.0 / 255.0,
	53.0 / 255.0,
	52.0 / 255.0,
	51.0 / 255.0,
	50.0 / 255.0,
	49.0 / 255.0,
	48.0 / 255.0,
	47.0 / 255.0,
	46.0 / 255.0,
	45.0 / 255.0,
	44.0 / 255.0,
	43.0 / 255.0,
	42.0 / 255.0,
	41.0 / 255.0,
	40.0 / 255.0,
	39.0 / 255.0,
	38.0 / 255.0,
	37.0 / 255.0,
	36.0 / 255.0,
	35.0 / 255.0,
	34.0 / 255.0,
	33.0 / 255.0,
	32.0 / 255.0,
	31.0 / 255.0,
	30.0 / 255.0,
	29.0 / 255.0,
	28.0 / 255.0,
	27.0 / 255.0,
	26.0 / 255.0,
	25.0 / 255.0,
	24.0 / 255.0,
	23.0 / 255.0,
	22.0 / 255.0,
	21.0 / 255.0,
	20.0 / 255.0,
	19.0 / 255.0,
	18.0 / 255.0,
	17.0 / 255.0,
	16.0 / 255.0,
	15.0 / 255.0,
	14.0 / 255.0,
	13.0 / 255.0,
	12.0 / 255.0,
	11.0 / 255.0,
	10.0 / 255.0,
	9.0 / 255.0,
	8.0 / 255.0,
	7.0 / 255.0,
	6.0 / 255.0,
	5.0 / 255.0,
	4.0 / 255.0,
	3.0 / 255.0,
	2.0 / 255.0,
	1.0 / 255.0,
	0.0 / 255.0};

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
	double *image_in_inverse = (double *)malloc(sizeof(double) * pixel_count);
	unsigned char *image_tmp = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);
	dst->cols = src.cols;
	dst->rows = src.rows;
	dst->data = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);

	//////////////////////////////////////////////////////////////////////////
	///// Initialization
	for (int p = 0; p < pixel_count; p++) {
		image_in_inverse[p] = Pixel_LUT[src.data[p]];
		image_tmp[p] = 255;
		dst->data[p] = 255;
	}

	//////////////////////////////////////////////////////////////////////////
	///// Find the number of Particle
	double CountParticle = 0;
	for (int p = 0; p < pixel_count; p++) {
		CountParticle = CountParticle + image_in_inverse[p];
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
						double t = image_in_inverse[p2] / (aa + b * b);
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

	free(image_in_inverse);
	free(image_tmp);

	return 1;
}
