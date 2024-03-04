#include "electrostatic.h"
#include "cv.hpp"
#include <Windows.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_THREAD 32

struct ES {
	int tmp;
	int current_iteration;
	int max_iterations;
	int enable_gridforce;
	int enable_shake;
	int color;
	int pixel_count;
	int particle_count;
	int cols;
	int rows;
	// struct CMat src;
	double *forcefield_y;
	double *forcefield_x;
	double *particles_y;
	double *particles_x;
	double *particles_new_y;
	double *particles_new_x;
	double *image_in_inverse;
	HANDLE semaphore;
	HANDLE mutex;
};

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

DWORD WINAPI forcefield_thread(struct ES *thread_data) {
	HANDLE mutex = thread_data->mutex;
	HANDLE semaphore = thread_data->semaphore;
	const int cols = thread_data->cols;
	const int rows = thread_data->rows;
	const double *image_in_inverse = thread_data->image_in_inverse;
	double *forcefield_y = thread_data->forcefield_y;
	double *forcefield_x = thread_data->forcefield_x;
	while (1) {
		WaitForSingleObject(mutex, INFINITE);
		int i = thread_data->tmp++;
		ReleaseMutex(mutex);
		if (i >= rows) {
			ReleaseSemaphore(semaphore, 1, NULL);
			return 1;
		}
		int p1 = i * cols;
		for (int j = 0; j < cols; j++) {
			forcefield_y[p1] = 0;
			forcefield_x[p1] = 0;
			int a = -i;
			for (int y = 0, p2 = 0; y < rows; y++) {
				int aa = a * a;
				int b = -j;
				for (int x = 0; x < cols; x++) {
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
		printf("%d\r", i);
	}
}

DWORD WINAPI iteration_thread(struct ES *thread_data) {
	HANDLE mutex = thread_data->mutex;
	HANDLE semaphore = thread_data->semaphore;

	const int max_iterations = thread_data->max_iterations;
	const int enable_gridforce = thread_data->enable_gridforce;
	const int enable_shake = thread_data->enable_shake;
	const int color = thread_data->color;
	const int cols = thread_data->cols;
	const int rows = thread_data->rows;
	const int particle_count = thread_data->particle_count;
	const int current_iteration = thread_data->current_iteration;
	double *forcefield_y = thread_data->forcefield_y;
	double *forcefield_x = thread_data->forcefield_x;
	double *particles_y = thread_data->particles_y;
	double *particles_x = thread_data->particles_x;
	double *particles_new_y = thread_data->particles_new_y;
	double *particles_new_x = thread_data->particles_new_x;

	while (1) {
		WaitForSingleObject(mutex, INFINITE);
		int current_particle = thread_data->tmp++;
		ReleaseMutex(mutex);
		if (current_particle >= particle_count) {
			ReleaseSemaphore(semaphore, 1, NULL);
			return 1;
		}
		double newpos_y = 0;
		double newpos_x = 0;
		double real_y = particles_y[current_particle] - (int)particles_y[current_particle];
		double real_x = particles_x[current_particle] - (int)particles_x[current_particle];

		// Attraction(by using bilinear interpolation)
		if (real_y == 0 && real_x == 0) {
			int p = (int)particles_y[current_particle] * cols + (int)particles_x[current_particle];
			newpos_y = forcefield_y[p];
			newpos_x = forcefield_x[p];
		} else {
			int bilinear_y1 = particles_y[current_particle];
			int bilinear_x1 = particles_x[current_particle];
			int bilinear_y2 = bilinear_y1 + 1;
			int bilinear_x2 = bilinear_x1 + 1;
			if (bilinear_y2 < rows && bilinear_x2 < cols) {
				double a = (double)bilinear_x2 - particles_x[current_particle];
				double b = (double)bilinear_y2 - particles_y[current_particle];
				double c = particles_x[current_particle] - (double)bilinear_x1;
				double d = particles_y[current_particle] - (double)bilinear_y1;
				int p11 = bilinear_y1 * cols + bilinear_x1;
				int p12 = bilinear_y1 * cols + bilinear_x2;
				int p21 = bilinear_y2 * cols + bilinear_x1;
				int p22 = bilinear_y2 * cols + bilinear_x2;
				newpos_y = forcefield_y[p11] * a * b + forcefield_y[p12] * c * b + forcefield_y[p21] * a * d + forcefield_y[p22] * c * d;
				newpos_x = forcefield_x[p11] * a * b + forcefield_x[p12] * c * b + forcefield_x[p21] * a * d + forcefield_x[p22] * c * d;
			}
		}

		// Repulsion
		for (int other_particle = 0; other_particle < particle_count; other_particle++) {
			if (current_particle != other_particle) {
				double instead_y = particles_y[other_particle] - particles_y[current_particle];
				double instead_x = particles_x[other_particle] - particles_x[current_particle];
				if (instead_y != 0 || instead_x != 0) {
					double t = color * (instead_y * instead_y + instead_x * instead_x);
					newpos_y -= instead_y / t;
					newpos_x -= instead_x / t;
				}
			}
		}

		// result (new position of particles)
		particles_new_y[current_particle] += 0.1 * newpos_y;
		particles_new_x[current_particle] += 0.1 * newpos_x;
		if (enable_gridforce) {
			// Add GridForce to find discrete particle locations
			// double real_y = particles_y[current_particle] - (int)particles_y[current_particle];
			// double real_x = particles_x[current_particle] - (int)particles_x[current_particle];
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
					particles_new_y[current_particle] += real_y * t;
				}
				if (real_x != 0) {
					// GridForce_X = 3.5 * real_x / (vector3 * (1 + pow(vector3, 8) * 10000));
					particles_new_x[current_particle] += real_x * t;
				}
			}
		}

		// Shake
		if (enable_shake && current_iteration % 10 == 0 && max_iterations > 64) {
			double t = (log10((double)max_iterations) / log10(2.0) - 6) * exp(current_iteration / 1000.0) / 10;
			particles_new_y[current_particle] += t;
			particles_new_x[current_particle] += t;
		}

		if (particles_new_y[current_particle] < 0) {
			particles_new_y[current_particle] = 0;
		} else if (particles_new_y[current_particle] >= rows) {
			particles_new_y[current_particle] = rows - 1;
		}
		if (particles_new_x[current_particle] < 0) {
			particles_new_x[current_particle] = 0;
		} else if (particles_new_x[current_particle] >= cols) {
			particles_new_x[current_particle] = cols - 1;
		}

		printf("%d\r", current_particle);
	}
}

int ElectrostaticHalftoning2010(struct CMat src, struct CMat *dst, int enable_initial_charge, int max_iterations, int enable_gridforce, int enable_shake, int enable_debug) {
	//////////////////////////////////////////////////////////////////////////
	///// exceptions
	// if(src.type()!=CV_8U){
	// 	CV_Error(CV_BadNumChannels,"[pixkit::halftoning::ElectrostaticHalftoning] image should be grayscale");
	// }
	if (enable_initial_charge != 0 && enable_initial_charge != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] enable_initial_charge should be 0 or 1");
		system("pause");
		exit(0);
	}
	if (max_iterations < 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] max_iterations should be bigger than 1");
		system("pause");
		exit(0);
	}
	if (enable_gridforce != 0 && enable_gridforce != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] enable_gridforce should be 0 or 1");
		system("pause");
		exit(0);
	}
	if (enable_shake != 0 && enable_shake != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] enable_shake should be 0 or 1");
		system("pause");
		exit(0);
	}
	if (enable_shake == 1 && max_iterations <= 64) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] max_iterations should be bigger than 64");
		system("pause");
		exit(0);
	}
	if (enable_debug != 0 && enable_debug != 1) {
		printf("[pixkit::halftoning::ElectrostaticHalftoning] Debug should be 0 or 1");
		system("pause");
		exit(0);
	}

	const int color = 3 - 1;
	unsigned char *particle_lut = (unsigned char *)malloc(sizeof(unsigned char) * (color + 1));
	char out_file[50];
	const int pixel_count = src.cols * src.rows;
	double *image_in_inverse = (double *)malloc(sizeof(double) * pixel_count);
	// unsigned char *image_tmp = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);
	int *image_particle = (int *)malloc(sizeof(int) * pixel_count);
	dst->cols = src.cols;
	dst->rows = src.rows;
	dst->data = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);

	//////////////////////////////////////////////////////////////////////////
	///// Initialization
	///// Find the number of Particle
	for (int i = 0; i <= color; i++) {
		particle_lut[i] = (double)i * 255 / color + 0.5;
	}
	int particle_count = pixel_count * 255;
	// memset(image_tmp, 255, sizeof(unsigned char) * pixel_count);
	memset(dst->data, 255, sizeof(unsigned char) * pixel_count);
	for (int p = 0; p < pixel_count; p++) {
		image_in_inverse[p] = Pixel_LUT[src.data[p]];
		// particle_count += image_in_inverse[p];
		particle_count -= src.data[p];
		image_particle[p] = color;
	}
	particle_count = particle_count * color / 255;
	printf("The number of black pixel(charge) = %d\n", (int)particle_count);

	//////////////////////////////////////////////////////////////////////////
	///// Initialize the Particle's position
	double *particles_y = (double *)malloc(sizeof(double) * particle_count);
	double *particles_x = (double *)malloc(sizeof(double) * particle_count);
	for (int i = particle_count; i > 0;) {
		// int RandY = rand() % src.rows;
		// int RandX = rand() % src.cols;
		int p = (rand() % src.rows) * src.cols + (rand() % src.cols);
		if (image_particle[p] <= 0) {
			continue;
		}
		// int RandNumber = rand() % 256;
		if (enable_initial_charge && rand() % 256 <= src.data[p]) {
			continue;
		}
		image_particle[p]--;
		i--;
	}
	if (enable_debug) {
		for (int p = 0; p < pixel_count; p++) {
			dst->data[p] = particle_lut[image_particle[p]];
		}
		sprintf(out_file, ".\\output\\0.bmp");
		cv_imwrite(out_file, *dst);
	}

	//////////////////////////////////////////////////////////////////////////
	///// Record the Particle's position
	for (int i = 0, p = 0, particle = 0; i < src.rows; i++) {
		for (int j = 0; j < src.cols; j++) {
			for (int k = image_particle[p]; k < color; k++) {
				particles_y[particle] = (double)i;
				particles_x[particle] = (double)j;
				particle++;
			}
			p++;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	///// Create Forcefield Table
	printf("Create Forcefield Table, \n");
	double *forcefield_y = (double *)malloc(sizeof(double) * pixel_count);
	double *forcefield_x = (double *)malloc(sizeof(double) * pixel_count);
	struct ES es;
	es.pixel_count = pixel_count;
	es.forcefield_y = forcefield_y;
	es.forcefield_x = forcefield_x;
	es.image_in_inverse = image_in_inverse;
	es.cols = src.cols;
	es.rows = src.rows;
	es.tmp = 0;
	es.mutex = CreateMutex(NULL, FALSE, TEXT("hMutex"));
	es.semaphore = CreateSemaphore(NULL, MAX_THREAD, MAX_THREAD, NULL);
	for (int i = 0; i < MAX_THREAD; i++) {
		WaitForSingleObject(es.semaphore, INFINITE);
		CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)forcefield_thread, &es, 0, NULL);
	}
	for (int i = 0; i < MAX_THREAD; i++) {
		WaitForSingleObject(es.semaphore, INFINITE);
	}
	ReleaseSemaphore(es.semaphore, MAX_THREAD, NULL);

	//////////////////////////////////////////////////////////////////////////
	///// process
	es.max_iterations = max_iterations;
	es.enable_gridforce = enable_gridforce;
	es.enable_shake = enable_shake;
	es.color = color;
	es.particle_count = particle_count;
	es.particles_y = particles_y;
	es.particles_x = particles_x;
	double *particles_new_y = (double *)malloc(sizeof(double) * particle_count);
	double *particles_new_x = (double *)malloc(sizeof(double) * particle_count);
	memcpy(particles_new_y, particles_y, sizeof(double) * particle_count);
	memcpy(particles_new_x, particles_x, sizeof(double) * particle_count);
	es.particles_new_y = particles_new_y;
	es.particles_new_x = particles_new_x;
	for (int iteration = 1; iteration <= max_iterations; iteration++) {
		printf("Iterations %d\n", iteration);

		WaitForSingleObject(es.mutex, INFINITE);
		es.tmp = 0;
		es.current_iteration = iteration;
		ReleaseMutex(es.mutex);
		for (int i = 0; i < MAX_THREAD; i++) {
			WaitForSingleObject(es.semaphore, INFINITE);
			CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)iteration_thread, &es, 0, NULL);
		}
		for (int i = 0; i < MAX_THREAD; i++) {
			WaitForSingleObject(es.semaphore, INFINITE);
		}
		ReleaseSemaphore(es.semaphore, MAX_THREAD, NULL);
		printf("\n");

		memcpy(particles_y, particles_new_y, sizeof(double) * particle_count);
		memcpy(particles_x, particles_new_x, sizeof(double) * particle_count);

		// Output
		for (int p = 0; p < pixel_count; p++) {
			dst->data[p] = 255;
			image_particle[p] = color;
		}
		for (int current_particle = 0; current_particle < particle_count; current_particle++) {
			int out_Y = particles_y[current_particle] + 0.5;
			int out_X = particles_x[current_particle] + 0.5;
			if (out_Y >= src.rows) {
				out_Y = src.rows - 1;
			}
			if (out_X >= src.cols) {
				out_X = src.cols - 1;
			}
			image_particle[out_Y * src.cols + out_X]--;
		}

		for (int p = 0; p < pixel_count; p++) {
			if (image_particle[p] < 0) {
				image_particle[p] = 0;
			}
			dst->data[p] = particle_lut[image_particle[p]];
		}

		if (enable_debug) {
			sprintf(out_file, ".\\output\\%d.bmp", iteration);
			cv_imwrite(out_file, *dst);
		}
	}

	// dst = dst->clone();

	free(image_in_inverse);
	// free(image_tmp);

	return 1;
}
