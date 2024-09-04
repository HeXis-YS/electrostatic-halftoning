
#include "cv.hpp"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

typedef unsigned char uint8_t;
typedef unsigned int uint_t;

using namespace cv;

static void copy_data(uint8_t *dst, const uint8_t *src, uint_t size, uint8_t inverse) {
	if (inverse) {
		for (int i = 0; i < size; i++) {
			dst[i] = ~src[i];
		}
	} else {
		memcpy(dst, src, sizeof(uint8_t) * size);
	}
}

int cv_imread(const char *path, CMat *csrc) {
	Mat src = imread(path, ImreadModes::IMREAD_GRAYSCALE);
	if (!src.empty()) {
		const uint_t pixel_count = src.rows * src.cols;
		uint_t strength = 0;
		csrc->data = (uint8_t *)malloc(sizeof(uint8_t) * pixel_count);
		csrc->cols = src.cols;
		csrc->rows = src.rows;
		for (int i = 0; i < pixel_count; i++) {
			strength += src.data[i];
		}
		if (strength <= (pixel_count * 255) / 2) {
			csrc->inverse = 1;
		} else {
			csrc->inverse = 0;
		}
		copy_data(csrc->data, src.data, pixel_count, csrc->inverse);
		return 1;
	}
	return 0;
}

int cv_imwrite(const char *path, const CMat *cdst) {
	Mat dst(cdst->rows, cdst->cols, CV_8UC1);
	copy_data(dst.data, cdst->data, cdst->rows * cdst->cols, cdst->inverse);
	if (imwrite(path, dst)) {
		return 1;
	}
	return 0;
}

void cv_imshow(const char *label, const CMat *cdst) {
	Mat dst(cdst->rows, cdst->cols, CV_8UC1);
	copy_data(dst.data, cdst->data, cdst->rows * cdst->cols, cdst->inverse);
	imshow(label, dst);
}

int cv_waitKey(int delay) {
	return waitKey(delay);
}