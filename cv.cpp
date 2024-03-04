#include "electrostatic.h"
#include "cv.hpp"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;
// using namespace std;

int cv_imread(const char *path, struct CMat *csrc) {
	cv::Mat src;
	src = cv::imread(path, cv::ImreadModes::IMREAD_GRAYSCALE);
	if (!src.empty()) {
		int pixel_count = src.rows * src.cols;
		csrc->data = (unsigned char *)malloc(sizeof(unsigned char) * pixel_count);
		csrc->cols = src.cols;
		csrc->rows = src.rows;
		for (int i = 0; i < pixel_count; i++) {
			csrc->data[i] = src.data[i];
		}
		return 1;
	}
	return 0;
}

int cv_imwrite(const char *path, struct CMat dst) {
	// cv::Mat dst;
	Mat real_dst(dst.rows, dst.cols, CV_8UC1);
	int pixel_count = dst.rows * dst.cols;
	for (int i = 0; i < pixel_count; i++) {
		real_dst.data[i] = dst.data[i];
	}
	if (cv::imwrite(path, real_dst)) {
		return 1;
	}
	return 0;
}

void cv_imshow(const char *label, struct CMat dst) {
	Mat real_dst(dst.rows, dst.cols, CV_8UC1);
	int pixel_count = dst.rows * dst.cols;
	for (int i = 0; i < pixel_count; i++) {
		real_dst.data[i] = dst.data[i];
	}
	cv::imshow(label, real_dst);
}

int cv_waitKey(int delay) {
	return cv::waitKey(delay);
}