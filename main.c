#include "electrostatic.h"
#include "cv.hpp"

void main() {

	struct CMat src, dst;

	// load image
	if (cv_imread("lena.bmp", &src)) {
		// process
		if (ElectrostaticHalftoning2010(src, &dst, 1, 8, 1, 0, 1)) {
			// write output
			cv_imwrite("../../output.bmp", dst);

			// show image
			cv_imshow("src", src);
			cv_imshow("dst", dst);
			cv_waitKey(0);
		}
	}
}