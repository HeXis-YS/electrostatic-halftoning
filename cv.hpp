#ifdef __cplusplus
extern "C" {
#endif
struct CMat {
	unsigned char *data;
	int rows;
	int cols;
};
int cv_imread(const char *, struct CMat *);
int cv_imwrite(const char *, struct CMat);
void cv_imshow(const char *label, struct CMat);
int cv_waitKey(int);
#ifdef __cplusplus
}
#endif
