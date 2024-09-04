#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
	unsigned char *data;
	unsigned int rows;
	unsigned int cols;
	unsigned char inverse;
} CMat;
int cv_imread(const char *, CMat *);
int cv_imwrite(const char *, const CMat *);
void cv_imshow(const char *label, const CMat *);
int cv_waitKey(int);
#ifdef __cplusplus
}
#endif
