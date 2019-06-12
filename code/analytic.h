
typedef struct {
	double x;
	double y;
	double _Complex En;
	int reim;
	int diagt; // 0 for diagonal, 1 for 12-type off diagonal, 2 for 21-type offdiagonal
}intg_params;

typedef struct {
	double *x;
	double *y;
	double _Complex En;
	int *diagt; // 0 for diagonal, 1 for 12-type off diagonal, 2 for 21-type offdiagonal
}intgac_params;



double graphenegfa(double _Complex En, int a1, int a2, int reim, int diagt);

double intga(double kx, void *p);
	
void getS1S2a (double _Complex *S1, double _Complex *S2, double x, double y, double kx, double _Complex cosky1, double _Complex cosky2, double _Complex En);



double _Complex *graphenegfac(double _Complex En, int *a1, int *a2, int *diagt, int dim);
void intgac(unsigned ndim, const double *kx, void *p, unsigned fdim, double *fval);

