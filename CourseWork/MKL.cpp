#include "pch.h"
#include"MKL.h"
inline void _dcopy(size_t n, const double *x, double *y) {
	MKL_INT m = (MKL_INT)n;
	MKL_INT inc = 1;
	//MKL_Complex16 * x1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(x));
	//MKL_Complex16 * y1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(y));
	dcopy(&m, x, &inc, y, &inc);
}
void _daxpby(int &N, const double a, const double *x, const int inc, const double b, const double *y, double *result) {
	_dcopy(N, y, result);
	daxpby(&N, &a, x, &inc, &b, result, &inc);
}
void _daxpbypcz(int &N, const double a, const double *x, const int inc, const double b, double *y, const double c, const double *z, double *result) {

	daxpby(&N, &a, x, &inc, &b, y, &inc);
	_dcopy(N, y, result);
	const double one = 1;
	daxpby(&N, &c, z, &inc, &one, result, &inc);
}
void dcopyv(const vector<double> &a, double *d, int threadNum) {
	int i = 0;

#pragma omp parallel shared(a, d) private(i)num_threads(threadNum) 
	{
#pragma omp for

		for (i = 0; i < int(a.size()); i++)
		{
			d[i] = a[i];
		}
	}
}

void icopyv(const vector<int> &a, int *d, int threadNum) {
	int i = 0;

#pragma omp parallel shared(a, d) private(i)num_threads(threadNum) 
	{
#pragma omp for

		for (i = 0; i < int(a.size()); i++)
		{
			d[i] = a[i];
		}
	}

}
// y := a*x + y
inline void _daxpy(size_t n, const double *a, const double *x, double * result) {
	MKL_INT m = (MKL_INT)n;
	MKL_INT inc = 1;
	//	MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(&a));
	//	MKL_Complex16 * x1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(x));
	//	MKL_Complex16 * y1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(y));
	daxpy(&m, a, x, &inc, result, &inc);
}

//c := (a*) * b
double _ddotc(size_t n, const double * a, const double * b)
{
	double result;
	const int n1 = n;
	const MKL_INT inc = 1;
	//	MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(a));
	//	MKL_Complex16 * b1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(b));
	//	MKL_Complex16 * pres = reinterpret_cast<MKL_Complex16 *>(&result);
	return ddot(&n1, a, &inc, b, &inc);

}
//c := |a * a|
double _dnrm2(size_t n, const double * a)
{
	MKL_INT n1 = (MKL_INT)n;
	MKL_INT inc = 1;
	return dnrm2(&n1, a, &inc);
}
//c := |a * b|
double scal_mul(size_t n, const double* a, const double * b)
{
	double d_p = 0.0;
#pragma omp parallel for reduction(+ : d_p)
	for (int i = 0; i < (int)n; i++)
		d_p += a[i] * b[i];
	return d_p;
}
//x := A * f
void mul_matrix(size_t n, const double* aa, const int *ia, const  int *ja, const double* f, double * x)
{
	MKL_INT m = (MKL_INT)n;
#pragma warning(disable : 4996)
	mkl_dcsrgemv("N", &m, aa, ia, ja, f, x);
}
