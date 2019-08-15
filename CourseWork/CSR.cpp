#include "pch.h"
#include"CSR.h"

vector<double> CSRmultiMV(const vector<double> &aelem, const vector<double> &x, const vector<int> &iptr, const vector<int> &jptr) {
	int rank = x.size();
	vector<double> result(rank, 0.);
	for (int i = 0; i < rank; i++) {
		for (int j = (int(iptr[i])); j<int(iptr[i + 1]); j++) {
			result[i] += x[jptr[j]] * aelem[j];
		}
	}
	return result;
}

vector<double> CSRmultiMV(const SLAE_CSR* A, const vector<double> &x ) {
	int rank = x.size();
	vector<double> result(rank, 0.);
	for (int i = 0; i < rank; i++) {
		for (int j = (int(A->iptr[i])); j<int(A->iptr[i + 1]); j++) {
			result[i] += x[A->jptr[j]] * A->A[j];
		}
	}
	return result;
}

void CSRmultiMV_OMP(const vector<double> &aelem, const vector<double> &x,
	const vector<int> &iptr, const vector<int> &jptr, const int threadNum,
	vector<double>& result) {
	const int rank = x.size();
	int i = 0, j = 0;

#pragma omp parallel shared(iptr,jptr,aelem,x,result) private(i,j)num_threads(threadNum) 
	{
#pragma omp for
		for (i = 0; i < rank; ++i) {
			//cout << "parallelized with " << omp_get_num_threads() << "threads,  i = " << i << endl;
			result[i] = 0;
			for (j = (int(iptr[i])); j<int(iptr[i + 1]); ++j) {
			//	cout << "parallelized with %d threads" << omp_get_num_threads() << ",  j = " << j << endl;
				result[i] += x[jptr[j]] * aelem[j];
			}
		}
	}

}

