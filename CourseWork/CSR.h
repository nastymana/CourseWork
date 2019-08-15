#ifndef _CSR_H_
#define _CSR_H_
//#pragma once
//#include "stdafx.h"

vector<double> CSRmultiMV(const vector<double> &aelem, const vector<double> &x, const vector<int> &iptr, const vector<int> &jptr);

vector<double> CSRmultiMV(const SLAE_CSR* A, const vector<double> &x);

void CSRmultiMV_OMP(const vector<double> &aelem, const vector<double> &x,
	const vector<int> &iptr, const vector<int> &jptr, const int threadNum, vector<double>& result);

void coutMatrixCSR(const SLAE_CSR* A, const char name[30]);
#endif // !_CSR_H_

