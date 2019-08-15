//#ifndef VECTORMATHOPERATIONS_H
//#define VECTORMATHOPERATIONS_H
//#include "stdafx.h"
#pragma once
// ФУНКЦИИ ЭЛЕМЕНТАРНЫХ МАТ ОПЕРАЦИЙ НАД ВЕКТОРАМИ И МАТРИЦАМИ
void innerProduct(const vector<vector<double>> &A, const vector<vector<double>> &B,
	vector<vector<double>> &result);
double innerProduct(const vector<double> &a, const vector<double> &b);
double inner_product_OMP(const vector<double> &a, const vector<double> &b, const int threadNum);
vector<double> innerProduct(const vector<vector<double>> &A, const vector<double> &B);
void innerProduct(const vector<vector<double>> &A, const vector<vector<double>> &B, vector<vector<double>> &result);
double summVectorElms(const vector<vector<double>> &A);
double summVectorElms(const vector<double> &A);
double vectorNorma(const vector<double> &a);
double vector_norma_OMP(const vector<double> &a, const int threadNum);
double angle_between_vectors(const vector<double> &a, const vector<double> &b);
vector<double> vector_normalyze(const vector<double> &a);
vector<vector<double>> operator+ (const vector<vector<double>> &A, const vector<vector<double>> &B);
vector<vector<int>> operator+ (const vector<vector<int>> &A, const vector<vector<int>> &B);
vector<double> operator+ (const vector<double> &A, const vector<double> &B);
void ax_plus_by_OMP(int a, const vector<double> &A, int b, const vector<double> &B, vector<double> &result, const int threadNum);
void ax_plus_by_plus_cz_OMP(int a, const vector<double> &A, int b, const vector<double> &B, 
	int c, const vector<double> &C, vector<double> &result, const int threadNum);

vector<vector<double>> operator- (const vector<vector<double>> &A,const vector<vector<double>> &B);
vector<vector<int>> operator- (const vector<vector<int>> &A, const vector<vector<int>> &B);
vector<double> operator- (const vector<double> &A, const vector<double> &B);

vector<double> operator* (const vector<double> &A, double c);
vector<vector<double>> operator* (const vector<vector<double>> &A, double c);
vector<vector<double>> operator* (const vector<vector<double>> &A, const vector<vector<double>> &b);
vector<double> operator* (double c, const vector<double> &A);
vector<vector<double>> operator* (double c, const vector<vector<double>> &A);

vector<double> operator* (const vector<double> &A, const vector<double> &B);
vector<vector<double>> operator* (const vector<vector<double>> &A, const vector<double> &c);

void sewvectors(vector<double> &a, const vector<double> &b);
vector<double> pushvector(const vector<double> &a, const double push);
void sewvectors(vector<int> &a, const vector<int> &b);

vector<double> multiply_vectorsn_to_nn(const  vector<double> &a, const  vector<double> &b);
vector<vector<double>> Uptomirror(const vector<vector<double>> &A);
vector<vector<double>> transposeMatrix(const vector<vector<double>>A);
double lenght(const double &x1, const  double &y1, const double &x2, const double &y2);
void sewvectors(const vector<vector<double>> &a, const vector<vector<double>>&b);
double average(const vector<double> &a);

void reshape_vector(const vector<vector<int>> &a, vector<int> &result);
// возвращает для исходного вектора вектор с членами в степени deg
void degreeVector(const vector<double> &a, vector<double> &result, const int &deg=2);
//#endif