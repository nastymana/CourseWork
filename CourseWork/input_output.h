#pragma once
//#ifndef INPUT_OUTPUT_H
//#define INPUT_OUTPUT_H
//#include"stdafx.h"

void outputf(ofstream &f,  double &a, const char s[100]);
void outputf_vectors(const vector<double> &A, const char namefile[100]);
void outputf_vectors(const vector<vector<int>> &A, const char namefile[100]);
void outputf_vectors(const vector<int> &A, const char namefile[100]);
void outputf_sets(const set<int> &A, const char namefile[100]);
void outputf_vectors(const vector<vector<double>> &A, const char namefile[100]);

vector<vector<vector<double>>> inputf_vector3d(ifstream &f, const char namefile[100]);
vector<vector<double>> inputf_vector2d(ifstream &f,const char namefile[100]);
vector<vector<int>> inputf_vector2i( ifstream &f, const char namefile[100]);
vector<double> inputf_vectord( ifstream &f, const char namefile[100]);
vector<int> inputf_vectori(ifstream &f, const char namefile[100]);
double inputf( ifstream &f, const char namefile[100]);

void outputf(ofstream &f, double &a, const char s, const wstring path);
void outputf_vectors(const vector<double> &A, const wstring &namefile, const wstring &path, const wstring newCatalog);
void outputf_vectors(const vector<string> &A, const wstring &namefile, const wstring &path, const wstring newCatalog);
void outputf_vectors(const vector<double> &A, const wstring namefile, const  wstring path);
void outputf_vectors(const vector<vector<int>> &A, const  wstring namefile, const wstring path);
void outputf_vectors(const vector<int> &A, const  wstring namefile, const wstring path);
void outputf_sets(const set<int> &A, const  wstring namefile, const wstring path);
void outputf_vectors(const vector<vector<double>> &A, const  wstring namefile, const  wstring path);

//#endif // !INPUT_OUTPUT_H