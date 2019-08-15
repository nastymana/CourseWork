#ifndef _COUTOBJECTS_H_
#define _COUTOBJECTS_H_


//#pragma once
//#ifndef COUTOBJECTS_H
//#define COUTOBJECTS_H


void coutVector(const vector<vector<double>> &x, const char name[30]);
void coutVector(const vector<double> &x, const char name[30], char direct, const int &precsn=2);
void coutVector(const vector<vector<int>> &x, const char name[30]);
void coutVector(const vector<int> &x, const char name[30], char direct);
void coutSet(const set<int> &x, const char name[30], char direct);
void coutVectorsets(const vector<set<int>> &A, const char name[50]);
void coutSet(const set<double> &A, const char name[50], char direct);

#endif // !COUTOBJECTS_H
