#include "pch.h"
#include"coutObjects.h"

////2

void coutVector(const vector<vector<double>> &x, const char name[30])
{
	for (int j = 0; j < int(x.size()); j++)
	{
		cout << name << " [" << j << "] = { ";
		for (int i = 0; i < int(x[0].size()) - 1; i++)
		{
			cout <<setprecision(2)<< x[j][i] << ",  ";
		}
		cout << setprecision(2)<<x[j].back() << " }" << endl;
	}
	cout << endl;
}
void coutVector(const vector<double> &x, const char name[30], char direct, const int &precsn)
{
	if (direct == 'H')
	{
		cout << name << " = { ";
		for (int i = 0; i <= int(x.size()) ; i++)
		{
			cout<< setprecision(precsn) << x[i] << ",  ";
		}

		cout << x.back() << " }" << endl;
	}
	else
	{
		for (int i = 0; i <= int(x.size()); i++)
		{
			cout << x[i] << endl;
		}
	}
}
void coutVector(const vector<vector<int>> &x, const char name[30])
{
	//int q = x.size();
	//int n = x[0].size();
	for (int j = 0; j < int(x.size()); j++)
	{
		cout << name << " [" << j << "] = { ";
		for (int i = 0; i < int(x[j].size()) - 1; i++)
		{
			cout << x[j][i] << ",  ";
		}
		cout << x[j].back() << " }" << endl;
	}
}
void coutVector(const vector<int> &x, const char name[30], char direct)
{
	if (direct == 'H')
	{
		cout << name << " = { ";
		for (int i = 0; i < int(x.size()) - 1; i++)
		{
			cout << x[i] << ",  ";
		}

		cout << x.back() << " }" << endl;
	}
	else
	{
		for (int i = 0; i < int(x.size()); i++)
		{
			cout << x[i] << endl;
		}
	}
}
void coutSet(const set<int> &x, const  char name[30], char direct)
{
	cout << name;
	if (direct == 'H')
	{
		cout << " = {";
		for (set<int>::iterator it = x.begin(); it != x.end(); ++it)
		{
			cout << setw(5) << *it;
		}

		cout << " }" << endl;
	}
	else
	{
		for (set<int>::iterator it = x.begin(); it != x.end(); ++it)
		{
			cout << *it << endl;
		}
	}
}
void coutVectorsets(const vector<set<int>> &A, const char name[50])
{
	cout << "\n" << name << "\n";
	for (int i = 0; i <int(A.size()); i++)
	{
		cout << i << "::  ";
		for (set<int>::iterator it = A[i].begin(); it != A[i].end(); ++it)
		{
			cout << *it << ' ';
		}
		cout << endl;
	}
}
void coutSet(const set<double> &A, const char name[50], char direct)
{
	cout << "\n";
	if (direct == 'H')
	{
		cout << name << " = { ";
		for (set<double>::iterator it = A.begin(); it != A.end(); ++it)
		{
			cout << *it << ' ';
		}
		cout << " }" << endl;
	}
	else
	{
		cout << name << endl;
		for (set<double>::iterator it = A.begin(); it != A.end(); ++it)
		{
			cout << *it << endl;
		}
	}
}
