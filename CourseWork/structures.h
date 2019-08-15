#pragma once

struct SLAE_CSR {
	// https://habr.com/ru/post/142662/ about struct in C/C++
	// A*x = b
	// ������� ����� ����� ����
	vector<double> A;
	// ������ ������ �����
	vector<double> b;
	// ������ ������� ( ���������� -0)
	vector<double> x;
	// ������� ��� �������� � � ������� CSr
	vector<int> iptr;
	vector<int> jptr;
	// ������� ��� ������ ������������ ������� � ���������� �
	void setX(const vector<double> &newX) { x = newX; }
	void coutSLAE(const string & name);

};

struct Mesh {
	// ���������� ��������� �����: ��� ������� N-�� �� ����� ����������� ����� �������� ����� (��� ����-�� ��-�� - 4 ����)
	vector<vector<int>> elemNodes;
	// ���������� ��������� ��: ��� ������� N-�� �� ����� ������������ ������ ��
	// (��� ����-�� ��-�� � ������� 1�� ������� - 4, 2�� - 9, 3�� - 16)
	vector<vector<int>> elemBF;
	// ������ ����� ��� DG - ��� ������� �������� ������ �������:) - �������������� � �������()
	// 1 - ����� �����, 2 - ����� ������, 3 - ����� �����, 4 - ����� ������
	vector<vector<int>> nghbrs;
	// ������ ��������� ��������� � ��, �� ����������� ��������� 
	vector<vector<vector<int>>> bndrElms;
	// ������ ������� (����������) �� �� �������� 
	// (��� ����-�� ������� ��������� ��������� ���������� ���������: 1 - �����, 2 - ������, 3 - ������, 4 - �������
	// � ������ ������� ������������ ��� �� �� ���� �������
	vector<vector<int>>  bndryBF;
	// ��� ������������� ����� ��� �������
	vector<double> X, Y;
	void coutMesh();
};

struct MeshDG :Mesh {
	vector<int> Elems;
	vector<vector<int>> neighbours;
	vector<double> X, Y;
};

struct FiniteElem {};

struct quadrElem{
	//: public FiniteElem {
	double x1, x2, y1, y2;
	double hx, hy;
	double middleX, middleY;
	double detJ;
	vector<vector<double>> J;
	void coutElem();
	vector<vector<double>> nVectors;
};

struct currentParamsValue {
	double lymbda;
	double gamma;
	vector<double> velocity;
};