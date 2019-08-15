#pragma once

struct SLAE_CSR {
	// https://habr.com/ru/post/142662/ about struct in C/C++
	// A*x = b
	// матрица левой части СЛАУ
	vector<double> A;
	// вектор правой части
	vector<double> b;
	// вектор решения ( изначально -0)
	vector<double> x;
	// массивы для хранения А в формате CSr
	vector<int> iptr;
	vector<int> jptr;
	// функция для явного присваивания решения в переменную х
	void setX(const vector<double> &newX) { x = newX; }
	void coutSLAE(const string & name);

};

struct Mesh {
	// глобальная нумерация узлов: для каждого N-го КЭ сетки представлен набор основных узлов (для прям-ых эл-ов - 4 узла)
	vector<vector<int>> elemNodes;
	// глобальная нумерация БФ: для каждого N-го КЭ сетки представлены наборы БФ
	// (для прям-ых эл-ов с базисом 1го порядка - 4, 2го - 9, 3го - 16)
	vector<vector<int>> elemBF;
	// больше нужно для DG - для каждого элемента список соседей:) - принадлежность к границе()
	// 1 - сосед слева, 2 - сосед справа, 3 - сосед снизу, 4 - сосед сверху
	vector<vector<int>> nghbrs;
	// список граничных элементов с БФ, кт принадлежат конретной 
	vector<vector<vector<int>>> bndrElms;
	// список номеров (глобальных) БФ по границам 
	// (для прям-ой области нумерация граничных плоскостей следующая: 1 - левая, 2 - правая, 3 - нижняя, 4 - верхняя
	// в каждом векторе перечисленны все БФ на этой границе
	vector<vector<int>>  bndryBF;
	// для прямоугольной сетки два вектора
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