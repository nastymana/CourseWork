#pragma once


Mesh generateQuadElementMesh(const string &DGorCG, const int basisRank, const vector<double> &indata);

void scaledGrid(int N_node, double coef, vector<double> &x);

// отработка алгоритма для полной матрицы (не CSR)
// на входе нужен файл с переписью элементов: номер элемента соотв-ет номеру ячейки. в ячейке хранится номер первой базисной функции на эл-те. 
void quadGridForQuadDomain(const vector<double> &indata, vector<double> &X, vector<double> &Y);

vector<vector<int>> elementListDG(const int &elemRank, const vector<double> &X, const vector<double> &Y);

vector<vector<int>> elemsNghbrsList(const int &Nelem, const int &NnodeX, const int &NnodeY);

vector<vector<double>> nodesForQuadGrid(const vector<double> &X, const vector<double>& Y);

vector<vector<int>> edgesForQuadGrid(const int &nodeX, const int &nodeY);

vector<vector<int>> elementBaseNodes(const int &nodeX, const int &nodeY);

vector<vector<int>> elementBF(const int &elemRank, const int& NElmRows,
	const vector<vector<int>> &elemNodes);

vector<vector<int>> elementEdges(const int &nodeX, const int &nodeY);

void elemsBndryList(const int &basisRank, const int &NnodeX, const int &NnodeY,
	const vector<vector<int>> &neighbours, const vector<vector<int>> elemsBF,
	vector<vector<vector<int>>> &bndrElems);

void CGbndryBF(const int &basisRank, const int& nX, const int& nY, const vector<vector<int>> &elemBF, vector<vector<int>> &bndrBF);

void DGelemsBndryList(const int &Nelem, const int &NBFon_elm, const int &NnodeX, const int &NnodeY,
	const vector<vector<int>> &neighbours, vector<vector<int>> &bndrElems, vector<vector<int>> &bndryBF);
