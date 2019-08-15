
#include "pch.h"
#include"Mesh.h"

//coutVector(indata, "indata",'H');
//vector<double>X, Y, Z;
Mesh generateQuadElementMesh(const string &DGorCGitem, const int basisRank, const vector<double> &indata)
{
	//case 2: 2d равномерна€ сетка в пр€моугольной области (пр€моугольные  Ё) 
	cout << "\n generateQuadElementMesh\n";
	Mesh mesh;
	
//	int innerboundarynZ, innerboundarynX, innerboundarynY;

	vector<vector<double>> nodes;
	int dx = indata[0], dy = indata[dx * 4 + 1];
	//cout << "dx= " << dx << "  dy = " << dy << endl;

	quadGridForQuadDomain(indata, mesh.X, mesh.Y);
	int nodeX = mesh.X.size(), nodeY = mesh.Y.size();
	int Nelem = (nodeX - 1)*(nodeY - 1);
	int Nnode = nodeX*nodeY;
	
	//cout << "nodeX = " << nodeX << "  nodeY = " << nodeY << endl;

	mesh.elemNodes = elementBaseNodes(nodeX, nodeY);
	if (DGorCGitem == "DG") { mesh.elemBF = elementListDG(basisRank, mesh.X, mesh.Y); }
	else if (DGorCGitem == "CG")	mesh.elemBF = elementBF(basisRank, nodeY - 1, mesh.elemNodes);
	mesh.nghbrs = elemsNghbrsList(Nelem, nodeX, nodeY);
	CGbndryBF(basisRank, nodeX, nodeY, mesh.elemBF, mesh.bndryBF);	
	elemsBndryList(basisRank, nodeX, nodeY, mesh.nghbrs,mesh.elemBF, mesh.bndrElms); // переделать функцию elemsBndryList
	
	//Domains.push_back(innerboundarynY*(nodeX - 1));
	//Domains.push_back(Nelem - 1);

	/*
	set<int> &innerBoundaryFaces, set<int> &innerBoundaryEdges, set<int> &innerBoundaryNodes, set<int> &innerboundaryElems,		
	vector<int> &Domains
	*/
	
	return mesh;
}

void scaledGrid(int N_node, double coef, vector<double> &x) {
	cout << "\n scaledGrid\n";
	cout << N_node << ", " << coef << endl;
	double h = 0, summ = 0;
	for (int i = 0; i < N_node - 1; i++)
	{
		summ += pow(coef, i);
	}
	h = 1 / summ;
	
	x = vector<double>(N_node, 0.);
//	coutVector(x, "X", 'H');
	for (int i = 1; i < N_node - 1; i++)
	{
		x[i] = (x[i - 1] + h * pow(coef, i - 1));
	//	coutVector(x, "X", 'H');
	}
	cout << x[N_node - 1] << ", " <<x.size()<<endl;
	x[N_node-1] = 1.;

	//coutVector(x, "X", 'H');
	//if (coef < 0.) {
	//	x_help.push_back(0.);
	//	for (int i = 1; i < n_node - 1; i++)
	//	{
	//		x_help.push_back(1 - x[n_node - i - 1]);
	//	}
	//	x_help.push_back(1);
	//	//int n = x.size();
	//	return x_help;
	//}
	//int n = x.size();
}

// отработка алгоритма дл€ полной матрицы (не CSR)
// на входе нужен файл с переписью элементов: номер элемента соотв-ет номеру €чейки. в €чейке хранитс€ номер первой базисной функции на эл-те. 
void quadGridForQuadDomain(const vector<double> &indata, vector<double> &X, vector<double> &Y) {
	// 2D согласованна€ сетка в пр€моугольной области (пр€моугольные  Ё) 
	cout << "\n quadGridForQuadDomain\n";
	int dx = int(indata[0]), dy = int(indata[dx * 4 + 1]);
	int innerbndry_nX;//, innerbndry_nY;
	std::cout << "dx= " << dx << "  dy = " << dy << std::endl;
	for (int i = 0; i < dx; i++)
	{
		const double &p = indata[i * 4 + 1];
		cout << *(&p) << std::endl;
		vector<double> help;
//		cout << "int(*(&p + 2)) + 1 = " << int(*(&p + 2)) + 1 << endl;
		scaledGrid(int(*(&p + 2)) + 1, *(&p + 3), help);
	//	coutVector(help, "help", 'H');
//		std::cout << "-p + *(&p + 1) = " << -p + *(&p + 1) << std::endl;
		help = help * (-p + *(&p + 1));
		help = pushvector(help, (p));
	//	coutVector(help, "help", 'H');
		if (dx != 0)
		{
			innerbndry_nX = X.size();
		}
		if (i != (dx - 1))
		{
			help.pop_back();
		}

		sewvectors(X, help);
	//	coutVector(X, "coordX", 'H');
	}
//	coutVector(X, "coordX", 'H');
	std::vector<int> Domains = { 0 };
	//	std::cout << "domains\n";
	for (int i = 0; i < dy; i++)
	{
		int nYelem = int(indata[dx * 4 + i * 4 + 4]);
		std::cout << nYelem << std::endl;
		int nXelem = X.size() - 1;
		int nst = nXelem * nYelem;
		Domains.push_back(Domains[i] + nst);
	}
	//coutVector(Domains, "domains", 'H');
	for (int i = 0; i < dy; i++)
	{
		const double &p = indata[dx * 4 + i * 4 + 2];
		//double &p = indata[dx * 4 + dy * 4 + i * 4 + 3];
//		std::cout << *(&p) << std::endl;
		std::vector<double> help;
		scaledGrid(int(*(&p + 2)) + 1, *(&p + 3), help);
		help = help * (-p + *(&p + 1));
		//	coutVector(help, "coordZ", 'H');
		help = pushvector(help, (p));
		//coutVector(help, "coordZ", 'H');
		if (dy != 0)
		{
			//		innerboundarynY = Y.size();
		}
		if (i != (dy - 1))
		{
			help.pop_back();
		}
		sewvectors(Y, help);
		//coutVector(Y, "coordY", 'H');
	}

}

vector<vector<int>> elementListDG(const int &elemRank,  const vector<double> &X, const vector<double> &Y) {
	//2d равномерна€ сетка в пр€моугольной области (пр€моугольные  Ё) 
	//int dx = indata[0], dy = indata[dx * 4 + 1];
	cout << "\n elementListDG\n";
	vector<vector<double>> nodes;
	vector<vector<int>>	Elems; // только узлы

	int Nelem = (X.size() - 1)*(Y.size() - 1),
		Nnode = X.size()*Y.size(),
		NnodeX = X.size(), NnodeY = Y.size(),
		Nedge = (NnodeX - 1)*NnodeY + (NnodeY - 1)*NnodeX;
	std::cout << "nodeX = " << NnodeX << "  nodeY = " << NnodeY << std::endl;
	/*for (int j = 0; j < int(Y.size()); j++)
	{
		for (int k = 0; k < int(X.size()); k++)
		{
			nodes.push_back({ X[k],Y[j] });
		}
	}*/
	// исправить

	//Domains.push_back(innerboundarynY*(nodeX - 1));
	//Domains.push_back(Nelem - 1);
	//std::vector<std::vector<int>> edges;
	//for (int j = 0; j < NnodeY; j++)
	//{
	//	for (int k = 0; k < NnodeX - 1; k++)
	//	{
	//		int n0 = NnodeX * j + k;
	//		edges.push_back({ n0, n0 + 1 });// грани в плоскости паралл XY кроме 	
	//	}
	//	if (j != int(Y.size() - 1))
	//	{
	//		for (int k = 0; k < NnodeX; k++)
	//		{
	//			int n0 = NnodeX * j + k;
	//			//int n2 = nodeX*(j + 1) + k;
	//			edges.push_back({ n0,n0 + NnodeX });
	//		}
	//	}
	//}
	
	//coutSet(boundaryEdges, "boundaryEdges", 'H');
	// простой случай, когда на всех элементах одинаковый базис - 4 билин функции
	Elems;
	int bf1 = 0, nElem;
	for (int j = 0; j < NnodeY - 1; j++)
	{
		for (int k = 0; k < NnodeX - 1; k++)
		{
			nElem = (NnodeX - 1)*j + k;
				if (elemRank == 1) {
					bf1 = nElem * 4;
					Elems.push_back({ bf1, bf1 + 1, bf1 + 2, bf1 + 3 });
				}
				else if (elemRank == 2) {
					bf1 = nElem * 9;
					Elems.push_back({ bf1, bf1 + 1, bf1 + 2, bf1 + 3, bf1 + 4, bf1 + 5, bf1 + 6, bf1 + 7, bf1 + 8 });
				}

		}
	}

	return Elems;
	//char charsize[5]; what is it?

}

vector<vector<int>> elemsNghbrsList(const int &Nelem, const int &NnodeX, const int &NnodeY) {
	//neighbours
	cout << "\n elemsNghbrsList\n";
	vector<vector<int>> neighbours(Nelem);
	int NelemX = NnodeX - 1, NelemY = NnodeY - 1;
	for (int i = 0; i < Nelem; i++){
		//std::cout << "nElem " << i << "\n";
		if (i / NelemX == 0) {// нижний р€д элементов
			if ((i != 0) & ((i + 1) % NelemX != 0)) {
				neighbours[i] = { 3, i - 1, i + 1, -1, i + NelemX };// остальные нижние элементы
				// neighbours.push_back({ 3, i - 1, i + 1, -1, i + NelemX });// остальные нижние элементы
			}
			else {
				if (i == 0) {
					neighbours[i] = { 13, -1, i + 1, -1, i + NelemX };
					// neighbours.push_back({ 13, -1, i + 1, -1, i + NelemX });
				} // левый нижний угол
				else if(i==NelemX-1) {
					neighbours[i] = { 23, i - 1, -1, -1, i + NelemX };
					// neighbours.push_back({ 23, i - 1, -1, -1, i + NelemX });
				}// правый нижний угол
			}
		}
		else {
			if (i%NelemX == 0) {//лева€ граница 
				if (i / NelemX + 1 == NelemY) { 
					neighbours[i] = { 14, -1, i + 1, i - NelemX, -1 };
					// neighbours.push_back({ 14, -1, i + 1, i - NelemX, -1 }); 
				}// верхний левый элемент
					
				else { 
					neighbours[i] = { 1, -1, i + 1, i - NelemX, i + NelemX };
					//	neighbours.push_back({ 1, -1, i + 1, i - NelemX, i + NelemX }); 
				}
			}
			else {
				if ((i + 1) % NelemX == 0) {// права€ граница
					if ((i + 1) / NelemX == NelemY) {
						neighbours[i] = { 24, i - 1, -1, i - NelemX, -1 };
						// neighbours.push_back({ 24, i - 1, -1, i - NelemX, -1 });
					}// верхний правый элемент
					else {
						neighbours[i] = { 2, i - 1, -1, i - NelemX, i + NelemX };
						// neighbours.push_back({ 2, i - 1, -1, i - NelemX, i + NelemX });
					}
				}
				else {
					if ((i / NelemX > 1) & (i / NelemX + 1 == NelemY)) {
						neighbours[i] = { 4, i - 1, i + 1, i - NelemX, -1 };
						// neighbours.push_back({ 4, i - 1, i + 1, i - NelemX, -1 });
					}
					else { 
						neighbours[i] = { 0, i - 1, i + 1, i - NelemX, i + NelemX };
						// neighbours.push_back({ 0, i - 1, i + 1, i - NelemX, i + NelemX }); 
					}
				}
			}
		}
			
	}
//	coutVector(neighbours, "neighbours");
	return neighbours;
}

vector<vector<double>> nodesForQuadGrid(const vector<double> &X, const vector<double>& Y){
	cout << "\n nodesForQuadGrid\n";
	vector<vector<double>> nodes = vector<vector<double>>(X.size()*Y.size());
	
	for (int j = 0; j < int(Y.size()); j++){
		for (int k = 0; k < int(X.size()); k++){
			nodes[j*Y.size()+k] = { X[k],Y[j] };
		}
	}
	return nodes;
}

vector<vector<int>> edgesForQuadGrid(const int &nodeX, const int &nodeY){
	cout << "\n edgesForQuadGrid\n";
	vector<vector<int>> edges;
	for (int j = 0; j < nodeY; j++){
		for (int k = 0; k < nodeX - 1; k++){
			int n0 = nodeX * j + k;
			edges.push_back({ n0, n0 + 1 });// грани в плоскости паралл XY кроме 	
		}

		if (j != int(nodeY - 1)){
			for (int k = 0; k < nodeX; k++){
				int n0 = nodeX * j + k;
				//int n2 = nodeX*(j + 1) + k;
				edges.push_back({ n0,n0 + nodeX });
			}
		}
	}
	return edges;
}

vector<vector<int>> elementBaseNodes(const int &nodeX, const int &nodeY){
	cout << "\n elementBaseNodes\n";
	vector<vector<int>> elemNodes;

	for (int j = 0; j < nodeY - 1; j++){
			for (int k = 0; k < nodeX - 1; k++){
				int n0 = nodeX * j + k;
				int e0 = ((nodeX - 1) + nodeX)*j + k;
				//		cout << "e0 = " << e0 << endl;
				//	int e1 = e0 + (nodeX - 1);// e2 = e1+1
				//int e3 = e0 + nodeX + (nodeX - 1); 
				//ƒќЅј¬»“№ в описание ЁЋ≈ћ≈Ќ“ј грани и ребра
				elemNodes.push_back({ n0,	n0 + 1,	n0 + nodeX, n0 + nodeX + 1 });
			}
		}
	return elemNodes;
}

vector<vector<int>> elementBF(const int &elemRank, const int& NElmRows, const vector<vector<int>> &elemNodes) {
	
	cout << "\n elementBF\n";
	int NBFonElm = 0;
	if (elemRank == 1) { return elemNodes; }
	else {
		/*if (elemRank == 2) NBFonElm = 9;
		else if (elemRank == 3) NBFonElm = 16;*/
		vector<vector<int>> elemBF(elemNodes.size());
		int NElmCols = elemNodes.size() / NElmRows; // NElmCol
		int nodeX = NElmRows + 1, nodeY = NElmRows + 1;
		int nElm = 0;
		int ad = elemRank - 1;
		int n0, n3, n6, n4, n8, n12;
		for (int j = 0; j < NElmCols; j++) {
			for (int k = 0; k < NElmRows; k++) {
				if (elemRank == 2) {
					nElm = j * NElmCols + k;
					n0 = elemNodes[nElm][0] + (NElmCols * ad + NElmCols + 1 + ad * NElmCols)*j + k * ad;
					n3 = n0 + NElmCols + 1 + ad * NElmCols;
					n6 = n3 + NElmCols + 1 + ad * NElmCols;
					//		cout << "e0 = " << e0 << endl;
					//	int e1 = e0 + (nodeX - 1);// e2 = e1+1
					//int e3 = e0 + nodeX + (nodeX - 1); 
					//ƒќЅј¬»“№ в описание ЁЋ≈ћ≈Ќ“ј грани и ребра
					elemBF[nElm]={ n0,		n0 + 1,		n0 + 2,
										n3,		n3 + 1,		n3 + 2,
										n6,		n6 + 1,		n6 + 2 };
				}
				else if (elemRank == 3) {
						nElm = j * NElmCols + k;
						n0 = elemNodes[nElm][0] + (NElmCols * ad + NElmCols + 1 + ad * NElmCols)*j + k * ad;
						n4 = n0 + NElmCols + 1 + ad * NElmCols;
						n8 = n4 + NElmCols + 1 + ad * NElmCols;
						n12 = n8 + NElmCols + 1 + ad * NElmCols;
						//		cout << "e0 = " << e0 << endl;
						//	int e1 = e0 + (nodeX - 1);// e2 = e1+1
						//int e3 = e0 + nodeX + (nodeX - 1); 
						//ƒќЅј¬»“№ в описание ЁЋ≈ћ≈Ќ“ј грани и ребра
						elemBF[nElm] = { n0,		n0 + 1,		n0 + 2,		n0 + 3,
										n4,		n4 + 1,		n4 + 2,		n4 + 3,
										n8,		n8 + 1,		n8 + 2,		n8 + 3,
										n12,	n12 + 1,		n12 + 2,		n12 + 3 };
					}
				}

			}
		
		return elemBF;
	}
}

vector<vector<int>> elementEdges(const int &nodeX, const int &nodeY) {
	
	cout << "\n elementEdges \n ";
	vector<vector<int>> elemEdges;
	for (int j = 0; j < nodeY - 1; j++){
		for (int k = 0; k < nodeX - 1; k++){
			int n0 = nodeX * j + k;
			int e0 = ((nodeX - 1) + nodeX)*j + k;
			//		cout << "e0 = " << e0 << endl;
			//	int e1 = e0 + (nodeX - 1);// e2 = e1+1
			//int e3 = e0 + nodeX + (nodeX - 1); 
			//ƒќЅј¬»“№ в описание ЁЋ≈ћ≈Ќ“ј грани и ребра
			
			elemEdges.push_back({ e0, e0 + nodeX - 1, e0 + nodeX, e0 + 2 * nodeX - 1 });
		}
	}
	return elemEdges;
}

void elemsBndryList(const int &basisRank, const int &NnodeX, const int &NnodeY,
		const vector<vector<int>> &neighbours, const vector<vector<int>> elemsBF,
		vector<vector<vector<int>>> &bndrElems) {
/*  bndrElems: 
	индекс границы (0/1/2/3)  это номер вектора на верхнем уровне
	в 0-€чейке каждой строки будут хранитьс€ номер элемента 
	ƒалее - только номера тех Ѕ‘, кт к данной границе принадлежат*/
	
	cout << "\n elemsBndryList \n ";
	int Nbndr_elms = 2 * ((NnodeX - 1) + (NnodeY - 1)),
		Nbndr_elmsX = (NnodeX - 1), Nbndr_elmsY = (NnodeY - 1);
	
	//  число Ѕ‘ на одной стороне  Ё (дл€ пр€м-го эл-та 1го пор€дка - 2
	int NBFinLine;
	if (basisRank == 1) { NBFinLine = 2; }
	else if (basisRank == 2) { NBFinLine = 3; }
	else if (basisRank == 3) { NBFinLine = 4; }
	else { cout << " basis Rank more then maximum option or wrong value, basisRank = " << basisRank << endl; }

	if (bndrElems.size() != Nbndr_elms) {
		bndrElems = { vector<vector<int>>(Nbndr_elmsY), vector<vector<int>>(Nbndr_elmsY),
			vector<vector<int>>(Nbndr_elmsX), vector<vector<int>>(Nbndr_elmsX) };
	}
	//	cout << Nbndr_elms << endl;
	// счетчики текущего номера граничного эл-та дл€ каждой граничной линии
	int k1 = 0, k2 = 0, k3 = 0, k4 = 0;
	vector<int> BndryBFofElem = vector<int>(NBFinLine+1);
	for (int i = 0; i < int(neighbours.size()); i++) {
		if (neighbours[i][0] != 0) {
			cout << "i=" << i << ":   ";
			BndryBFofElem[0] = i;// нулевой элемент строки хранит глобальный номер граничного элемента
			if (neighbours[i][0] == 1 || neighbours[i][0] == 13 || neighbours[i][0] == 14) {
				for (size_t m = 0; m < NBFinLine; m++) {
					BndryBFofElem[m+1] = elemsBF[i][m*NBFinLine];
				}
				bndrElems[0][k1] = BndryBFofElem;
				//			cout << "BND 1,  k = "<< k1<<",  ";
				k1++;
			}

			if (neighbours[i][0] == 2 || neighbours[i][0] == 24 || neighbours[i][0] == 23) {
				BndryBFofElem[0] = i;// нулевой элемент строки хранит глобальный номер граничного элемента	
				for (size_t m = 0; m < NBFinLine; m++) {
					BndryBFofElem[m+1] = elemsBF[i][m*NBFinLine + (NBFinLine - 1)];
				}
				bndrElems[1][k2] = BndryBFofElem;
				//		cout << "BND 2,  k = " << k2 << ",  ";
				k2++;
			}
			if (neighbours[i][0] == 3 || neighbours[i][0] == 23 || neighbours[i][0] == 13) {
				for (size_t m = 0; m < NBFinLine; m++) {
					BndryBFofElem[m+1] = elemsBF[i][m];
				}
				bndrElems[2][k3] = BndryBFofElem;
				//	cout << "BND 3  k = " << k3 << ",  ";
				k3++;
			}
			if (neighbours[i][0] == 4 || neighbours[i][0] == 14 || neighbours[i][0] == 24) {
				for (size_t m = 0; m < NBFinLine; m++) {
					BndryBFofElem[m+1] = elemsBF[i][m + (NBFinLine*(NBFinLine - 1))];
				}
				bndrElems[3][k4] = BndryBFofElem;
				//	cout << "BND 4,  k = " << k4 << ",  ";
				k4++;
			}
			//	cout << endl;
		}
		//neighbours[i]
		//cout<<"we work on it\n";
	}
}

void CGbndryBF(const int &basisRank, const int& nX, const int& nY, const vector<vector<int>> &elemBF, vector<vector<int>> &bndrBF) {
	cout << "\n CGbndryBF \n ";
	int nBFinRow = (nX - 1)*(basisRank - 1),
		nBFinCol = (nY - 1)*(basisRank - 1);
	bndrBF = { vector<int>(nBFinRow), vector<int>(nBFinRow), vector<int>(nBFinCol), vector<int>(nBFinCol) };
	for (int i = 0; i < nBFinRow; i++) {
		bndrBF[0][i] = i;
		bndrBF[1][i] = i + (nBFinCol - 1)*nBFinRow;
	}

	for (int i = 0; i < nBFinCol; i++) {
		bndrBF[2][i] = i * nBFinRow;
		bndrBF[3][i] = (i + 1)*(nBFinRow)-1;
	}
}

void DGelemsBndryList(const int &Nelem, const int &NBFon_elm, const int &NnodeX, const int &NnodeY,
	const vector<vector<int>> &neighbours, vector<vector<int>> &bndrElems, vector<vector<int>> &bndryBF) {
	// Checking
	cout << "\n DGelemsBndryList \n ";
	int Nbndr_pnts = 2 * ((NnodeX - 1)*pow(NBFon_elm, 0.5) + (NnodeY - 1)*pow(NBFon_elm, 0.5)),
		Nbndr_elms = 2 * ((NnodeX - 1) + (NnodeY - 1)),
		Nbndr_elmsX = (NnodeX - 1), Nbndr_elmsY = (NnodeY - 1);

	if (bndryBF.size() != Nbndr_pnts) {
		bndryBF = vector<vector<int>>(Nbndr_pnts, vector<int>(2, 0));
	}

	if (bndrElems.size() != Nbndr_elms) {
		// в первых двух €чейках каждой строки будут хранитьс€ номер элемента и индекс границы (1/2/3/4) 
		// ƒалее - только номера тех Ѕ‘, кт на данной границе
		bndrElems = vector<vector<int>>(Nbndr_elms); // , vector<int>(2+pow(NBFon_elm, 0.5), 0)); 
	}
	int k1 = 0, k2 = 0, k3 = 0, k4 = 0;
	for (int i = 0, j = 0, k = 0; i < int(neighbours.size()); i++) {
		// bndrElems
		if (neighbours[i][0] != 0) {
			cout << "i=" << i << ":   ";
			if (neighbours[i][0] == 1 || neighbours[i][0] == 13 || neighbours[i][0] == 14) {
				bndrElems[k1] = { 1, i, i*NBFon_elm + 0, i*NBFon_elm + 2 };
				//			cout << "BND 1,  k = "<< k1<<",  ";
				k1++;
			}

			if (neighbours[i][0] == 2 || neighbours[i][0] == 24 || neighbours[i][0] == 23) {
				bndrElems[Nbndr_elmsY + k2] = { 2,i, i*NBFon_elm + 1, i*NBFon_elm + 3 };
				//		cout << "BND 2,  k = " << k2 << ",  ";
				k2++;
			}
			if (neighbours[i][0] == 3 || neighbours[i][0] == 23 || neighbours[i][0] == 13) {
				bndrElems[Nbndr_elmsX + Nbndr_elmsY + k3] = { 3,i, i*NBFon_elm + 0, i*NBFon_elm + 1 };
				//	cout << "BND 3  k = " << k3 << ",  ";
				k3++;
			}
			if (neighbours[i][0] == 4 || neighbours[i][0] == 14 || neighbours[i][0] == 24) {
				bndrElems[2 * Nbndr_elmsX + Nbndr_elmsY + k4] = { 4,i, i*NBFon_elm + 2, i*NBFon_elm + 3 };
				//	cout << "BND 4,  k = " << k4 << ",  ";
				k4++;
			}
		}
	}
}
