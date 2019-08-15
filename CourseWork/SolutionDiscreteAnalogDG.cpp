#include "pch.h"
#include "SolutionDiscreteAnalogDG.h"


SolutionDiscreteAnalogDG::SolutionDiscreteAnalogDG()
{
	cout << "\nSolutionDiscreteAnalogDG::default Constructor()\n";
}

SolutionDiscreteAnalogDG::SolutionDiscreteAnalogDG(AbstractModel * newModel, Mesh * newMesh, SDAFactory * newDAfitches, const vector<double>& timeParams)
{
}

SolutionDiscreteAnalogDG::SolutionDiscreteAnalogDG(AbstractModel * newModel, Mesh * newMesh, SDAFactory * newDAfitches)
{
	cout << "/SolutionDiscreteAnalogDG:  Constructor/n";
	sdaFactory = newDAfitches;
	mesh = newMesh;
	model = newModel;

}

SolutionDiscreteAnalogDG::~SolutionDiscreteAnalogDG()
{
	cout << "\nSolutionDiscreteAnalogDG::Destructor()\n";
}

void SolutionDiscreteAnalogDG::defineSolutionProcess(const string & basisItem, const int & rankBasis)
{
	cout << "\nSolutionDiscreteAnalogDG::defineSolutionProcess()\n";
	basis = sdaFactory->createBasis(basisItem, rankBasis);
	slae = sdaFactory->createSlaePortret(*mesh);
	//slae.coutSLAE("slae initial");
	

}

void SolutionDiscreteAnalogDG::assembleSlae(const int &nTime){
	
	cout << "\nSolutionDiscreteAnalogDG::assembleSlae()\n";
	cout << "\nSolutionDiscreteAnalogCG: assembleSlae\n";
	int Nelem = mesh->elemBF.size(); // slae.iptr.size(); // или jptr??????????
	double scaleMu = 0.;
	quadrElem curElem;
	currentParamsValue curModelParams;
	vector<double>localVectorB, localTimeVector, localIC;
	vector<vector<double>> localMatrixA, localTimeMatrix;
	cout << "Nelem = " << Nelem << endl;
	vector<vector<double>> _xy;
	vector<vector<double>> FI = basis->getFiValuesInGaussPoints();
	int stStr, endStr, nelem, NBF;
	if (nTime > 0) slae.b = vector<double>(slae.b.size(), 0.);
	for (int i = 0; i < Nelem; i++) {
		curElem = calcLocalElem(i); // X, Y достанет из mesh
		cout << "elem " << i << ":: ";
		//	curElem.coutElem();
		NBF = basis->getbasisSize();
		curModelParams = model->getCurParamsValue(curElem.middleX, curElem.middleY);

		calcLocalMatrices(curModelParams, curElem, localMatrixA, localVectorB);
		coutVector(localMatrixA, "localMatrixA");
		coutVector(localVectorB, "localVectorB", 'H');
		if (nTime >= 0) {
			basis->getXYforGaussPoints(curElem, _xy);
			localIC = model->getInitCondValues(_xy);
			calcLocalTimeMatrix(curElem, FI, localTimeMatrix);
			//		coutVector(localTimeMatrix, "localTimeMatrix");
			calcLocalTimeVector(curElem, localIC, FI, localTimeVector);
			//		coutVector(localTimeVector, "localTimeVector", 'H');
			localVectorB = timeStep * localVectorB + localTimeVector;
			localMatrixA = timeStep * localMatrixA + localTimeMatrix;
		}
		//	coutVector(localMatrixA, "localMatrixA");

		for (int j = 0; j < NBF; j++) {
			slae.b[mesh->elemBF[i][j]] += localVectorB[j];
			stStr = slae.iptr[mesh->elemBF[i][j]];  // начало строки
			endStr = slae.iptr[mesh->elemBF[i][j] + 1]; // конец
			for (int k = 0; k < NBF; k++) {
				nelem = findElement(slae.jptr, stStr, endStr, mesh->elemBF[i][k]);// Elem[i][k] OR Elem[i][j] ORElem[j][k]  ???????????????????
				slae.A[nelem] += localMatrixA[j][k];
			}
		}
	}
	
	int nedge_elem = 4, rank_basis = basis->getbasisRank();
	double normaV = pow(innerProduct(curModelParams.velocity, curModelParams.velocity ),2);

	scaleMu = 0.01;//(nedge_elem*rank_basis*gamma*lymbda*normaV) / 1.0;
	//slae.coutSLAE(" just diag elements");

		addNumFluxToSLAE();
		//slae.coutSLAE(" added num fluxes");
}

void SolutionDiscreteAnalogDG::BoundCondDerichlet(const int &nBndry, const int &nTimeStep)
{
	cout << "\nSolutionDiscreteAnalogDG::BoundCondDerichlet()\n";
	int NBF = basis->getbasisSize();
	vector<vector<double>> localCondMatrix;
	vector<vector<double>> FI = basis->getFiValuesInGaussPoints(nBndry);
	vector<vector<vector<double>>> gradFi = basis->getGradFiValuesInGaussPoints(nBndry);
	
	// boundary condition is some functin depended by coordinates
	int stx = 0, sty,
	nElem, curBF,
		nCurElm = 0, stStr, endStr,
		NBFinLine = mesh->bndrElms[0][0].size() - 1;
	quadrElem elem;
	currentParamsValue params;
	vector<double> bndryCoor, localCondVector;
	vector<double> nleft = { -1., 0. }, nright = { 1., 0. }, ndoun = { 0., -1. }, nup = { 0., 1. }, n;
	int curBFplace;
	for (int i = 0; i < mesh->bndrElms[nBndry - 1].size(); i++) {
		nElem = mesh->bndrElms[nBndry-1][i][0];
		cout << "nelem " << nElem << endl;
		elem = calcLocalElem(nElem);
		params = model->getCurParamsValue(elem.middleX, elem.middleY);
		calcLocalBcDerichletMatrices(nBndry, FI, gradFi, elem, params, localCondMatrix, localCondVector);		
		
		coutVector(localCondVector, "local Derich;let vector", 'H');
		coutVector(localCondMatrix, "local Derichlet Matrix");
		for (size_t k = 0; k < NBF; k++) {
			curBF = mesh->elemBF[nElem][k];//mesh->bndrElms[nBndry - 1][i][k + 1]; // текущая строка в СЛАУ
			
			
			stStr = slae.iptr[curBF];
			endStr = slae.iptr[curBF + 1];
			curBFplace = findElement(slae.jptr, stStr, endStr, curBF-k); // or Elems[i]+j ???
			//SLAU[p] = help;
			cout << "current BF (string) " << curBF <<", cur BF place"<<curBFplace<< endl;
			//cout << "start = " << st_str << "  end = " << end_str << endl;
			//if (nTimeStep <= 0) {		
			slae.b[curBF] += localCondVector[k];
			for (int j = 0; j < NBF; j++) {
				slae.A[curBFplace+j] += localCondMatrix[k][j];
				
			}
		}
	}
	
}

void SolutionDiscreteAnalogDG::calcLocalBcDerichletMatrices(const int nBndry,
	const vector<vector<double>>& Fi, const vector<vector<vector<double>>>& gradFi,
	quadrElem & elem, currentParamsValue & params, vector<vector<double>>& localAd, vector<double>& localBd)
{
	// Вывести всю процедуру заполнения векторов в отдельную функцию, в кт на вход подается номер кейса
	
	int NBF = basis->getbasisSize();
	localAd = vector<vector<double>>(NBF, vector<double>(NBF));
	localBd = vector<double>(NBF);
	vector<vector<double>> G_nx (NBF, vector<double>(NBF)),
		G_ny (NBF, vector<double>(NBF)), G_IP (NBF, vector<double>(NBF));
	
	vector<double> b_ny = vector<double>(NBF), b_nx = vector<double>(NBF), b_IP = vector<double>(NBF);
	vector<double> gaussWghts = basis->getGaussWeights(nBndry);
	vector<vector<double>> _xy; 
	double revDetJ;
	vector<double> Jdotn = innerProduct(elem.J, elem.nVectors[nBndry - 1]);

	if (nBndry == 1 || nBndry == 2) { revDetJ = 1. / elem.J[1][1]; }// левая или правая границы	
	
	if (nBndry == 3 || nBndry == 4) { revDetJ = 1. / elem.J[0][0]; }// нижняя или верхняя границы
	
	basis->getXYforGaussPoints(nBndry, elem, _xy);
	double mu = IPstabParamMU(nBndry ,elem, params);
	coutVector(_xy, "xy");
	vector<double> BCDerValues = model->getBC1Values(_xy);
	coutVector(BCDerValues, "BC Derichlet Values ", 'H');

	for (size_t i = 0; i < NBF; i++) {
		for (size_t j = 0; j < NBF; j++) {
			G_nx[i][j] = summVectorElms((-1.*gradFi[i][0] * Fi[j] - gradFi[j][0]*Fi[i])* gaussWghts); //(JdotV[0] * (1. / curElem.detJ));				// Aconss [ij]K1.K1)
			G_ny[i][j] = summVectorElms((-1.*gradFi[i][1] *Fi[j] - gradFi[j][1] * Fi[i])*gaussWghts);
			G_IP[i][j] = summVectorElms((Fi[i] * Fi[j])*gaussWghts);
		}
		localBd[i] = summVectorElms(BCDerValues* ((mu*revDetJ)*Fi[i] 
			- gradFi[i][0] *(revDetJ*Jdotn[0]) 
			- gradFi[i][1] *(revDetJ*Jdotn[1]))*gaussWghts);
	}
	localAd = (revDetJ* params.lymbda)*(G_nx *Jdotn[0]+ G_ny*Jdotn[1] + G_IP);
	localBd = params.lymbda*localBd;
}

void SolutionDiscreteAnalogDG::addNumFluxToSLAE(){
	
	vector<double> nleft = { -1., 0. }, nright = { 1., 0. }, ndoun = { 0., -1. }, nup = { 0., 1. };
	
	int nin, curBFleft=0, curBFright= 0, curBFdoun= 0, curBFup = 0,
		nexleft = 0, nexright = 0, nexdoun = 0, nexup = 0,
		nElinRow = 0, NBF = basis->getbasisSize();
	vector<vector<double>>		
		Eex_left(NBF, vector<double>(NBF)),	Ein_left(NBF, vector<double>(NBF)),	
		Eex_right(NBF, vector<double>(NBF)), Ein_right(NBF, vector<double>(NBF)), 
		Eex_doun(NBF, vector<double>(NBF)),	Ein_doun(NBF, vector<double>(NBF)),	
		Eex_up(NBF, vector<double>(NBF)), Ein_up(NBF, vector<double>(NBF));
	vector<vector<double>> J;
	double detJ = 0, mu;
	quadrElem elem;
	currentParamsValue params;
	for (int i = 0; i < mesh->elemBF.size(); i++) {
		elem = calcLocalElem(i);//detJacobian(hx, hy, 2);// type_ME = 2 - ME [-1, 1]x[-1, 1] type_ME=1 [0,1][0,1]

		nin = mesh->elemBF[i][0];
		params = model->getCurParamsValue(elem.middleX, elem.middleY);
		double ch_l = 0., ch_r = 0., ch_u = 0., ch_d = 0.;

		if (mesh->nghbrs[i][1] != -1) {
			caclLocalNumerFluxMatrices(1, elem, params, Ein_left, Eex_left);
			curBFleft = mesh->elemBF[mesh->nghbrs[i][1]][0];
			ch_l = 1.;
		}
		if (mesh->nghbrs[i][2] != -1){
			caclLocalNumerFluxMatrices(2, elem, params, Ein_right, Eex_right);
			curBFright = mesh->elemBF[mesh->nghbrs[i][2]][0];
			ch_r = 1.;
		}
		if (mesh->nghbrs[i][3] != -1) {
			caclLocalNumerFluxMatrices(3, elem, params, Ein_doun, Eex_doun);
			curBFdoun = mesh->elemBF[mesh->nghbrs[i][3]][0];
			ch_d = 1.;
		}
		if (mesh->nghbrs[i][4] != -1) {
			caclLocalNumerFluxMatrices(4, elem, params, Ein_up, Eex_up);
			curBFup = mesh->elemBF[mesh->nghbrs[i][4]][0];
			ch_u = 1.;
		}
	
		if (i == 0) {
			coutVector(Eex_left, "Eexleft");
			coutVector(Ein_left, "Ein_left");
			coutVector(Eex_right, "Eexright");
			coutVector(Ein_right, "Ein_right");
			coutVector(Eex_doun, "Eex_doun");
			coutVector(Ein_doun, "Ein_doun");
			coutVector(Eex_up, "Eex_up");
			coutVector(Ein_up, "Ein_up");
		}
		nexleft = 0, nexright = 0, nexdoun = 0, nexup = 0;
		NBF = basis->getbasisSize();
			for (int j = 0; j < NBF; j++) {
				int st_str = slae.iptr[mesh->elemBF[i][j]],  // начало строки
					end_str = slae.iptr[mesh->elemBF[i][j] + 1];
				cout << "nin = " << nin+j << ", nexLeft=" << curBFleft +j<< ", nexRight=" 
					<< curBFright+j << " nexdoun=" << curBFdoun+j << ", nexUp=" << curBFup+j << endl;
				nin = findElement(slae.jptr, st_str, end_str, mesh->elemBF[i][0]); // or Elems[i]+j ???
				nexleft = findElement(slae.jptr, st_str, end_str, curBFleft);
				nexdoun = findElement(slae.jptr, st_str, end_str, curBFdoun);
				nexright = findElement(slae.jptr, st_str, end_str, curBFright);
				nexup = findElement(slae.jptr, st_str, end_str, curBFup);
				cout << "nin = " << nin << ", nexLeft=" << nexleft << " nexdoun=" << nexdoun <<", nexRight=" << nexright <<  ", nexUp=" << nexup << endl;
				for (int k = 0; k <NBF; k++) {// сделать проверку сколько граничных элементов есть
					slae.A[nexleft + k] += ch_l * Eex_left[j][k];
					slae.A[nexright + k] += ch_r * Eex_right[j][k];
					slae.A[nexdoun + k] += ch_d * Eex_doun[j][k];
					slae.A[nexup + k] += ch_u * Eex_up[j][k];
					slae.A[nin + k] += ch_l * Ein_left[j][k] + ch_r * Ein_right[j][k] +
						ch_u * Ein_up[j][k] + ch_d * Ein_doun[j][k];
				}
			}
		//	slae.coutSLAE("added num fluxes for elem " + to_string(i));
	}
}

void SolutionDiscreteAnalogDG::caclLocalNumerFluxMatrices(int nBnd, const quadrElem &elem,
	const currentParamsValue& params,
	vector<vector<double>> &resultEin, vector<vector<double>> &resultEex) {

	// надо объеденить в один массив размера NBFx2NBF Eex, Ein - получим два массива E(nx), E(ny) и
	// или три + E(nz)
	char item = 'x';
	std::vector<double> coord2 = basis->getGaussPoints(item);
	// интегрирования вдоль левой границы по оси Y
	int NGauss = basis->getbasisRank(), NBF = basis->getbasisSize();		// dimension
	int nBndNgh = 0;
	if (nBnd == 1) nBndNgh = 2;
	else if (nBnd == 2) nBndNgh = 1;
	else if (nBnd == 3) nBndNgh = 4;
	else if (nBnd == 4) nBndNgh = 3;
	
	vector<vector<double>> FI_K1 = basis->getFiValuesInGaussPoints(nBnd),
		FI_K2 = basis->getFiValuesInGaussPoints(nBndNgh);
	coutVector(FI_K1, "Fi in");
	vector<vector<vector<double>>> gradFI_K1 = basis->getGradFiValuesInGaussPoints(nBnd),
		gradFI_K2 = basis->getGradFiValuesInGaussPoints(nBndNgh);
	
	// Вывести всю процедуру заполнения векторов в отдельную функцию, в кт на вход подается номер кейса
	vector<vector<double>> EinK1_nx = vector<vector<double>>(NBF, vector<double>(NBF)),
		EinK1_ny = vector<vector<double>>(NBF, vector<double>(NBF)),
		EexK2_K1_nx = vector<vector<double>>(NBF, vector<double>(NBF)),
		EexK2_K1_ny = vector<vector<double>>(NBF, vector<double>(NBF)),
		EexIP = vector<vector<double>>(NBF, vector<double>(NBF)),
		EinIP = vector<vector<double>>(NBF, vector<double>(NBF));

	vector<double> gaussWghts = basis->getGaussWeights(nBnd);
	double mu = IPstabParamMU(nBnd, elem, params), revDetJ;

	if (nBnd == 1 || nBnd== 2) { revDetJ = 1. / elem.J[1][1]; }// левая или правая границы	

	if (nBnd == 3 || nBnd == 4) { revDetJ = 1. / elem.J[0][0]; } //
	for (size_t i = 0; i < NBF; i++) {
		for (size_t j = 0; j < NBF; j++) {
			
			EinK1_nx[i][j] = summVectorElms(FI_K1[j] * (-1.)*gradFI_K1[i][0]		// Asymm[ij]K1.K1 
				+ FI_K1[i] * (-1.)*(gradFI_K1[j][0])* gaussWghts); //(JdotV[0] * (1. / curElem.detJ));				// Aconss [ij]K1.K1)
			
			EinK1_ny [i][j]= summVectorElms(FI_K1[j] * (-1.)*gradFI_K1[i][1]   // 
				+ (-1.)*FI_K1[i] * gradFI_K1[j][1] * gaussWghts);
			EexK2_K1_nx[i][j] = summVectorElms(FI_K2[j] * gradFI_K1[i][0]	// Asymm[ij] K1.K2 
				+ FI_K1[i] * (-1.)*gradFI_K2[j][0]*gaussWghts);			// Acons[ij] K1.K2 
			EexK2_K1_ny[i][j] = summVectorElms(FI_K2[j] * gradFI_K1[i][1]
				+ FI_K1[i] * (-1.)*gradFI_K2[j][1]*gaussWghts);// т.к. она будет затем умножаться также на нормаль внешнюю для рассматриваемого эл-та
		
			EinIP[i][j] = summVectorElms(FI_K1[i] * FI_K1[j]*gaussWghts); 		// A_IP	K1.K1																				// а она равна (-n) внешнего элемента
			EexIP[i][j] = summVectorElms(-1.*FI_K1[i] * FI_K2[j] * gaussWghts); 	// A_IP	K1.K2																						// а она равна (-n) внешнего элемента
		}
	}	

	vector<double> Jdot_n = innerProduct(elem.J, elem.nVectors[nBnd-1]);
	resultEin = (EinK1_nx) * (Jdot_n[0] *revDetJ*params.lymbda*0.5) +
		(EinK1_ny) * (Jdot_n[1] * params.lymbda*0.5 *revDetJ) + (EinIP) * (mu *revDetJ*params.lymbda);

	resultEex = EexK2_K1_nx * (Jdot_n[0] * revDetJ*params.lymbda*0.5) +
		EexK2_K1_ny * (Jdot_n[1] * params.lymbda*0.5*revDetJ) +
		EexIP * (mu *revDetJ*params.lymbda);
	std::cout << "Good Buy!";

}


double SolutionDiscreteAnalogDG::IPstabParamMU(const int nBnd, const quadrElem &elem, const currentParamsValue& params) {

	double mu = 4.;
	if (nBnd == 1 || nBnd == 2) mu = mu / elem.hy;
	else if (nBnd == 3 || nBnd == 4) mu = mu / elem.hx;
	return mu;
}

