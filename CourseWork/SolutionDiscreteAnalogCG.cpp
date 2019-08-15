#include "pch.h"
#include "SolutionDiscreteAnalogCG.h"


SolutionDiscreteAnalogCG::SolutionDiscreteAnalogCG()
{
	cout << "\nSolutionDiscreteAnalogCG::default Constructor()\n";
}

SolutionDiscreteAnalogCG::SolutionDiscreteAnalogCG(AbstractModel * newModel, Mesh * newMesh, SDAFactory * newDAfitches, const vector<double>& timeParams)
{
	T0 = timeParams[0];
	Tn = timeParams[1];
	NtimeSteps = timeParams[2];
	timeStep = (Tn - T0) / NtimeSteps;
	sdaFactory = newDAfitches;
	mesh = newMesh;
	model = newModel;
}

SolutionDiscreteAnalogCG::SolutionDiscreteAnalogCG(AbstractModel * newModel, Mesh * newMesh, SDAFactory * newDAfitches)
{
	cout << "/SolutionDiscreteAnalogCG:  Constructor/n";
	sdaFactory = newDAfitches;
	mesh = newMesh;
	model = newModel;
}

SolutionDiscreteAnalogCG::~SolutionDiscreteAnalogCG()
{
	cout << "\nSolutionDiscreteAnalogCG::Destructor()\n";
}

void SolutionDiscreteAnalogCG::defineSolutionProcess(const string &basisItem, const int& rankBasis)
{
	cout << "\nSolutionDiscreteAnalogCG::defineSolutionProcess()\n";	
		slae = sdaFactory->createSlaePortret(*mesh);
		basis = sdaFactory->createBasis(basisItem, rankBasis);
}

void SolutionDiscreteAnalogCG::assembleSlae(const int& nTime){
	cout << "\nSolutionDiscreteAnalogCG: assembleSlae\n";
	int Nelem = mesh->elemBF.size(); // slae.iptr.size(); // или jptr??????????
	quadrElem curElem;
	currentParamsValue curModelParams;
	vector<double>localVectorB, localTimeVector, localIC;
	vector<vector<double>> localMatrixA, localTimeMatrix;
	cout << "Nelem = "<<Nelem << endl;
	vector<vector<double>> _xy;
	vector<vector<double>> FI = basis->getFiValuesInGaussPoints();
	int stStr, endStr, nelem;
	if (nTime > 0) slae.b = vector<double>(slae.b.size(), 0.);
	for (int i = 0; i < Nelem; i++) {
		curElem = calcLocalElem(i); // X, Y достанет из mesh
		cout << "elem " << i << ":: ";
	//	curElem.coutElem();

		curModelParams = model->getCurParamsValue(curElem.middleX, curElem.middleY);
		
		calcLocalMatrices(curModelParams, curElem, localMatrixA, localVectorB);
		if (nTime >= 0) {	
			basis->getXYforGaussPoints(curElem, _xy);
			localIC = model->getInitCondValues(_xy);
			calcLocalTimeMatrix(curElem, FI, localTimeMatrix);
			calcLocalTimeVector(curElem, localIC, FI, localTimeVector);
			localVectorB = timeStep*localVectorB + localTimeVector;
			localMatrixA = timeStep*localMatrixA + localTimeMatrix;
		}
	
		for (int j = 0; j < 4; j++) {
			slae.b[mesh->elemBF[i][j]] += localVectorB[j];			
				stStr = slae.iptr[mesh->elemBF[i][j]];  // начало строки
				endStr = slae.iptr[mesh->elemBF[i][j] + 1]; // конец
				for (int k = 0; k < 4; k++) {
					nelem = findElement(slae.jptr, stStr, endStr, mesh->elemBF[i][k]);// Elem[i][k] OR Elem[i][j] ORElem[j][k]  ???????????????????
					slae.A[nelem] += localMatrixA[j][k];
			}
		}
	}
	
}

void SolutionDiscreteAnalogCG::BoundCondDerichlet(const int &nBoundary, const int&nTimeStep) {
	// boundary condition is some constant function nondepended by coordinates
	// nBoundary = {1, 2, 3, 4}
	cout << "\nSDA CG::BoundCondDerichlet() on boundary "<< nBoundary <<"\n";
	int nDiagElm, curBF, nCurElm = 0, stStr, endStr, NBFinLine = mesh->bndrElms[0][0].size()-1;
	vector<double> bndryCoord;
	quadrElem curElem;
	
	for (size_t i = 0, n = mesh->bndrElms[nBoundary - 1].size(); i < n; i++) {
		
		nCurElm = mesh->bndrElms[nBoundary-1][i][0];
		curElem = calcLocalElem(nCurElm);
		cout << "bndry elem = " << nCurElm << ":: ";
	//	curElem.coutElem();
		if (nBoundary == 1 || nBoundary == 2) {
			basis->getXorYforCarrier('Y', curElem, bndryCoord);		// лева€ или права€ границы
		}

		if (nBoundary == 3 || nBoundary == 4) {
			basis->getXorYforCarrier('X', curElem, bndryCoord);		// нижн€€ или верхн€€ границы
		}
	//	coutVector(bndryCoord, "bndryCoord", 'H');
		for (size_t k = 0; k < NBFinLine; k++) {
			curBF = mesh->bndrElms[nBoundary-1][i][k+1]; // текуща€ строка в —Ћј”
			stStr = slae.iptr[curBF];
			endStr = slae.iptr[curBF + 1];
			//SLAU[p] = help;
			//cout << "start = " << st_str << "  end = " << end_str << endl;
			if (nTimeStep <= 0) {
				for (int j = stStr; j < endStr; j++) {
					slae.A[j] = 0.0;
				}
				nDiagElm = findElement(slae.jptr, stStr, endStr, curBF);
				//	cout << "p = " << p << "  int(p) = " << int(p) << "  nelem = " << nelem << endl;
				slae.A[nDiagElm] = 1.;
			}
			if (nBoundary == 1) slae.b[curBF] = model->BC1(curElem.x1, bndryCoord[k]);
			else if (nBoundary == 2) slae.b[curBF] = model->BC1(curElem.x2, bndryCoord[k]);
			else if (nBoundary == 3) slae.b[curBF] = model->BC1(bndryCoord[k], curElem.y1);
			else if (nBoundary == 4) slae.b[curBF] = model->BC1(bndryCoord[k], curElem.y2);

		//	slae.coutSLAE("slae BC1");
		}
	}
}


