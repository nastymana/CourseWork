#include "pch.h"
#include "SolutionDiscreteAnalog.h"





SolutionDiscreteAnalog::SolutionDiscreteAnalog()
{
	cout << "/nSolutionDiscreteAnalog: default Constructor/n";
}

SolutionDiscreteAnalog::SolutionDiscreteAnalog(AbstractModel * newModel, Mesh * newMesh, SDAFactory * newDAfitches, const vector<double>& timeParams)
{
}

SolutionDiscreteAnalog::SolutionDiscreteAnalog(AbstractModel *newModel, Mesh *newMesh, SDAFactory *newDAfitches)
{
	cout << "/nSolutionDiscreteAnalog:  Constructor/n";
	sdaFactory = newDAfitches;
	mesh = newMesh;
	model = newModel;
}


SolutionDiscreteAnalog::~SolutionDiscreteAnalog() {
	cout << "/nSolutionDiscreteAnalog: Destructor/n";
}

SLAE_CSR* SolutionDiscreteAnalog::getSLAE()
{
	cout << "\n SolutionDiscreteAnalog::getSLAE  \n ";
	return &slae;
}

Mesh * SolutionDiscreteAnalog::getMesh()
{
	cout << "\n SolutionDiscreteAnalog::getMesh  \n ";
	return mesh;
}

AbstractModel* SolutionDiscreteAnalog::getModel()
{
	cout << "\n SolutionDiscreteAnalog::getModel  \n ";
	return model;
}

void SolutionDiscreteAnalog::calcLocalMassMatrix(const double & curGamma, const quadrElem & curElem,
	const vector<vector<double>> &FI, vector<vector<double>> &resultM)
{
	cout << "/nSolutionDiscreteAnalog: calcLocalMassMatrix/n";
	size_t NBF = basis->getbasisSize();
	// if M.size !=NBF
	resultM = vector<vector<double>>(NBF, std::vector<double>(NBF,0.));

	vector<double> gaussWghts = basis->getGaussWeights();
	for (size_t i = 0; i < NBF; i++){
		for (size_t j = 0; j < NBF; j++){
			resultM[i][j] = summVectorElms(FI[j] * FI[i] * gaussWghts)*((1./curElem.detJ)*curGamma);
		}
	}
}

void SolutionDiscreteAnalog::calcLocalStiffnessMatrix(const double &curLambda, const quadrElem &curElem,
	const vector<vector<vector<double>>> &gradFI, vector<vector<double>> &resultG)
{
	cout << "/nSolutionDiscreteAnalog: calcLocalStiffnessMatrix/n";
	size_t NBF = basis->getbasisSize();
	vector<vector<double>>  Gx(NBF, std::vector<double>(NBF, 0.)), Gy(NBF, std::vector<double>(NBF, 0.));

	vector<double> gaussWghts = basis->getGaussWeights();
	for (size_t i = 0; i < NBF; i++) {
		for (size_t j = 0; j < NBF; j++) {
			Gx[i][j] = summVectorElms(gradFI[j][0] * gradFI[i][0] * gaussWghts)*((curLambda / curElem.detJ)*pow(curElem.J[0][0], 2));
			Gy[i][j] = summVectorElms(gradFI[j][1] * gradFI[i][1] * gaussWghts)* ((curLambda / curElem.detJ)*pow(curElem.J[1][1], 2));
		}
	}
//	coutVector(Gx, "Gx");
//	coutVector(Gy, "Gy");
	resultG = Gx + Gy;
}

void SolutionDiscreteAnalog::calcLocalAdvectiveMatrix(const vector<double>& curVelocity, const quadrElem & curElem,
	const vector<vector<double>> &FI, const vector<vector<vector<double>>> &gradFI, vector<vector<double>> &resultK){
	cout << "/nSolutionDiscreteAnalog: calcLocalAdvectiveMatrix/n";
	size_t NBF = basis->getbasisSize();
	vector<vector<double>> Kx(NBF, std::vector<double>(NBF,0.)), Ky(NBF, std::vector<double>(NBF, 0.));
	vector<double> JdotV = innerProduct(curElem.J, curVelocity);
	vector<double> gaussWghts = basis->getGaussWeights();
	for (size_t i = 0; i < NBF; i++){
		for (size_t j = 0; j < NBF; j++){
			Kx[i][j] = summVectorElms(FI[j] * gradFI[i][0] * gaussWghts)*(JdotV[0] * (1. / curElem.detJ));
			Ky[i][j] = summVectorElms(FI[j] * gradFI[i][1] * gaussWghts)*(JdotV[1] * (1. / curElem.detJ));
		}
	}
resultK = Kx  + Ky  ;
}

void SolutionDiscreteAnalog::calcLocalRightVector(const quadrElem & curElem, 
	const vector<vector<double>> &FI,  vector<double>& resultB){
	cout << "/nSolutionDiscreteAnalog: calcLocalRightVector/n";
	// для общего случая, когда базис неоднородный по области. нужно сдлеать так, 
	// чтобы все функции для basis_function принимали на вход элемент и по координатам определеяли какой базис вернуть
	size_t NBF = basis->getbasisSize();
	vector<vector<double>> _xy;
	resultB = vector<double>(NBF, 0.);
	vector<double> w = basis->getGaussWeights();
	
	basis->getXYforGaussPoints(curElem, _xy);
	//coutVector(_xy, "_xy");
	vector<double> sourceValues = model->getSourceValues(_xy);
	 coutVector(sourceValues, "sourceValues ", 'H');
	for (size_t i = 0; i < NBF; i++) {
		resultB[i] = summVectorElms(w*sourceValues*FI[i])*(1. / curElem.detJ);
	}
	//coutVector(resultB, "b", 'H');
	
}

void SolutionDiscreteAnalog::calcLocalTimeMatrix( 
	const quadrElem & curElem, const vector<vector<double>>& FI,
	vector<vector<double>>& resultM)
{
	calcLocalMassMatrix(1., curElem, FI, resultM);
}
vector<double> SolutionDiscreteAnalog::calcPointSolutionOnElmInGP(int nelem){
	cout << "\n SolutionDiscreteAnalog::calcPointSolutionOnElmInGP  \n ";
	vector<double> basisWeights(basis->getbasisSize());
	vector<vector<double>> fiInGP = basis->getFiValuesInGaussPoints(); // fiInGP = matrix(basisSize*GP.size)
	fiInGP = transposeMatrix(fiInGP); // fiInGP = matrix(GP.size*basisSize)
	vector<double> solution(fiInGP.size());
//	cout << "basis->getbasisSize() = " << basis->getbasisSize() << endl;
	for (size_t i = 0; i < basis->getbasisSize(); i++) {
		cout << "mesh->elemBF[nelem][i] = " << mesh->elemBF[nelem][i] << endl;
	//	coutVector(mesh->elemBF[nelem], "elemBF[nelem]", 'H');
		basisWeights[i] = slae.x[mesh->elemBF[nelem][i]];
	}
//	coutVector(fiInGP, "fiInGP");
//	coutVector(basisWeights, "basisWeights", 'H');

	for (size_t i = 0; i < fiInGP.size(); i++){
		solution[i] = summVectorElms(basisWeights*fiInGP[i]);
	}

	return solution;
}

void SolutionDiscreteAnalog::calcRightVectorNewTimeStep(const vector<double>& prevStepU)
{
	slae.b = vector<double>(slae.b.size(), 0.);
		cout << "\n SolutionDiscreteAnalog: calcRightVectorNewTimeStep \n";
		int Nelem = mesh->elemBF.size(); // slae.iptr.size(); // или jptr??????????
		quadrElem curElem;
		vector<double>localVectorB, localVectorIC;
		vector<vector<double>> localMatrixA;
		vector<double> localPrevSol;
		vector<vector<double>> FI = basis->getFiValuesInGaussPoints();
		cout << "Nelem = " << Nelem << endl;
		for (int i = 0; i < Nelem; i++) {
			curElem = calcLocalElem(i); // X, Y достанет из mesh
			cout << "elem " << i << ":: ";

			calcLocalRightVector(curElem, FI, localVectorB);
			vector<vector<double>> _xy;
			basis->getXYforGaussPoints(curElem, _xy);
			localPrevSol = calcPointSolutionOnElmInGP(i);
			calcLocalTimeVector(curElem, localPrevSol, FI, localVectorIC);
			coutVector(localVectorB, "localVectorB", 'H');
			coutVector(localPrevSol, "localIC", 'H');
			coutVector(localVectorIC, "localVectorIC", 'H');

			
			for (int j = 0; j < 4; j++) {
				slae.b[mesh->elemBF[i][j]] += localVectorB[j]+ localVectorIC[j];
		}
	}
}



quadrElem SolutionDiscreteAnalog::calcLocalElem(int nElem){
	cout << "\n SolutionDiscreteAnalog::calcLocalElem  \n ";
	int stx = (nElem % (mesh->X.size() - 1)),
	sty = int(nElem / (mesh->X.size() - 1));
	double x1 = mesh->X[stx], x2 = mesh->X[stx + 1],
		y1 = mesh->Y[sty], y2 = mesh->Y[sty + 1];
	double hx = abs(x2 - x1),
		hy = abs(y2 - y1);
	double middleX = x1 + hx/2. , middleY= y1 + hy / 2.;
	double detJ = basis->getDetJacobian({ x1, x2 }, { y1, y2 });
	vector<vector<double>> J = basis->getJacobian({ x1, x2 }, { y1, y2 });
	vector<vector<double>> n = { {-1., 0.},{ 1., 0. },{0., -1.},{0., 1.} };

	return { x1, x2, y1, y2, hx, hy, middleX, middleY, detJ, J, n };
}

void SolutionDiscreteAnalog::coutSDA()
{
	cout << "\n SolutionDiscreteAnalog::coutSDA  \n ";
	model->coutModelParams();
	mesh->coutMesh();
	slae.coutSLAE("slae");
	basis->coutBasisParams();
	//SDAFactory* sdaFactory;
}

double SolutionDiscreteAnalog::L2Norma(){
	cout << "\n SolutionDiscreteAnalog::L2Norma  \n ";
	// (J(|U-Uh|^2 dW)^0.5
	// норму будем считать как сумму интегралов Гаусса по элементам
	vector<vector<double>> GP = basis->getPairsOfGaussPoints(), XYonGPinElm;
	vector<double> GW = basis->getGaussWeights(),
	 realSolOnElm(GP.size(), 0.), findedSolOnElm(GP.size(), 0.);
	quadrElem curElem;
	double sumDeltaU_Uh = 0, L2norma; // сумма разностей решения на 
	vector<double> squaredUdif;
	for (size_t i = 0; i < mesh->elemBF.size(); i++){
		// считаем параметры текущего КЭ
		curElem = calcLocalElem(i);
		cout << "elem:: " << i << endl;
		curElem.coutElem();
		// для текущего КЭ получаем значения координат в точках Гаусса
		basis->getXYforGaussPoints(curElem, XYonGPinElm);
//		coutVector(XYonGPinElm, "XY");
		// получаем вектор значений аналитического решения в точках Гаусса на реальном элементе
		model->getRealSol(XYonGPinElm, realSolOnElm);
	//	coutVector(realSolOnElm, "realSolOnElm", 'H');
		findedSolOnElm = calcPointSolutionOnElmInGP(i);
	//	coutVector(findedSolOnElm, "findedSolOnElm", 'H');
		degreeVector((realSolOnElm - findedSolOnElm), squaredUdif);
	//	coutVector(squaredUdif, "squaredUdif", 'H');
		// для текущего значения получаем вектор значений решения в точках Гаусса
		
		sumDeltaU_Uh+= summVectorElms(squaredUdif*GW);
	}
	L2norma = pow(sumDeltaU_Uh, 0.5);

	return L2norma;
}

void SolutionDiscreteAnalog::calcLocalTimeVector(const quadrElem & curElem, const vector<double>& curF, const vector<vector<double>>& FI, vector<double>& resultB){
	size_t NBF = basis->getbasisSize();
	vector<vector<double>> _xy;
	resultB = vector<double>(NBF, 0.);
	vector<double> w = basis->getGaussWeights();

	for (size_t i = 0; i < NBF; i++) {
		resultB[i] = summVectorElms(w*curF*FI[i])*(1. / curElem.detJ);
	}

}

int SolutionDiscreteAnalog::findPointInROw(const double &p, const vector<double> &row) {
	int count = 0;
	for (size_t i = 0; i < row.size()-1; i++)
	{
		if (p >= row[i] && p <= row[i + 1]) {
			count = i;
			i = row.size();
		}
	}
	return count;
}

double SolutionDiscreteAnalog::solutionInPoint(const double& x, const double& y) {
	/*double x0 = mesh->X[0], xn = mesh->X.back(),
	y0 = mesh->Y[0], yn = mesh->Y.back();*/
	int nx = findPointInROw(x, mesh->X),
		ny = findPointInROw(y, mesh->Y),
		nElem = nx+ ny*(mesh->X.size()-1);
	cout << "nx=" << nx<<", X="<<x << ", ny = " << ny << ", nelm = " << nElem << endl;
	quadrElem curElem = calcLocalElem(nElem);
	double eps = basis->XtoEps(mesh->X[nx], mesh->X[nx + 1], x),
		eta = basis->XtoEps(mesh->Y[ny], mesh->Y[ny + 1], y);;
	vector<double> basisWeights(basis->getbasisSize());
	cout << "eps=" << eps << ", eta=" << eta << endl;
	vector<double> fiInP = basis->getFiValuesInPoint(eps,eta); // fiInGP = matrix(basisSize*GP.size)
	coutVector(fiInP, "fi", 'H');
	double solution; 
	for (size_t i = 0; i < basis->getbasisSize(); i++) {
		cout << "mesh->elemBF[nelem][i] = " << mesh->elemBF[nElem][i] << endl;
		basisWeights[i] = slae.x[mesh->elemBF[nElem][i]];
	}
	coutVector(basisWeights, "basisWeights", 'H');
		solution = summVectorElms(basisWeights*fiInP);
		cout <<  "solution=" << solution << endl;
	return solution;

}

void SolutionDiscreteAnalog::calcPointSolution(const vector<double>&X, const vector<double>&Y, vector<double> &pointSol)
{
	// hQ   делитель для шага сетки
	cout << "\n SolutionDiscreteAnalog::calcPointSolution  \n ";
	cout << "\n we work on it \n";
	
	int nNodeX = X.size(),
		nNodeY = Y.size() ;
	pointSol = vector<double>(nNodeX*nNodeY);
	for (size_t i = 0; i < nNodeY; i++){
		for (size_t j = 0; j < nNodeX; j++){
			pointSol[i*nNodeX+j] = solutionInPoint(X[j], Y[i]);
			cout << "pointSol[" << i << "][" << j << "]="<<pointSol[i*nNodeX + j] << endl;
		}
	}
	
}

void SolutionDiscreteAnalog::addDescriptionToModel(const string & newLine)
{
	model->addDescription(newLine);
}

vector<string> SolutionDiscreteAnalog::getDescription()
{
	return model->getDescription();
}


