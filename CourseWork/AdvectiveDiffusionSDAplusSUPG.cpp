#include "pch.h"
#include "AdvectiveDiffusionSDAplusSUPG.h"


AdvecDifSDAplusSUPG::AdvecDifSDAplusSUPG(AbstractModel * newModel, Mesh * newMesh, SDAFactory * newDAfitches, const vector<double>& timeParams)
{
	T0 = timeParams[0];
	Tn = timeParams[1];
	NtimeSteps = timeParams[2];
	timeStep = (Tn - T0) / NtimeSteps;
	sdaFactory = newDAfitches;
	mesh = newMesh;
	model = newModel;
}

AdvecDifSDAplusSUPG::AdvecDifSDAplusSUPG(AbstractModel *newModel, Mesh* newMesh, SDAFactory *newDAfitches){
	model = newModel;
	mesh = newMesh;
	sdaFactory = newDAfitches;
	cout << "\n AdvectiveDiffusionSDAplusSUPG ::Constructor \n";
}

AdvecDifSDAplusSUPG::~AdvecDifSDAplusSUPG()
{
	cout << "\n AdvectiveDiffusionSDAplusSUPG ::Destructor \n";
}


void AdvecDifSDAplusSUPG::calcLocalMatrices(const currentParamsValue &Params, const quadrElem &curElem,
	vector<vector<double>> &resultA, vector<double> &resultB) {

	cout << "\n AdvectiveDiffusionSDAplusSUPG ::calcLocalMatrices \n";
	vector<vector<double>> localStiffMatrix, localAdvectiveMatrix, localSUPGMatrix;
	vector<double> localbVector, localSUPGVector;
	vector<vector<double>> FI = basis->getFiValuesInGaussPoints();
	vector<vector<vector<double>>> gradFI = basis->getGradFiValuesInGaussPoints();

	SolutionDiscreteAnalog::calcLocalAdvectiveMatrix(Params.velocity, curElem, FI, gradFI, localAdvectiveMatrix);
	SolutionDiscreteAnalog::calcLocalStiffnessMatrix(Params.lymbda, curElem, gradFI, localStiffMatrix);
	if(abs(vectorNorma(Params.velocity))>=10e-10){
		if (basis->getbasisRank() == 1) {
			calcLocalSUPGMatrix1stOrderBasis(0, Params, curElem, gradFI, localSUPGMatrix, localSUPGVector);

		}
	else if (basis->getbasisRank() > 1 && basis->getbasisRank() < 3) {
		vector<vector<vector<double>>> laplacFi = basis->getLaplacFiValuesInGaussPoints();
		calcLocalSUPGMatrix2ndOrderBasis(0, Params, curElem, gradFI, laplacFi, localSUPGMatrix, localSUPGVector);
	}
	}
	SolutionDiscreteAnalog::calcLocalRightVector(curElem, FI, localbVector);
	
	resultA = localStiffMatrix + localAdvectiveMatrix + localSUPGMatrix;
	resultB = localbVector+localSUPGVector;
	if (curElem.x1 == 0 && curElem.y1 == 0) {
		coutVector(localStiffMatrix, "localStiffMatrix");
		coutVector(localAdvectiveMatrix, "localAdvectiveMatrix");
		coutVector(localStiffMatrix + localAdvectiveMatrix, "G+K");
		coutVector(localSUPGMatrix, "localSUPGMatrix");
		coutVector(resultA, "localA");
		coutVector(resultB, "localB", 'H');
	}
}

void AdvecDifSDAplusSUPG::calcLocalSUPGMatrix1stOrderBasis(const int& nTime, const currentParamsValue &Params, 
	const quadrElem &Elem,	const vector<vector<vector<double>>> &gradFi,
	vector<vector<double>> &resultSUPGmatrix, vector<double> &resultSUPGvector){
	cout << "\n AdvectiveDiffusionSDAplusSUPG ::calcLocalSUPGMatrix1stOrderBasis \n";

	vector<double> bGrad(gradFi.size());
	vector<vector<double>> xyPairs;
	basis->getXYforGaussPoints(Elem, xyPairs);
	vector<double> sourceOnElm = model->getSourceValues(xyPairs);

	vector<vector<double>> gradGrad = vector<vector<double>>(gradFi.size(),	vector<double>(gradFi.size(), 0.));

	vector<double> gaussWghts = basis->getGaussWeights();
	double tau = calcSUPGtau(Elem, Params, 6.);


	for (size_t i = 0; i < basis->getbasisSize(); i++) {
		bGrad[i] = summVectorElms(gaussWghts*sourceOnElm*(gradFi[i][0] * (Elem.J[0][0] * Params.velocity[0]) +
			gradFi[i][1] * (Elem.J[1][1] * Params.velocity[1])))*(1. / Elem.detJ);
		}
	resultSUPGvector = bGrad * tau;

	if (nTime <= 0) {
	for (size_t i = 0; i < basis->getbasisSize(); i++) {
			for (size_t j = 0; j < basis->getbasisSize(); j++) {
				// std::cout << " j" << j << "]\n";
				gradGrad[i][j] =
					summVectorElms(gradFi[i][0] * gradFi[j][0] * gaussWghts)*(pow((Elem.J[0][0] * Params.velocity[0]), 2)*(1. / Elem.detJ))
					+ summVectorElms(gradFi[i][1] * gradFi[j][1] * gaussWghts)*(pow((Elem.J[1][1] * Params.velocity[1]), 2)* (1. / Elem.detJ))
					+ summVectorElms((gradFi[i][1] * gradFi[j][0] + gradFi[i][0] * gradFi[j][1])* gaussWghts)
					* (Elem.J[1][1] * Elem.J[0][0] * (1. / Elem.detJ) * Params.velocity[1] * Params.velocity[0]);
			}
		}
	resultSUPGmatrix = gradGrad * tau;
	}
		

}

void AdvecDifSDAplusSUPG::calcLocalSUPGMatrix2ndOrderBasis(const int& nTime, const currentParamsValue &Params, const quadrElem &Elem,
		const vector<vector<vector<double>>> &gradFi, const vector<vector<vector<double>>> &laplacFi,
	vector<vector<double>> &resultSUPGmatrix, vector<double> &resultSUPGvector){

	cout << "\n AdvectiveDiffusionSDAplusSUPG ::calcLocalMatrices \n";
	// добавить якобиан

	vector<vector<double>> gradGrad(gradFi.size(), vector<double>(gradFi.size(), 0.)),
		laplaclaplac(gradFi.size(), vector<double>(gradFi.size(), 0.)),
		laplacGrad(gradFi.size(), vector<double>(gradFi.size(), 0.)),
		gradlaplac(gradFi.size(), vector<double>(gradFi.size(), 0.));
	
	vector<double> gaussWghts = basis->getGaussWeights(), 
		bGrad(gradFi.size()), blaplac(gradFi.size());
	vector<vector<double>> xyPairs;
	basis->getXYforGaussPoints(Elem, xyPairs);

	vector<double> sourceOnElm = model->getSourceValues(xyPairs);

	for (int i = 0; i < gradFi.size(); i++) {
		for (int j = 0; j < gradFi.size(); j++) {
			std::cout << " j" << j << "]\n";
			laplaclaplac[i][j] = (
				summVectorElms(laplacFi[i][0] * laplacFi[j][0] * gaussWghts)*pow(Elem.J[0][0], 4) +
				summVectorElms(laplacFi[i][0] * laplacFi[j][1] * gaussWghts)*pow(Elem.J[0][0], 2)*pow(Elem.J[1][1],2) +
				summVectorElms(laplacFi[i][1] * laplacFi[j][0] * gaussWghts)*pow(Elem.J[0][0], 2)*pow(Elem.J[1][1], 2)+
				summVectorElms(laplacFi[i][1] * laplacFi[j][1] * gaussWghts)*pow(Elem.J[1][1], 4))*Params.lymbda; // lyambda может и зависеть от координат
			
			laplacGrad[i][j] =(
				summVectorElms(laplacFi[i][0] * gradFi[j][0] * gaussWghts)*(Params.velocity[0] * pow(Elem.J[0][0], 3))
				+ summVectorElms(laplacFi[i][0] * gradFi[j][1] * gaussWghts)*(Params.velocity[1] * pow(Elem.J[0][0], 2)*Elem.J[1][1])
				+ summVectorElms(laplacFi[i][1] * gradFi[j][1] * gaussWghts)*(Params.velocity[1] * pow(Elem.J[1][1], 3))
				+ summVectorElms(laplacFi[i][1] * gradFi[j][0] * gaussWghts)*(Params.velocity[0] * pow(Elem.J[1][1], 2)*Elem.J[0][0]))*Params.lymbda;

			gradlaplac[i][j] = (
				summVectorElms(laplacFi[j][0] * gradFi[i][0] * gaussWghts)*(Params.velocity[0] * pow(Elem.J[0][0], 3))
				+ summVectorElms(laplacFi[j][0] * gradFi[i][1] * gaussWghts)*(Params.velocity[1] * pow(Elem.J[0][0], 2)*Elem.J[1][1])
				+ summVectorElms(laplacFi[j][1] * gradFi[i][1] * gaussWghts)*(Params.velocity[1] * pow(Elem.J[1][1], 3))
				+ summVectorElms(laplacFi[j][1] * gradFi[i][0] * gaussWghts)*(Params.velocity[0] * pow(Elem.J[1][1], 2)*Elem.J[0][0]))*Params.lymbda;


			gradGrad[i][j] = summVectorElms(gradFi[i][0] * gradFi[j][0] * gaussWghts)*pow(Elem.J[0][0]*Params.velocity[0],2)
						+ summVectorElms(gradFi[i][1] * gradFi[j][1] * gaussWghts)*pow(Elem.J[1][1]*Params.velocity[1], 2)
						+ summVectorElms((gradFi[i][1] * gradFi[j][0] + gradFi[i][0] * gradFi[j][1]) 
							* gaussWghts)*Elem.J[1][1]*Elem.J[0][0]*Params.velocity[1]* Params.velocity[0];
		}
		bGrad[i] = summVectorElms(gaussWghts*sourceOnElm*(gradFi[i][0]*(Elem.J[0][0]*Params.velocity[0]) +
									gradFi[i][1]*(Elem.J[1][1]*Params.velocity[1])));
		blaplac[i]= summVectorElms(gaussWghts*sourceOnElm*(laplacFi[i][0]*(Elem.J[0][0]*Elem.J[0][0]*Params.lymbda) +
									laplacFi[i][1]*(Elem.J[1][1]*Elem.J[1][1]*Params.lymbda)));

	}
	double dim = sqrt(pow(Elem.hx, 2) + pow(Elem.hy, 2)),
		tau = 1.,
		Pecle = vectorNorma(Params.velocity)*dim / (6.0 * Params.lymbda);

	if (Pecle < 1.0) { tau = (dim*Pecle) / (2.0 * vectorNorma(Params.velocity)); }
	else { tau = dim / (2.0 * vectorNorma(Params.velocity)); }
	// добавить якобиан
	resultSUPGmatrix= (laplaclaplac - laplacGrad + gradlaplac - gradGrad)*((1. / Elem.detJ)*tau);
	resultSUPGvector = (blaplac + bGrad)*(tau*(1. / Elem.detJ));
	
}

void AdvecDifSDAplusSUPG::calcRightVectorNewTimeStep(const vector<double>& prevStepU)
{
	slae.b = vector<double>(slae.b.size(), 0.);
	cout << "\n SolutionDiscreteAnalog: calcRightVectorNewTimeStep \n";
	int Nelem = mesh->elemBF.size(); // slae.iptr.size(); // или jptr??????????
	quadrElem curElem;
	currentParamsValue params;
	vector<double>localVectorB, localVectorIC;
	vector<vector<double>> localMatrixA, localSUPGMatrix;
	vector<double> localIC, localSUPGVector;
	vector<vector<double>> FI = basis->getFiValuesInGaussPoints();
	vector<vector<vector<double>>> gradFI = basis->getGradFiValuesInGaussPoints();
	cout << "Nelem = " << Nelem << endl;
	for (int i = 0; i < Nelem; i++) {
		curElem = calcLocalElem(i); // X, Y достанет из mesh
		cout << "elem " << i << ":: ";
		params = model->getCurParamsValue(curElem.middleX, curElem.middleY);
		calcLocalRightVector(curElem, FI, localVectorB);
		vector<vector<double>> _xy;
		basis->getXYforGaussPoints(curElem, _xy);
		localIC = calcPointSolutionOnElmInGP(i);
		calcLocalTimeVector(curElem, localIC, FI, localVectorIC);
		if (abs(vectorNorma(params.velocity)) >= 10e-10) {
			if (basis->getbasisRank() == 1) {
				calcLocalSUPGMatrix1stOrderBasis(0, params, curElem, gradFI, localSUPGMatrix, localSUPGVector);

			}
			else if (basis->getbasisRank() > 1 && basis->getbasisRank() < 3) {
				vector<vector<vector<double>>> laplacFi = basis->getLaplacFiValuesInGaussPoints();
				calcLocalSUPGMatrix2ndOrderBasis(0, params, curElem, gradFI, laplacFi, localSUPGMatrix, localSUPGVector);
			}
		}
		
		for (int j = 0; j < 4; j++) {
			slae.b[mesh->elemBF[i][j]] += (localVectorB[j] + localVectorIC[j]+localSUPGVector[j]);
		}
	}
}

double AdvecDifSDAplusSUPG::calcSUPGtau(quadrElem Elem, currentParamsValue Params, double Ctau) {
	double dim = sqrt(pow(Elem.hx, 2) + pow(Elem.hy, 2)),
		tau = 1., timeParam = 1,
		Pecle = vectorNorma(Params.velocity)*dim / (2.0 * Params.lymbda);
		
	if (timeStep > 1e-20) timeParam = timeStep;

	cout << "tanh(Pecle)= " << tanh(Pecle) << ", 1/tanh(Pecle)= " << 1. / tanh(Pecle) << ",   (1. / Pecle)) = " << (1. / Pecle) << endl;

	tau = timeParam* Ctau*(dim / (2. * vectorNorma(Params.velocity)))*((1. / tanh(Pecle)) - (1. / Pecle));
	
	cout << "Pecle = " << Pecle << ",  tau = " << tau;
	return tau;
}
