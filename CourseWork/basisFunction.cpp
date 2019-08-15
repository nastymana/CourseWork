#include "pch.h"
#include "basisFunction.h"


basisFunction::basisFunction(){
	cout << "\n basisFunction:: default constructor\n";
}

basisFunction::basisFunction(const int & rankItem){
	cout << "\n basisFunction:: Constructor\n";
	basisRank = 0;
	GaussRank = 0;
}

basisFunction::~basisFunction(){ cout << "\n basisFunction:: Destructor\n"; }

// для треугольников возможно придется переопределить эту функцию
vector<vector<double>> basisFunction::getFiValuesInGaussPoints(){
	cout << "\n basisFunction:: getFiValuesInGaussPoints\n";
	vector<vector<double>> Fi(basisSize, vector<double>(pow(GaussRank,2), 0.));
	
	for (size_t s = 0; s < basisSize; s++){
		for (size_t i = 0; i < GaussRank; i++){
			for (size_t j = 0; j < GaussRank; j++) {
				Fi[s][i*GaussRank+j] = fi[s](GaussPoints[j], GaussPoints[i]);
			}
		}	
	}
	return Fi;
}

vector<vector<double>> basisFunction::getFiValuesInGaussPoints(int nBndry)
{
vector<vector<double>> Fi(basisSize, vector<double>(GaussRank, 0.));

	for (size_t s = 0; s < basisSize; s++){
		if (nBndry == 1 || nBndry == 2) {
			for (size_t i = 0; i < GaussRank; i++) {
				Fi[s][i] = fi[s](carrierX[nBndry-1], GaussPoints[i]);
			}
		}
		else if (nBndry == 3 || nBndry == 4) {
			for (size_t i = 0; i < GaussRank; i++) {
				Fi[s][i] = fi[s](GaussPoints[i], carrierY[nBndry - 3]);
			}
		}
	}
	return Fi;

}

vector<double> basisFunction::getFiValuesInPoint(const double &eps, const double &eta) {
	cout << "\n basisFunction:: getFiValuesInGaussPoints\n";
	vector<double> Fi(basisSize);

	for (size_t s = 0; s < basisSize; s++) {
				Fi[s] = fi[s](eps, eta);
	}
	return Fi;
}

vector<vector<vector<double>>> basisFunction::getGradFiValuesInGaussPoints(){
	cout << "\n basisFunction:: getGradFiValuesInGaussPoints\n";
	vector<vector<vector<double>>> gradFi(basisSize, vector<vector<double>>(2, vector<double>(pow(GaussRank, 2), 0.)));
	vector<double> gradi;
	for (size_t s = 0; s < basisSize; s++) {
		for (size_t i = 0; i < GaussRank; i++) {
			for (size_t j = 0; j < GaussRank; j++) {
				gradi = grad[s](GaussPoints[j], GaussPoints[i]);
				gradFi[s][0][i*GaussRank + j] = gradi[0] ;
				gradFi[s][1][i*GaussRank + j] = gradi[1];
			}
		}
	}
	return gradFi;
}

vector<vector<vector<double>>> basisFunction::getGradFiValuesInGaussPoints(int nBndry)
{
	vector<vector<vector<double>>> gradFi(basisSize, vector<vector<double>>(2, vector<double>(GaussRank, 0.)));
	vector<double> gradi;

	
	for (size_t s = 0; s < basisSize; s++) {
		if (nBndry == 1 || nBndry == 2) {
			for (size_t i = 0; i < GaussRank; i++) {
				gradi = grad[s](carrierX[nBndry - 1], GaussPoints[i]);
				gradFi[s][0][i] = gradi[0];
				gradFi[s][1][i] = gradi[1];
			}
		}
		else if (nBndry == 3 || nBndry == 4) {
			for (size_t i = 0; i < GaussRank; i++) {
				gradi = grad[s](GaussPoints[i], carrierY[nBndry - 3]);
				gradFi[s][0][i] = gradi[0];
				gradFi[s][1][i] = gradi[1];
			}
		}
	}
	return gradFi;

}

vector<vector<vector<double>>> basisFunction::getLaplacFiValuesInGaussPoints()
{
	cout << "\n basisFunction:: getLaplacFiValuesInGaussPoints\n";
	vector<vector<vector<double>>> laplacFi(basisSize, vector<vector<double>>(2, vector<double>(pow(GaussRank, 2), 0.)));
	vector<double> laplacian_i;
	for (size_t s = 0; s < basisSize; s++) {
		for (size_t i = 0; i < GaussRank; i++) {
			for (size_t j = 0; j < GaussRank; j++) {
				laplacian_i= laplacian[s](GaussPoints[j], GaussPoints[i]);
				laplacFi[s][0][i*GaussRank + j] = laplacian_i [0];
				laplacFi[s][1][i*GaussRank + j] = laplacian_i[1];
			}
		}
	}
	return laplacFi;
}

vector<double> basisFunction::getGaussWeights(){
	cout << "\n basisFunction:: getGaussWeights\n";
	vector<double> GPquadrElem(pow(GaussRank, 2));
	for (size_t i = 0; i < GaussRank; i++) { // по У
		for (size_t j = 0; j < GaussRank; j++) {// по Х
			GPquadrElem[i*GaussRank + j] = (GaussWeights[j] * GaussWeights[i]);
		};
		return GPquadrElem;
	}
}

vector<double> basisFunction::getGaussWeights(const int &nBndry) {
	cout << "\n basisFunction:: getGaussWeights\n";
		return GaussWeights;
	
}

vector<double> basisFunction::getCarrier(const char & XorYitem)
{
	if (XorYitem == 'X' || XorYitem == 'x') { return carrierX; }
	else if (XorYitem == 'Y' || XorYitem == 'y') { return carrierY; }
	else cout << " \n WRONG input parameter \n";
}

vector<vector<double>> basisFunction::getPairsOfGaussPoints(){
	cout << "\n basisFunction:: getPairsOfGaussPoints\n";
	vector<vector<double>> GPpairs(pow((GaussRank+1), 2), vector<double>(2));
	for (size_t i = 0; i < GaussRank; i++){
		for (size_t j = 0; j < GaussRank; j++){
			GPpairs[i*GaussRank + j] = { GaussPoints[j], GaussPoints[i] };
		}
	}
	return GPpairs;
}

vector<double> basisFunction::getGaussPoints(const char & XorYitem)
{
	return GaussPoints;
}


int basisFunction::getbasisSize()
{
	return basisSize;
}

int basisFunction::getbasisRank()
{
	return basisRank;
}

int basisFunction::getGaussPointNumber()
{
	return pow(GaussPoints.size(), 2);
}

void basisFunction::getXYforGaussPoints(const quadrElem &elem, vector<vector<double>>& XYpairs){
	cout << "\n basisFunction:: getXYforGaussPoints\n";
	XYpairs = vector<vector<double>>(pow(GaussRank, 2), vector<double>(2));
	for (int i = 0; i < GaussRank; i++) {
		for (size_t j = 0; j < GaussRank; j++){
			XYpairs[i*GaussRank + j][0] = epsToX(elem.x1, elem.x2, GaussPoints[j]);// Х
			XYpairs[i*GaussRank + j][1] = epsToX(elem.y1, elem.y2, GaussPoints[i]); // У
		}
	}
}

void basisFunction::getXorYforCarrier(const char &XorYitem, const quadrElem & elem, vector<double>& coord)
{
	coord = vector<double>(carrierX.size());
	if (XorYitem == 'X' || XorYitem == 'x') {
		coord = vector<double>(carrierX.size());
		for (size_t i = 0; i < carrierX.size(); i++) {
			coord[i] = elem.x1 + elem.hx*(carrierX[i] - carrierX[0]);
		}
	}
	else if (XorYitem == 'Y' || XorYitem == 'y') {
		for (size_t i = 0; i < carrierY.size(); i++) {
			coord[i] = elem.y1 + elem.hy*(carrierY[i] - carrierY[0]);
		}
	}
}

void basisFunction::getXorYforGaussPoints(const char &XorYitem, const quadrElem & elem, vector<double>& coord)
{
	coord = vector<double>(GaussPoints.size());
	if (XorYitem == 'X' || XorYitem == 'x') { 
		for (size_t j = 0; j < GaussRank; j++) {
			coord[j] = epsToX(elem.x1, elem.x2, GaussPoints[j]);	// Х
		}
	}
	else if (XorYitem == 'Y' || XorYitem == 'y') {
		for (size_t j = 0; j < GaussRank; j++) {
			coord[j] = epsToX(elem.y1, elem.y2, GaussPoints[j]);
		}
	}
	else cout << " \n WRONG input parameter \n";
}

void basisFunction::getXYforGaussPoints(const int &nBnd, const quadrElem & elem, vector<vector<double>>& coord)
{
	coord = vector<vector<double>>(GaussPoints.size(), vector<double>(2));
	if (nBnd == 3 || nBnd == 4) {
		for (size_t j = 0; j < GaussRank; j++) {
			coord[j][0] = epsToX(elem.x1, elem.x2, GaussPoints[j]);	// Х
			if (nBnd == 3)	coord[j][1] = elem.y1;
			else coord[j][1] = elem.y2;
		}
	}
	else if (nBnd == 1 || nBnd== 2) {
		for (size_t j = 0; j < GaussRank; j++) {
			if (nBnd == 1)	coord[j][0] = elem.x1;
			else coord[j][0] = elem.x2;
			coord[j][1] = epsToX(elem.y1, elem.y2, GaussPoints[j]);
		}
	}
	else cout << " \n WRONG input parameter \n";
}

double basisFunction::getDetJacobian(const vector<double> &x, const vector<double> &y) {
	// по умолчанию на мастер - элементе[-1, 1]x[-1, 1]
	cout << "\n basisFunction::getDetJacobian \n";
	// double hx_m = 2., hy_m = 2.;
	double detJ = ((carrierX.back() - carrierX[0])*(carrierY.back() - carrierY[0]))/ ((x[1] - x[0]) * (y[1] - y[0]));
	return detJ;
}

vector<vector<double>>  basisFunction::getJacobian(const vector<double>& x, const vector<double>& y){
	// по умолчанию на мастер - элементе[-1, 1]x[-1, 1] для ПРЯМОУГОЛЬНОГО эл-та 
	cout << "\n basisFunction:: getJacobian\n";
	/* 
		dKsi/dx,		0
			0,		dEta/dy	
	*/
	return {
		{((carrierX.back() - carrierX[0]) / (x[1] - x[0])),						0.							  },
		{				0.,								((carrierY.back() - carrierY[0]) / (y[1] - y[0])) } };
}

double epsToX(const double &x1, const double &x2, const double &eps) {
	// eps_to_x_m11 по умолчанию на мастер-элементе [-1,1]x[-1,1]
	cout << "\n basisFunction:: epsToX\n";
	return  ((x2 + x1) *0.5 + eps * (x2 - x1) *0.5);
}

// видимо для произвольного 4хУгольника
double epsilon_to_x(double x1, double x2, double x3, double x4, double eps, double eta){
	//cout << "x1 = " << x1 << "  x2 = " << x2 << "  x3 = " << x3 << "  x4 = " << x4 << endl;
	return  (1 - eps)*(1 - eta)*x1 + eps * (1 - eta)*x2 + (1 - eps)*eta*x3 + eps * eta*x4;
	
}

