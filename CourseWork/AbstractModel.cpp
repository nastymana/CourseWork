#include"pch.h"
#include "AbstractModel.h"


AbstractModel::AbstractModel()
{// cout<<"\n  :: \n";
	cout << "\n AbstractModel :: default Constructor\n";
}


AbstractModel::AbstractModel(const vector<double>& inData, const vector<double>& domainParams)
{
	cout << "\n AbstractModel :: abstract Constructor\n";
	
}

AbstractModel::~AbstractModel()
{
	cout << "\n AbstractModel :: Destructor\n";
}

currentParamsValue AbstractModel::getCurParamsValue(const double & x, const double & y)
{
	cout << "\n AbstractModel :: getCurParamsValues\n";
	
	return {	lymbda(x,y),
				gamma(x,y), 
				velocity(x,y) };
}

double AbstractModel::getDomain(const double & x, const double & y)
{
	cout << "\n AbstractModel :: getDomain\n";
	return 0.0;
}

vector<double> AbstractModel::getSourceValues(const vector<vector<double>>& xyPairs)
{
	cout << "\n AbstractModel :: getSourceValues\n";
	vector<double> SV(xyPairs.size());
	for (size_t i = 0; i < xyPairs.size(); i++){
	
		SV[i] = source(xyPairs[i][0], xyPairs[i][1]);
	//	coutVector(xyPairs[i], "xy", 'H');
		cout << SV[i] << endl;
	}
	return SV;
}
vector<double> AbstractModel::getBC1Values(const vector<vector<double>>& xyPairs)
{
	cout << "\n AbstractModel :: getSourceValues\n";
	vector<double> SV(xyPairs.size());
	for (size_t i = 0; i < xyPairs.size(); i++) {

		SV[i] = BC1(xyPairs[i][0], xyPairs[i][1]);
		//	coutVector(xyPairs[i], "xy", 'H');
		cout << SV[i] << endl;
	}
	return SV;
}

vector<double> AbstractModel::getInitCondValues(const vector<vector<double>>& xy)
{
	cout << "\n AbstractModel :: getSourceValues\n";
	vector<double> SV(xy.size());
	for (size_t i = 0; i < xy.size(); i++) {

		SV[i] = initialCond(xy[i][0], xy[i][1]);
	}
	return SV;
}

void AbstractModel::addDescription(const string & newLine)
{
	cout << "\n AbstractModel :: addDescription\n";
	discription.push_back(newLine);
}

vector<string> AbstractModel::getDescription()
{ 
	return discription;
}

void AbstractModel::getRealSol2D(const vector<double>& X, const vector<double>& Y, vector<vector<double>> Sol)
{
	cout << "\n AbstractModel :: getRealSol2D\n";
	Sol = vector<vector<double>>(Y.size(), vector<double>(X.size()));
	for (size_t i = 0; i < Y.size(); i++){
		for (size_t j = 0; j < X.size(); j++){
			Sol[i][j] = realSol(X[j], Y[i]);
		}
	}
}
void  AbstractModel::getRealSol(const vector<vector<double>>& XYpairs, vector<double>&Sol)
{
	cout << "\n AbstractModel :: getRealSol\n";
	Sol = vector<double>(XYpairs.size());
	for (size_t i = 0; i < XYpairs.size(); i++) {
			Sol[i] = realSol(XYpairs[i][0], XYpairs[i][1]);
		}
}

void AbstractModel::getRealSol(const vector<double>& X, const vector<double>& Y, vector<double>&Sol)
{
	cout << "\n AbstractModel :: getRealSol\n";
	Sol = vector<double>(X.size()*Y.size());
	for (size_t i = 0; i < X.size(); i++) {
		for (size_t j = 0; j < Y.size(); j++) {
			Sol[i] = realSol(X[j], Y[i]);
		}
	}

}

void AbstractModel::coutModelParams()
{
	// coutVector( domain, "domain",'H');
	// vector<string> discription;
	//vector<double> domainData;
	cout << "lyambdaParam = " << lyambdaParam << ", gammaParam = " << gammaParam << ",\n velocityParams = {"
		<< velocityParam1 << ", " << velocityParam2 << "}, sourceParam = " << sourceParam << endl;
}
