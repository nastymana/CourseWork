#pragma once

class AbstractModel
{
protected:
	vector<double> domain;
	vector<string> discription;
	vector<double> domainData;
	double lyambdaParam;
	double gammaParam;
	double velocityParam1, velocityParam2;
	double sourceParam;
	double timeParam;
public:
	AbstractModel();
	AbstractModel(const vector<double> &inData, const vector<double>& domainParams);
	virtual	~AbstractModel();
	virtual double realSol(const double &x, const double &y) = 0;
	virtual double BC1(const double &x, const double &y) = 0;
	virtual double BC2(const double &x, const double &y) = 0;
	virtual double BC3(const double &x, const double &y) = 0;
	virtual double initialCond(const double &x, const double &y) = 0;

	virtual double source(const double &x, const double &y) = 0;

	virtual double lymbda(const double &x, const double &y) = 0;
	virtual double gamma(const double &x, const double &y) = 0;
	virtual vector<double> velocity(const double &x, const double &y) = 0;

	
	
	currentParamsValue getCurParamsValue(const double &x, const double &y);
	double getDomain(const double &x, const double &y);
	vector<double> getSourceValues(const vector<vector<double>>& xy);
	vector<double> getInitCondValues(const vector<vector<double>>& xy);
	vector<double> getBC1Values(const vector<vector<double>>& xyPairs);
	void addDescription(const string &newLine);
	vector<string> getDescription();
	void getRealSol(const vector<vector<double>> &XYpairs, vector<double>&result);
	void getRealSol(const vector<double>& X, const vector<double>& Y, vector<double>& result);
	void getRealSol2D(const vector<double> &X, const vector<double> &Y, vector<vector<double>> SOl);
	void coutModelParams();

// protected:

	

};

