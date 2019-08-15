#pragma once
void createClassicSolutionforCG_MP1(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams);

void createDGSolutionforCG_MP1(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams);

void createClassicSolutionforCGplusSUPG_MP1(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams);
void createClassicSolutionforMP2(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams);
void createClassicSolutionforCG_MP4(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams);
void createClassicSolutionforCGplusSUPG_MP4(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams);

void createCGSolutionforCG_MP7(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams);
void createClassicSolutionforCG_MP8(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams);

void createClassicSolutionforCG_RP1(const int Ntest, const vector<double>&indataGrid,
	const vector<double> &indataParams, const vector<double> &timeParams);

void outputfSolutionResults(const wstring path, const wstring newCatalogName,
	const vector<double> &analytSol, const vector<double> &solution, 
	const vector<double> &indataMesh, const vector<double> &meshX, const vector<double> &meshY,
	const vector<double> modelParams, const vector<string>& discription,
	const vector<double> &normaR1, const vector<double> &normaRes, const vector<double> &normaExRes);