#pragma once

class Solver
{
public:
	Solver();
	virtual ~Solver();
	vector<double> solution;
	vector<double> res;
	vector<double> err;
	double stop_eps;
	double pivot;

};


// ************************************************************************************************************************************************************************
// ************************************************************* BiCgSTAB solvers family **********************************************************************************
// ************************************************************************************************************************************************************************

// Origin BiCGSTAB for dense SLAU matrix
vector<double> BiCGSTAB_dense(const vector<vector<double>> &A, const vector<double> &b, const vector<double> &newx0, const vector<double> &RealSol,
	int &niter, double &eps, int &stopres, vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit);

// Origin BiCGSTAB for CSR SLAU matrix
vector<double> BiCGSTAB_CSR(const vector<double> &A, const vector<double> &b, const vector<double> &newx0, const vector<double> &RealSol,
	const vector<int> &iptr, const vector<int> &jptr, int &niter, double &eps, int &stopres,
	vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit);


void BiCGSTAB_CSR(SLAE_CSR* slae, const vector<double> &newx0, const vector<double> &RealSol,
	int &niter, double &eps, int &pivot,
	vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit);

// Origin BiCGSTAB for CSR SLAU matrix parallized with OpenMP
vector<double> BiCGSTAB_CSR_OMP(const vector<double> &A, const vector<double> &b, const vector<double> &newx0, const vector<double> &RealSol,
	const vector<int> &iptr, const vector<int> &jptr, int &niter, double &eps, int &stopres,
	vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit, const int threadNum);

// Origin BiCGSTAB for CSR SLAU matrix parallized with MKL
bool BiCGSTAB_CSR_MKL(int N, const double *A, const double *b, const double *newx0, const double *RealSol,
	const int *iptr, const int *jptr, int &niter, double &eps, int &stopres,
	vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit, const int threadNum, double* result);

// Origin BiCGSTAB for CSR SLAU matrix with Yacobi preconditioner
vector<double> BiCGSTAB_Yacobi_CSR(const vector<double> &A, const vector<double> &b, const vector<double> &newx0, const vector<double> &RealSol,
	const vector<int> &iptr, const vector<int> &jptr, int &niter, double &eps, int &stopres,
	vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit);


// ************************************************************************************************************************************************************************
// *********************************************************************** Other Solvers **********************************************************************************
// ************************************************************************************************************************************************************************

//метод  Простой Итерации
vector<vector<double>> SimpleIter(const vector<vector<double>> &A, const vector<double> &f, const  vector<double> &x0, double eps);

//41метод гаусса прямой и обратный ход
vector<double> Gauss_straight(vector<vector<double>> &A, vector<double> &b);
//42
void GaussSolver(vector<vector<double>> &SLAU, vector<double> &right, vector<double> &futureresult);


// ************************************************************************************************************************************************************************
// ********************************************************************* Preconditioners **********************************************************************************
// ************************************************************************************************************************************************************************
void diag_normalize(const int &rank, vector<double> &SLAU,
	const vector<int> &iptr, const  vector<int> &jptr, vector<double> &right);

void Yacobi_precondition(const int &rank, const vector<double> &SLAU,
	const vector<int> &iptr, const  vector<int> &jptr, vector<double> &in, vector<double> &out);

void reverseMatrixByIterShultze(const vector<vector<double>> &A, vector<vector<double>> &revA);