// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.

#ifndef PCH_H
#define PCH_H

#define PI 3.1415926535
#include "locale.h"
#include "math.h"
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <locale>
#include <omp.h>
#include <set>
#include <stdio.h>
#include <cstdlib>
#include<stdlib.h>
#include<string.h>
#include <tchar.h>
#include <time.h>
#include <vector>


#include "profile.h"

// #include <mkl.h>
//#include<mkl_cblas.h>
//#include <mkl_blas.h>
//#include <mkl_spblas.h>
#include <corecrt_wstring.h>
#include <direct.h>
#include<filesystem>


typedef double(*scalar_function) (const double &, const double&);
typedef std::vector<double>(*vector_function) (const double &, const double&);
typedef std::vector<double>(*vect_function) (const double &, const double&);
typedef double(*parameters_function) (const double &, const double&, const double&);

typedef double(*spec_scalar_function) (const double &, const double &,
	const vector<double> &, const vector_function &);



#include"VectorMathOperations.h"

#include "targetver.h"
#include "coutObjects.h"
#include "input_output.h"


#include"structures.h"
#include"Mesh.h"
#include"Solvers.h"

#include"CSR.h"

#include"mesh_trg.h"


#include "basisFunction.h"
#include"basisLagrangeQuad.h"
#include "basisLagrangeTrg.h"


#include "AbstractModel.h"

#include "SolutionDiscreteAnalogFactory.h"
#include "SolutionDiscreteAnalogFactoryCG.h"
#include "SolutionDiscreteAnalogFactoryDG.h"

#include "SolutionDiscreteAnalog.h"
#include "SolutionDiscreteAnalogCG.h"
#include "AdvectiveDiffusionSDA.h"
#include"AdvectiveDiffusionSDAplusSUPG.h"


#include "SolutionDiscreteAnalogDG.h"
#include"AdvecDiffusDG_SDA.h"

#include "SolutionCreater.h"
#include "SolutionCreaterCG.h"
#include "SolutionCreaterDG.h"
#include"TimeDepSolutionCreaterCG.h"
#include"MP1biquadAS.h"
#include"MP2withBndryLayers.h"
#include"MP3SinglrtyOnBndry.h"
#include"MP4withConstBC1.h"
#include"MP5withConstBC1on2bndrs.h"
#include"MP6DGwithSinglrtyAndDiscAS.h"
#include"MP7withBndryLayer.h"
#include"MP8constSource.h"
#include"RP1_filterInPorousMedia.h"

#include"launchOOP.h"

using std::vector;
using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::set;

int findElement(const vector<int> &a, int st, int end, int value);

vector<double> BiCGSTAB_CSR(const vector<double> &A, const vector<double> &b,
	const vector<double> &newx0, const vector<double> &RealSol,
	const vector<int> &iptr, const vector<int> &jptr, int &niter, double &eps, int &stopres,
	vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit);
vector<double> BiCGSTAB_CSR_OMP(const vector<double> &A, const vector<double> &b, const vector<double> &newx0, const vector<double> &RealSol,
	const vector<int> &iptr, const vector<int> &jptr, int &niter, double &eps, int &stopres,
	vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit, const int threadNum);

vector<double> BiCGSTAB_Yacobi_CSR(const vector<double> &A, const vector<double> &b, const vector<double> &newx0, const vector<double> &RealSol,
	const vector<int> &iptr, const vector<int> &jptr, int &niter, double &eps, int &stopres, vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit);

void Yacobi_precondition(const int &rank, const vector<double> &SLAU,
	const vector<int> &iptr, const  vector<int> &jptr, vector<double> &in, vector<double> &out);

void diag_normalize(const int &rank, vector<double> &SLAU, const vector<int> &iptr, const  vector<int> &jptr, vector<double> &right);


// TODO: add headers that you want to pre-compile here

#endif //PCH_H
