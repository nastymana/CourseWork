//

#include "pch.h"
#include"VectorMathOperations.h"

//26
void innerProduct(const vector<vector<double>> &A, const vector<vector<double>> &B,
	vector<vector<double>> &result) {
	int Ncolumn = B.size(), Nstring = A.size();
	result = vector<vector<double>>(Nstring, vector<double>(Ncolumn));
	vector<vector<double>> Bt = transposeMatrix(B); // transposed B
	for (size_t i = 0; i < Nstring; i++)
	{
		for (size_t j = 0; j < Ncolumn; j++)
		{
			result[i][j] = innerProduct(A[i], Bt[j]);
		}
	}
}

double innerProduct(const vector<double> &a, const vector<double> &b)
{
	double result = 0.;
	int n = a.size();
	if (a.size() == b.size())
		//cout << "in prod  " << n<<endl;
	{
		for (int i = 0; i < int(a.size()); i++)
		{
			result += (a[i] * b[i]);
		}
		return result;
	}
	else
	{
		cout << "size of matrices does not equal\n";
		//break;
	}

}

double inner_product_OMP(const vector<double> &a, const vector<double> &b, const int threadNum)
{
	double result = 0.;
	int n = a.size(), i=0;
	if (a.size() == b.size())
		//cout << "in prod  " << n<<endl;
	{
		
#pragma omp parallel for shared (a, b) private(i) reduction(+:result) num_threads(threadNum)
		for (i = 0; i < int(a.size()); i++)
		{
			result += (a[i] * b[i]);
		}
		return result;
	}
	else
	{
		cout << "size of matrices does not equal\n";
		//break;
	}

}

vector<double> innerProduct (const vector<vector<double>> &A, const vector<double> &B)
{
	vector<double> result;
	if (A[0].size() == B.size())
	{
		for (int i = 0; i < int(A.size()); i++)
			{
				result.push_back(innerProduct(A[i], B));
			}
	return result;
	}
	else
	{
		cout << "size of matrices does not equal\n";
		//break;
	}
}

//void inner_product(const vector<vector<double>>& A, const vector<vector<double>>& B, vector<vector<double>>& result)
//{
//	if (result.size() != A.size()) {
//		result = vector<vector<double>>(B.size(), vector<double> (A.size(), 0.));
//}
//	if(A[0].size()==B.size()){
//	vector<vector<double>> Btransp = transposeMatrix(B);
//	for (int i = 0; i < int(A.size()); i++)
//	{
//		for (int j = 0; j < int(B[].size()); j++) {
//			result[i][j] = (inner_product(A[i], Btransp[i]));
//		}
//	}
//	}
//	else
//	{
//		cout << "size of matrices does not equal\n";
//		//break;
//	}
//}

double summVectorElms(const vector<vector<double>> &A)
{
	double result = 0;
	for (int i = 0; i<int(A.size()); i++)
	{
		for (int j = 0; j <int(A[i].size()); j++)
		{
			result += A[i][j];
		}
	}
	return result;
}
double summVectorElms(const vector<double> &A)
{
	double result = 0;
	for (int i = 0; i<int(A.size()); i++)
	{
		result += A.at(i);
	}
	return result;
}
//8
double vectorNorma(const vector<double> &a)
{
	double norma = innerProduct(a,a);
	norma = sqrt(abs(norma));
	return norma;
}
double vector_norma_OMP(const vector<double> &a, const int threadNum)
{
	double norma = inner_product_OMP(a, a, threadNum);
	norma = sqrt(abs(norma));
	return norma;
}

double average(const vector<double> &a)
{
	double aver = 0.;
	for (int i = 0; i < int(a.size()); i++)
	{
		aver += a[i];
	}
	return aver/double(a.size());
}
//10
double angle_between_vectors(const vector<double> &a, const vector<double> &b)
{
	double angle;
	double norma_a = sqrt(innerProduct(a, a));
	double norma_b = sqrt(innerProduct(b, b));
	double ab = innerProduct(a, b);
	//cout << "inner_product = " << ab << endl;
	angle = ab / (norma_a*norma_b);
	angle = acos(angle) * 180 / PI;
	return angle;
}
vector<double> vector_normalyze(const vector<double> &a)
{
	int n = a.size();
	double norma = sqrt(vectorNorma(a));
	vector<double> result;
	for (int i = 0; i < n; i++)
	{
		result[i] = a[i] / norma;
	}
	return result;
}
vector<vector<double>> operator+ (const vector<vector<double>> &A, const vector<vector<double>> &B)
{
	vector<vector<double>> result(int(A.size()));// = A;
	for (int i = 0; i < int(A.size()); i++)
	{
		//cout << "A[i] " << A[i].size() << "B[i] = " << B[i].size() << endl;
		for (int j = 0; j < int(A[i].size()); j++)
		{
			result[i].push_back(A[i][j] + B[i][j]);
		}
	}
	return result;
}
vector<vector<int>> operator+ (const vector<vector<int>> &A, const vector<vector<int>> &B)
{
	vector<vector<int>> result;// = A;
	for (int i = 0; i < int(A.size()); i++)
	{
		//cout << "A[i] " << A[i].size() << "B[i] = " << B[i].size() << endl;
		for (int j = 0; j < int(A[i].size()); j++)
		{
			result[i][j] = A[i][j] + B[i][j];
		}
	}
	return result;
}
vector<double> operator+ (const vector<double> &A, const vector<double> &B)
{
	vector<double> result;

	//cout << "A[i] " << A[i].size() << "B[i] = " << B[i].size() << endl;
	for (int j = 0; j < int(A.size()); j++)
	{
		result.push_back(A[j] + B[j]);
	}

	return result;
}

void ax_plus_by_OMP (int a, const vector<double> &A, int b, const vector<double> &B, vector<double> &result, const int threadNum)
{
	int rank = int (A.size()), j=0;
	//cout << "A[i] " << A[i].size() << "B[i] = " << B[i].size() << endl;
#pragma omp parallel for shared (a, A, b, B) private(j) num_threads(threadNum)
	for (j = 0; j < rank; j++)
	{
		result[j] = a*A[j] + b*B[j];
	}
	
}
void ax_plus_by_plus_cz_OMP(int a, const vector<double> &A, int b, const vector<double> &B, 
	int c, const vector<double> &C, vector<double> &result, const int threadNum)
{
	int j = 0;
	//cout << "A[i] " << A[i].size() << "B[i] = " << B[i].size() << endl;
#pragma omp parallel for shared (a, A, b, B, c, C) private(j) num_threads(threadNum)
	for (j = 0; j < int(A.size()); j++)
	{
		result[j]=a*A[j] + b*B[j]+c*C[j];
	}

}


vector<vector<double>> operator- (const vector<vector<double>> &A, const vector<vector<double>> &B)
{
	vector<vector<double>> result(int(A.size()));// = A;
	for (int i = 0; i < int(A.size()); i++)
	{
		//cout << "A[i] " << A[i].size() << "B[i] = " << B[i].size() << endl;
		for (int j = 0; j < int(A[i].size()); j++)
		{
			result[i].push_back(A[i][j] - B[i][j]);
		}
	}
	return result;
}
vector<vector<int>> operator- (const vector<vector<int>> &A, const vector<vector<int>> &B)
{
	vector<vector<int>> result;// = A;
	for (int i = 0; i < int(A.size()); i++)
	{
		//cout << "A[i] " << A[i].size() << "B[i] = " << B[i].size() << endl;
		for (int j = 0; j < int(A[i].size()); j++)
		{
			result[i][j] = A[i][j] - B[i][j];
		}
	}
	return result;
}
vector<double> operator- (const vector<double> &A, const vector<double> &B)
{
	vector<double> result;

	//cout << "A[i] " << A[i].size() << "B[i] = " << B[i].size() << endl;
	for (int j = 0; j < int(A.size()); j++)
	{
		result.push_back(A[j] - B[j]);
	}

	return result;
}

vector<double> operator* (const vector<double> &A, double c)
{
	vector<double> result;
	for (int i = 0; i < int(A.size()); i++)
	{
		result.push_back( A[i] * c);
	}
	return result;
}
vector<double> operator* ( double c, const vector<double> &A)
{
	vector<double> result;
	for (int i = 0; i < int(A.size()); i++)
	{
		result.push_back(A[i] * c);
	}
	return result;
}
vector<vector<double>> operator* (const vector<vector<double>> &A, double c)
{
	vector<vector<double>> result(int(A.size()));
	for (int i = 0; i < int(A.size()); i++)
	{
		for (int j = 0; j < int(A[i].size()); j++)
		{
			result[i].push_back( A[i][j] * c);
		}
	}
	return result;
}
vector<vector<double>> operator* (double c,const vector<vector<double>> &A)
{
	vector<vector<double>> result(int(A.size()));
	for (int i = 0; i < int(A.size()); i++)
	{
		for (int j = 0; j < int(A[i].size()); j++)
		{
			result[i].push_back(A[i][j] * c);
		}
	}
	return result;
}
vector<vector<double>> operator* (const vector<vector<double>> &A, const vector<vector<double>> &B)
{
	vector<vector<double>> result(A.size(), vector<double> (B.size()));
	for (int i = 0; i < int(A.size()); i++)
	{
		for (int j = 0; j < int(A[i].size()); j++)
		{
			result[i][j] = (A[i][j] * B[i][j]);
		}
	}
	return result;
}

vector<double> operator* (const vector<double> &A, const vector<double> &B)
{
	vector<double> result(A.size());
	for (int i = 0; i < int(A.size()); i++)
	{
		result[i] = (A[i] * B[i]);
	}
	return result;
}
vector<vector<double>> operator* (const vector<vector<double>> &A, const vector<double> &c)
{
	vector<vector<double>> result;
	for (int i = 0; i < int(A.size()); i++)
	{
		for (int j = 0; j < int(A[i].size()); j++)
		{
			result.push_back(A[i] *c );
		}
	}
	return result;
}

void sewvectors(vector<double> &a, const  vector<double> &b)
{
	int a_size = a.size();
	a.resize(a_size + b.size());
//	coutVector(a, "a resized", 'H');
	for (int i = a_size, j=0; j < int(b.size()); i++, j++){
	//	cout << "i=" << i << ", j=" << j << endl;
		a[i] = b[j];
	}	
}
vector<double> pushvector(const vector<double> &a, double push)
{
	vector<double> result=a;
	for (int i = 0; i < a.size(); i++)
	{
		result[i] += push;
	}
	return result;
}
void sewvectors(vector<int> &a, const vector<int> &b)
{
	int a_size = a.size();
	a.resize(a_size + b.size());
	for (int i = a_size, j = 0; j < int(b.size()); i++, j++) {
		a[i] = b[j];
	}
}
vector<double> multiply_vectorsn_to_nn(const vector<double> &a, const vector<double> &b)
{
	vector<double> result;
	for (int i = 0; i < int(a.size()); i++)
	{
		//cout << "a->size" << a.size() << endl;
		for (int j = 0; j < int(b.size()); j++)
		{
			result.push_back(a.at(i) * b.at(j));
		}
	}
	return result;
}
//35
vector<vector<double>> Uptomirror(const vector<vector<double>> &A)
{//ИСПОЛЬЗОВАТЬ ФУНКЦИИ КЛАССА VECTOR, а не костыли. 
 //либо создать пп для  заполнения 0 МНОГОМЕРНЫХ МАССИВОВ
	int N = int(A.size());
	vector<double> help(N, 0);
	vector<vector<double>> result(N, help);
	for (int i = 0; i <N; i++)
	{
		result[i] = help;
	}
	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j <int(A[i].size()); j++)
		{
			result[i][i + j] = A[i][j];
		}
	}
	for (int i = N - 1; i >= 0; i--)
	{
		for (int j = 0; j < i; j++)
		{
			result[i][j] = result[j][i];
		}
	}
	return result;
}

vector<vector<double>> transposeMatrix(const vector<vector<double>>A)
{//N*M -> M*N
	vector<vector<double>> result;
	for (int i = 0; i < int(A[0].size()); i++)
	{
		vector<double> help;
		for (int j = 0; j <int(A.size()); j++)
		{
			help.push_back(A[j][i]);
		}
		result.push_back(help);
	}
	return result;
}
//48
double lenght(const double &x1, const  double &y1, const double &x2, const double &y2)
{
	double norma = pow((x1 - x2), 2) + pow((y1 - y2), 2);
	norma = sqrt(norma);
	return norma;
}
vector<vector<double>>sewvectors( vector<vector<double>> &a, const vector<vector<double>>&b)
{
	vector<vector<double>> result=a;
	for (int i = 0; i < int(b.size()); i++)
	{
		result.push_back(b[i]);

	}
	return result;
}

int min_abs(const std::vector<int> inbox, int start, int end) {
	int Number=0;
	return Number;
}
int max_abs(const std::vector<int> inbox, int start, int end) {
	int Number=0;
	return Number;
}
// меняем местами столбцы i, k
void permute_column(std::vector<std::vector<int>> &inbox, int i, int k) {
	int M = inbox.size(), N = inbox[0].size();

	int copy_i = 0, copy_k = 0;
	for (int j = 0; j < M; j++)
	{
		copy_i = inbox[j][i];
		inbox[j][i] = inbox[j][k];
		inbox[j][k] = copy_i;
	}
}

// меняем местами строки i, k
void permute_string(std::vector<std::vector<int>> &inbox, int i, int k) {}
// вычитаем из столбца i k-ый*factor
void subtract_column(std::vector<std::vector<int>> &inbox, int i, int k, int factor) {
	int M = inbox.size(), N = inbox[0].size();

	for (int j = 0; j < M; j++)
	{
		inbox[j][i] -= inbox[j][k] * factor;
	}
}




// На входе получаем матрицу размерности (Nex+Nx) x (Nx+1)
void to_doun_triangle(std::vector<std::vector<int>> &inbox, const int Nx, const int Nequat) {
	// Сначала приведем к трапеция-виду матрицу А.
	int  Nstring_min, Nstring_max,
		dev_int = 1, rest_int = 0, start, end;

	for (int i = 0; i < Nequat - 1; i++)
	{
		int Null_count = 0;
		for (int j = 0; j < Nx; j++)
		{// проверка на нули
			if (inbox[i][j] == NULL) {
				permute_column(inbox, i, (Nx - 1 - Null_count)); //Проверить ограничение
				Null_count++;
			}
		}
		while (Null_count != (Nx - 1 - i))// Проверить ограничение
		{
			start = i + 1;
			end = Nx - 1 - i;  /// NB!!!  - просто случайно сюда поставила. проверить!!!
			Nstring_min = min_abs(inbox[i], start, end);// 
			Nstring_max = max_abs(inbox[i], start, end);
			dev_int = inbox[i][Nstring_max] / inbox[i][Nstring_min];
			rest_int = inbox[i][Nstring_max] % inbox[i][Nstring_min]; // Для проверки что появился новый ноль в строке
			subtract_column(inbox, Nstring_max, Nstring_min, dev_int);
			if (rest_int == NULL) {
				permute_column(inbox, Nstring_max, Nx - Null_count - 1);
				Null_count++;
			}
		}
	}

	// ПОСЛЕ - занулим правую часть
	int Null_count = 0, i = 0;
	{
		Nstring_min = min_abs(inbox[i], start, end);// 
		dev_int = inbox[i][Nx] / inbox[i][Nstring_min];
		rest_int = inbox[i][Nx] % inbox[i][Nstring_min]; // Для проверки что появился новый ноль в строке
		while (Null_count != NULL) {
			subtract_column(inbox, Nstring_max, Nx - Null_count - 1, 1);
			Null_count++;
		}
		i++;
	}


}

void reshape_vector(const vector<vector<double>> &a, vector<double> &result) {
	result = a[0];
	for (size_t i = 1; i < a.size(); i++)
	{
		sewvectors(result, a[i]);
	}
}

void reshape_vector(const vector<vector<int>> &a, vector<int> &result) {
	result = a[0];
	for (size_t i = 1; i < a.size(); i++)
	{
		sewvectors(result, a[i]);
	}
}

void degreeVector(const vector<double>& a, vector<double>& result, const int &deg)
{
	result = vector<double>(a.size());
	for (size_t i = 0; i < a.size(); i++){
		result[i] = pow(a[i], deg);
	}
}
