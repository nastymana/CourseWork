#pragma once

class basisFunction
{
public:
	typedef double(*ScalarFunction) (const double &, const double&);
	basisFunction();
	basisFunction(const int &rankItem);
	virtual ~basisFunction();
protected:	
	int GaussRank;
	// ранг базиса: 1 - линейный(билинейный), 2 - квадратичный (биквадратичный), 3 - кубический(бикубический
	int basisRank; 
	// кол-во БФ  в базисе (билинейны - 4; биквадратичный - 9; бикубический - 16)
	int basisSize; 
	vector<scalar_function> fi; // БФ скалярного базиса
	vector<vector_function> psi;// для векторных базисов
	vector<vector_function> grad; // градиент для функций скалярного базиса.
	
	// Лапласиан функции. равен (d2Fi/dx2 + d2Fi/dy2). для удобаства, как и grad  разделены:
	// laplacia[0] = d2Fi/dx2; laplacia[1] = d2Fi/dy2; (для прямоугольных элементов
	vector<vector_function> laplacian; 
	vector<double> carrierX, carrierY;
	vector<double> GaussPoints;
	vector<double> GaussWeights;
public:
	// возвращает значения every BF в точках гаусса на МЭ. 
	// размерность возвращ вектора равна числу БФ, размерность внутр строк - числу точек Гаусса на МЭ
	vector<vector<double>> getFiValuesInGaussPoints();
	vector<vector<double>> getFiValuesInGaussPoints(int nBndry);
	// 
	vector<double> getFiValuesInPoint(const double &eps, const double &eta);
	// возвращает значения every gradBF в точках гаусса на МЭ. ;  // Размерность: NBF*dim*NGP (dim=2 (2D))
	vector<vector<vector<double>>> getGradFiValuesInGaussPoints();
	vector<vector<vector<double>>> getGradFiValuesInGaussPoints(int nBndry);
	// возвращает значения every LaplacianBF в точках гаусса на МЭ.  // Размерность: NBF*dim*NGP (dim=2 (2D))
	vector<vector<vector<double>>> getLaplacFiValuesInGaussPoints(); 
	// для 2D эл-ов возвращает перемножение весов (те размер вектора дб равен числу точек гаусса для 2D  эл-та
	vector<double> getGaussWeights(); 
	vector<double> getGaussWeights(const int& nBndry);
	vector<double> getCarrier(const char &XorYitem);
	virtual double epsToX(const double &x1, const double &x2, const double &eps)=0;
	virtual double XtoEps(const double & x1, const double & x2, const double & x)=0;
	/* мастер элемент:
	х1,у2 ...... х2,у2
	.				.
	.				.
	.				.
	х1,у1 ...... х2,у1

	___________________________

	х1,у3 . . х2,у3	. . х3,у3
	.						.
	.						.
	х1,у2 . . х2,у2	. . х3,у2
	.						.
	.						.
	х1,у1 . . х2,у1	. . х3,у1
	*/
	// возвращает пары точек Гаусса для МЭ в порядке  слева направо снизу вверх {{х1,у1}, {x2,y1}, ...{x1,y2},{x2,y2}...{xn, yn}}
	vector<vector<double>> getPairsOfGaussPoints();
	vector<double> getGaussPoints(const char &XorYitem);

	int getbasisSize();
	int getbasisRank();
	int getGaussPointNumber();
	void getXYforGaussPoints(const quadrElem &elem, vector<vector<double>>& XYpairs);
	void getXYforGaussPoints(const int &nBnd, const quadrElem & elem, vector<vector<double>>& coord);
	void getXorYforCarrier(const char &XorYitem, const quadrElem &elem, vector<double>& coord);
	void getXorYforGaussPoints(const char &XorYitem, const quadrElem &elem, vector<double>& coord);
	virtual double getDetJacobian(const vector<double> &x, const vector<double> &y);

	virtual vector<vector<double>>  getJacobian(const vector<double> &x, const vector<double> &y);
	virtual void coutBasisParams()=0;

};


