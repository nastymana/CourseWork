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
	// ���� ������: 1 - ��������(����������), 2 - ������������ (��������������), 3 - ����������(������������
	int basisRank; 
	// ���-�� ��  � ������ (��������� - 4; �������������� - 9; ������������ - 16)
	int basisSize; 
	vector<scalar_function> fi; // �� ���������� ������
	vector<vector_function> psi;// ��� ��������� �������
	vector<vector_function> grad; // �������� ��� ������� ���������� ������.
	
	// ��������� �������. ����� (d2Fi/dx2 + d2Fi/dy2). ��� ���������, ��� � grad  ���������:
	// laplacia[0] = d2Fi/dx2; laplacia[1] = d2Fi/dy2; (��� ������������� ���������
	vector<vector_function> laplacian; 
	vector<double> carrierX, carrierY;
	vector<double> GaussPoints;
	vector<double> GaussWeights;
public:
	// ���������� �������� every BF � ������ ������ �� ��. 
	// ����������� ������� ������� ����� ����� ��, ����������� ����� ����� - ����� ����� ������ �� ��
	vector<vector<double>> getFiValuesInGaussPoints();
	vector<vector<double>> getFiValuesInGaussPoints(int nBndry);
	// 
	vector<double> getFiValuesInPoint(const double &eps, const double &eta);
	// ���������� �������� every gradBF � ������ ������ �� ��. ;  // �����������: NBF*dim*NGP (dim=2 (2D))
	vector<vector<vector<double>>> getGradFiValuesInGaussPoints();
	vector<vector<vector<double>>> getGradFiValuesInGaussPoints(int nBndry);
	// ���������� �������� every LaplacianBF � ������ ������ �� ��.  // �����������: NBF*dim*NGP (dim=2 (2D))
	vector<vector<vector<double>>> getLaplacFiValuesInGaussPoints(); 
	// ��� 2D ��-�� ���������� ������������ ����� (�� ������ ������� �� ����� ����� ����� ������ ��� 2D  ��-��
	vector<double> getGaussWeights(); 
	vector<double> getGaussWeights(const int& nBndry);
	vector<double> getCarrier(const char &XorYitem);
	virtual double epsToX(const double &x1, const double &x2, const double &eps)=0;
	virtual double XtoEps(const double & x1, const double & x2, const double & x)=0;
	/* ������ �������:
	�1,�2 ...... �2,�2
	.				.
	.				.
	.				.
	�1,�1 ...... �2,�1

	___________________________

	�1,�3 . . �2,�3	. . �3,�3
	.						.
	.						.
	�1,�2 . . �2,�2	. . �3,�2
	.						.
	.						.
	�1,�1 . . �2,�1	. . �3,�1
	*/
	// ���������� ���� ����� ������ ��� �� � �������  ����� ������� ����� ����� {{�1,�1}, {x2,y1}, ...{x1,y2},{x2,y2}...{xn, yn}}
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


