#include"pch.h"
#include"Solvers.h"



// ************************************************************************************************************************************************************************
// ************************************************************* BiCgSTAB solvers family **********************************************************************************
// ************************************************************************************************************************************************************************

// 6/02/19  - сделать, но не сегодня интерфейс для решателя, где можно было бы подавать параметр решателя. (dense. csr, preconditioner, etc.) - интерфейс. 
// Либо отложить это до этапа реализации комплекса в ООП парадигме
// Origin BiCGSTAB for dense SLAU matrix
vector<double> BiCGSTAB_dense(const vector<vector<double>> &A, const vector<double> &b, const vector<double> &newx0, const vector<double> &RealSol,
	int &niter, double &eps, int &stopres, vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit) {
	/*"""Программа принимает на вход систему Ax=b,
	причем А = квадратная матрица NxN элементов, хранимая в CSR формате.
	Программа возвращает вектор решений х;
	номер итерации на которой провесс расчет остановился;
	вектор содержащий норму невязок :*/
	std::cout << "BiCGSTAB Solution start" << endl;
	vector <double> x0 = newx0;
	int N = b.size();
	std::cout << "N" << N << endl;
	//	coutVector(A, "A", 'H');
	//	coutVector(b, "b", 'H');
	vector<double>	r0 = (b - innerProduct(A, x0));
	//	coutVector(CSRmultiMV(A, x0, iptr, jptr), "CSRmultiMV(A, x0, iptr, jptr)", 'H');
	//coutVector(r0, "r0", 'H');
	vector<double>	_r(N, 1.),
		p0 = r0,
		x1,
		Ap,
		s,
		As,
		res,
		r1;
	double	normab = vectorNorma(b),
		r0_r = 0.,
		Apr = 0.,
		alfa = 0.,
		sqrAs = 0.,
		w = 0.,
		betta = 0.;

	for (int i = 0; i < niter; i++) {
		cout << "i::" << i;
		Ap = innerProduct(A, p0);
		//	coutVector(Ap, "Ap", 'H');
		//	cout << Ap.size() << "  " << endl;
		r0_r = innerProduct(r0, _r);
		//	cout << "r0_r "<< r0_r << endl;
		Apr = innerProduct(Ap, _r);
		//	cout << "Apr " << Apr << endl;
		if (Apr == 0)
		{
			std::cout << "devision on (Ap,r) = 0)" << endl;
			std::cout << "сошлось при nitr=" << i << ",  ||res|| = " <<
				//			normaResit[i-1] << " ,  ||R1|| = " 
				//			<< normaR1it[i-1] <<
				endl;
			return x0;// , i - 1, normaR1it, normaResit, exResit)
		}
		alfa = r0_r / Apr;
		s = r0 - alfa * Ap;
		As = innerProduct(A, s);
		sqrAs = innerProduct(As, As);
		if (sqrAs <= 10e-20) {
			//std::cout << "devision on (As,As) = 0)" << endl;
			//std::cout << "сошлось при ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x0;// , i - 1, normaR1it, normaResit, exResit)
		}
		w = innerProduct(As, s) / sqrAs;
		x1 = x0 + alfa * p0 + w * s;
		res = b - innerProduct(A, x1);
		if (i >= stopres) { r1 = res; }
		else { r1 = s - w * As; }
		normaR1it.push_back(vectorNorma(r1));
		normaResit.push_back(vectorNorma(res));
		cout << ", res = " << normaResit[i] << endl;
		//	std:cout << "RealSol.size() " << RealSol.size() << "  x1.size() " << x1.size() << endl;
	//	exResit.push_back(vector_norma(RealSol - x1));
		if (normaR1it[i] / normab <= eps) {
			std::cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x1;// , i, normaR1it, normaResit, exResit)
		}
		if (r0_r == 0) {
			std::cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x1;// , i, normaR1it, normaResit, exResit)
		}
		if (w == 0) {
			std::cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x1;// , i, normaR1it, normaResit, exResit)
		}
		betta = (innerProduct(r1, _r) / r0_r)*(alfa / w);
		p0 = r1 + betta * (p0 - w * Ap);
		r0 = r1;
		x0 = x1;
	//	coutVector(x1, "sol", 'H');

	}
	int n = normaResit.size() - 1;
	std::cout << "сошлось при nitr=" << n << ",  ||res|| = " << normaResit[n] << " ,  ||R1|| = " << normaR1it[n] << endl;
	return x1;// , niter, normaR1it, normaResit, exResit)
}

// Origin BiCGSTAB for CSR SLAU matrix
vector<double> BiCGSTAB_CSR(const vector<double> &A, const vector<double> &b, const vector<double> &newx0, const vector<double> &RealSol,
	const vector<int> &iptr, const vector<int> &jptr, int &niter, double &eps, int &stopres,
	vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit) {
	/*"""Программа принимает на вход систему Ax=b,
	причем А = квадратная матрица NxN элементов, хранимая в CSR формате.
	Программа возвращает вектор решений х;
	номер итерации на которой провесс расчет остановился;
	вектор содержащий норму невязок :*/
	cout << "BiCGSTAB Solution start" << endl;
	vector <double> x0 = newx0;
	int N = b.size();
	cout << "N" << N << endl;
	//	coutVector(A, "A", 'H');
	//	coutVector(b, "b", 'H');
	vector<double>	r0 = (b - CSRmultiMV(A, x0, iptr, jptr));
	//	coutVector(CSRmultiMV(A, x0, iptr, jptr), "CSRmultiMV(A, x0, iptr, jptr)", 'H');
	//	coutVector(r0, "r0", 'H');
	vector<double>	_r(N, 1.),
		p0 = r0,
		x1,
		Ap,
		s,
		As,
		res,
		r1;
	double	normab = vectorNorma(b),
		r0_r = 0.,
		Apr = 0.,
		alfa = 0.,
		sqrAs = 0.,
		w = 0.,
		betta = 0.;

	for (int i = 0; i < niter; i++) {
		Ap = CSRmultiMV(A, p0, iptr, jptr);
		//	coutVector(Ap, "Ap", 'H');
		//	cout << Ap.size() << "  " << endl;
		r0_r = innerProduct(r0, _r);
		//	cout << "r0_r "<< r0_r << endl;
		Apr = innerProduct(Ap, _r);
		//	cout << "Apr " << Apr << endl;
		if (Apr == 0)
		{
			cout << "devision on (Ap,r) = 0)" << endl;
			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x0;// , i - 1, normaR1it, normaResit, exResit)
		}
		alfa = r0_r / Apr;
		s = r0 - alfa * Ap;
		As = CSRmultiMV(A, s, iptr, jptr);
		sqrAs = innerProduct(As, As);
		if (sqrAs <= 10e-20) {
			//cout << "devision on (As,As) = 0)" << endl;
			//cout << "сошлось при ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x0;// , i - 1, normaR1it, normaResit, exResit)
		}
		w = innerProduct(As, s) / sqrAs;
		x1 = x0 + alfa * p0 + w * s;
		res = b - CSRmultiMV(A, x1, iptr, jptr);
		if (i >= stopres) { r1 = res; }
		else { r1 = s - w * As; }
		normaR1it.push_back(vectorNorma(r1));
		normaResit.push_back(vectorNorma(res));
		//	std:cout << "RealSol.size() " << RealSol.size() << "  x1.size() " << x1.size() << endl;
		exResit.push_back(vectorNorma(RealSol - x1));
		if (normaR1it[i] / normab <= eps) {
			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x1;// , i, normaR1it, normaResit, exResit)
		}
		if (r0_r == 0) {
			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x1;// , i, normaR1it, normaResit, exResit)
		}
		if (w == 0) {
			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x1;// , i, normaR1it, normaResit, exResit)
		}
		betta = (innerProduct(r1, _r) / r0_r)*(alfa / w);
		p0 = r1 + betta * (p0 - w * Ap);
		r0 = r1;

		x0 = x1;

	}
	int n = normaResit.size() - 1;
	cout << "сошлось при nitr=" << n << ",  ||res|| = " << normaResit[n] << " ,  ||R1|| = " << normaR1it[n] << endl;
	return x1;// , niter, normaR1it, normaResit, exResit)
}


void BiCGSTAB_CSR(SLAE_CSR* slae,	const vector<double> &newx0, const vector<double> &RealSol,
	int &niter, double &eps, int &pivot,
	vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit) {
	/*"""Программа принимает на вход систему Ax=b,
	причем А = квадратная матрица NxN элементов, хранимая в CSR формате.
	Программа возвращает вектор решений х;
	номер итерации на которой провесс расчет остановился;
	вектор содержащий норму невязок :*/
	cout << "BiCGSTAB Solution start" << endl;
	vector <double> x0 = newx0;
	int N = slae->b.size();
	cout << "N" << N << endl;
	//	coutVector(A, "A", 'H');
	//	coutVector(b, "b", 'H');
	vector<double>	r0 = (slae->b - CSRmultiMV(slae, x0 ));
	//	coutVector(CSRmultiMV(A, x0, iptr, jptr), "CSRmultiMV(A, x0, iptr, jptr)", 'H');
	//	coutVector(r0, "r0", 'H');
	vector<double>	_r(N, 1.),
		p0 = r0,
		x1,
		Ap,
		s,
		As,
		res,
		r1;
	double	normab = vectorNorma(slae->b),
		r0_r = 0.,
		Apr = 0.,
		alfa = 0.,
		sqrAs = 0.,
		w = 0.,
		betta = 0.;

	for (int i = 0; i < niter; i++) {
		
		Ap = CSRmultiMV(slae, p0);
		//	coutVector(Ap, "Ap", 'H');
		//	cout << Ap.size() << "  " << endl;
		r0_r = innerProduct(r0, _r);
		//	cout << "r0_r "<< r0_r << endl;
		Apr = innerProduct(Ap, _r);
		cout << "Apr " << Apr << endl;
		if (abs(Apr) <= 10e-20){
			cout << "devision on (Ap,r) = 0)" << endl;
		//	cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i-1] << " ,  ||R1|| = " << normaR1it[i-1] << endl;
			slae->setX(x0);// , i - 1, normaR1it, normaResit, exResit)
			i = niter;
		}		
		else {
			alfa = r0_r / Apr;
			s = r0 - alfa * Ap;
			As = CSRmultiMV(slae, s);
			sqrAs = innerProduct(As, As);
			if (abs(sqrAs) <= 10e-20) {
				cout << "devision on (As,As) = 0)" << endl;
				//cout << "сошлось при ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
				slae->setX(x0);// , i - 1, normaR1it, normaResit, exResit)
				i = niter;
			}
			else {
				w = innerProduct(As, s) / sqrAs;

				// coutVector(alfa * p0 + w * s, "alfa * p0 + w * s", 'H', 5);
				x1 = x0 + alfa * p0 + w * s;
				res = slae->b - CSRmultiMV(slae, x1);
				if (i == pivot) { r1 = res; }
				else { r1 = s - w * As; }
				normaR1it.push_back(vectorNorma(r1));
				normaResit.push_back(vectorNorma(res));
				//	std:cout << "RealSol.size() " << RealSol.size() << "  x1.size() " << x1.size() << endl;
				exResit.push_back(vectorNorma(RealSol - x1));
				if (normaR1it[i] / normab <= eps) {
					cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
					slae->setX(x1);// , i, normaR1it, normaResit, exResit)
					i = niter;
				}
				else {
					if ((abs(r0_r) <= 10e-20) || (abs(w )<= 10e-20)) {
						cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
						slae->setX(x1);// , i, normaR1it, normaResit, exResit)
						i = niter;

					}
					else {
						betta = (innerProduct(r1, _r) / r0_r)*(alfa / w);
						p0 = r1 + betta * (p0 - w * Ap);
						r0 = r1;
						x0 = x1;
						int n = normaResit.size() - 2;
						cout << "iter:: " << n << ", eps=" << normaR1it[n] / normab << ", normaR1[i]=" << normaR1it[n]
							<< ", normaRes[i]=" << normaResit[n] << ", exRes[i]=" << exResit[n] << endl;
						//coutVector(x1, "x", 'H', 5); 
					}
				}
			}
		}
	}
	int n = normaResit.size() - 1;
	cout << "сошлось при nitr=" << n << ",  ||res|| = " << normaResit[n] << " ,  ||R1|| = " << normaR1it[n] << endl;
	slae->setX(x1);// , niter, normaR1it, normaResit, exResit)
}


vector<double> BiCGSTAB_CSR_OMP(const vector<double> &A, const vector<double> &b, const vector<double> &newx0, const vector<double> &RealSol,
	const vector<int> &iptr, const vector<int> &jptr, int &niter, double &eps, int &stopres,
	vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit, const int threadNum) {
	/*"""Программа принимает на вход систему Ax=b,
	причем А = квадратная матрица NxN элементов, хранимая в CSR формате.
	Программа возвращает вектор решений х;
	номер итерации на которой провесс расчет остановился;
	вектор содержащий норму невязок :*/
	cout << "BiCGSTAB Solution start" << endl;
	vector <double> x0 = newx0;
	int N = b.size();
	cout << "N" << N << endl;
	//	coutVector(A, "A", 'H');
	//	coutVector(b, "b", 'H');

	vector<double>	r0(N, 0.), Ax0(N, 0.);
	CSRmultiMV_OMP(A, x0, iptr, jptr, threadNum, Ax0);
	ax_plus_by_OMP(1., b, -1., Ax0, r0, threadNum);

	//	coutVector(CSRmultiMV(A, x0, iptr, jptr), "CSRmultiMV(A, x0, iptr, jptr)", 'H');
	//	coutVector(r0, "r0", 'H');
	vector<double>	_r(N, 1.),
		p0 = r0,
		x1(N, 0.),
		Ap(N, 0.),
		s(N, 0.),
		As(N, 0.),
		res(N, 0.),
		r1(N, 0.);
	double	normab = vector_norma_OMP(b, threadNum),
		r0_r = 0.,
		Apr = 0.,
		alfa = 0.,
		sqrAs = 0.,
		w = 0.,
		betta = 0.;

	for (int i = 0; i < niter; i++) {
		CSRmultiMV_OMP(A, p0, iptr, jptr, threadNum, Ap);
		//	coutVector(Ap, "Ap", 'H');
		//	cout << Ap.size() << "  " << endl;
		r0_r = inner_product_OMP(r0, _r, threadNum);
		//	cout << "r0_r "<< r0_r << endl;
		Apr = inner_product_OMP(Ap, _r, threadNum);
		//	cout << "Apr " << Apr << endl;
		if (Apr == 0)
		{
			cout << "devision on (Ap,r) = 0)" << endl;
			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x0;// , i - 1, normaR1it, normaResit, exResit)
		}
		alfa = r0_r / Apr;
		ax_plus_by_OMP(1., r0, -1., alfa*Ap, s, threadNum);
		CSRmultiMV_OMP(A, s, iptr, jptr, threadNum, As);
		sqrAs = inner_product_OMP(As, As, threadNum);
		if (sqrAs <= 10e-20) {
			//cout << "devision on (As,As) = 0)" << endl;
			//cout << "сошлось при ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x0;// , i - 1, normaR1it, normaResit, exResit)
		}
		w = inner_product_OMP(As, s, threadNum) / sqrAs;
		ax_plus_by_plus_cz_OMP(alfa, p0, w, s, 1., x0, x1, threadNum);
		vector<double> Ax1(N, 0.);
		CSRmultiMV_OMP(A, x1, iptr, jptr, threadNum, Ax1);
		ax_plus_by_OMP(1., b, -1., Ax1, res, threadNum);
		if (i >= stopres) { r1 = res; }
		else { ax_plus_by_OMP(1., s, -w, As, r1, threadNum); }
		normaR1it.push_back(vectorNorma(r1));
		normaResit.push_back(vectorNorma(res));
		//	std:cout << "RealSol.size() " << RealSol.size() << "  x1.size() " << x1.size() << endl;
		vector<double>dExres(N, 0.0);
		ax_plus_by_OMP(1., RealSol, -1., x1, dExres, threadNum);
		exResit.push_back(vector_norma_OMP(dExres, threadNum));
		if (normaR1it[i] / normab <= eps) {
			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x1;// , i, normaR1it, normaResit, exResit)
		}
		if (r0_r == 0) {
			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x1;// , i, normaR1it, normaResit, exResit)
		}
		if (w == 0) {
			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x1;// , i, normaR1it, normaResit, exResit)
		}
		betta = (inner_product_OMP(r1, _r, threadNum) / r0_r)*(alfa / w);
		ax_plus_by_plus_cz_OMP(1., r1, betta, p0, -betta * w, Ap, p0, threadNum);
		r0 = r1;
		x0 = x1;

	}
	int n = normaResit.size() - 1;
	cout << "сошлось при nitr=" << n << ",  ||res|| = " << normaResit[n] << " ,  ||R1|| = " << normaR1it[n] << endl;
	return x1;// , niter, normaR1it, normaResit, exResit)
}

//bool BiCGSTAB_CSR_MKL(int N, const double *A, const double *b, const double *newx0, const double *RealSol,
//	const int *iptr, const int *jptr, int &niter, double &eps, int &stopres,
//	vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit, const int threadNum, double* result) {
//	/*"""Программа принимает на вход систему Ax=b,
//	причем А = квадратная матрица NxN элементов, хранимая в CSR формате.
//	Программа возвращает вектор решений х;
//	номер итерации на которой провесс расчет остановился;
//	вектор содержащий норму невязок :*/
//	cout << "BiCGSTAB Solution start" << endl;
//
//	double *x0 = new double[N, 0.];
//	_dcopy(N, newx0, x0);
//
//	double*	r0 = new double[N, 0.], *Ax0 = new double[N, 0.];
//	//CSRmultiMV_(A, x0, iptr, jptr, threadNum, Ax0);
//	mul_matrix(N, A, iptr, jptr, x0, Ax0);
//	_daxpby(N, -1., Ax0, 1, 1., b, r0);
//
//	//	coutVector(CSRmultiMV(A, x0, iptr, jptr), "CSRmultiMV(A, x0, iptr, jptr)", 'H');
//	//	coutVector(r0, "r0", 'H');
//	double*	_r = new double[N, 1.],
//		*p0 = r0,
//		*Ap = new double[N, 0.],
//		*s = new double[N, 0.],
//		*r1 = new double[N, 0.],
//		*As = new double[N, 0.],
//		*res = new double[N, 0.];
//	double	normab = _dnrm2(N, b),
//		r0_r = 0.,
//		Apr = 0.,
//		alfa = 0.,
//		sqrAs = 0.,
//		w = 0.,
//		betta = 0.;
//
//	for (int i = 0; i < niter; i++) {
//		mul_matrix(N, A, iptr, jptr, p0, Ap);
//		//	coutVector(Ap, "Ap", 'H');
//		//	cout << Ap.size() << "  " << endl;
//		r0_r = _ddotc(N, r0, _r);
//		//	cout << "r0_r "<< r0_r << endl;
//		Apr = _ddotc(N, Ap, _r);
//		//	cout << "Apr " << Apr << endl;
//		if (Apr == 0)
//		{
//			cout << "devision on (Ap,r) = 0)" << endl;
//			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
//			return x0;// , i - 1, normaR1it, normaResit, exResit)
//		}
//		alfa = r0_r / Apr;
//		_daxpby(N, 1., r0, 1, -alfa, Ap, s);
//		mul_matrix(N, A, iptr, jptr, s, As);
//		sqrAs = _ddotc(N, As, As);
//		if (sqrAs <= 10e-20) {
//			//cout << "devision on (As,As) = 0)" << endl;
//			//cout << "сошлось при ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
//			return 1;// , i - 1, normaR1it, normaResit, exResit)
//		}
//		w = _ddotc(N, As, s) / sqrAs;
//		double *x1 = new double[N, 0.0];
//		_daxpbypcz(N, alfa, p0, 1, w, s, 1., x0, x1);
//		double* Ax1 = new double[N, 0.];
//		mul_matrix(N, A, iptr, jptr, x1, Ax1);
//		_daxpby(N, 1., b, 1, -1., Ax1, res);
//		if (i >= stopres) { r1 = res; }
//		else { _daxpby(N, 1., s, 1, -w, As, r1); }
//		normaR1it.push_back(_dnrm2(N, r1));
//		normaResit.push_back(_dnrm2(N, res));
//		//	std:cout << "RealSol.size() " << RealSol.size() << "  x1.size() " << x1.size() << endl;
//		double *dExres = new double[N, 0.0];
//		_daxpby(N, 1., RealSol, 1, -1., x1, dExres);
//		exResit.push_back(_dnrm2(N, dExres));
//		if (normaR1it[i] / normab <= eps) {
//			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
//			return 1;// , i, normaR1it, normaResit, exResit)
//		}
//		if (r0_r == 0) {
//			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
//			return 1;// , i, normaR1it, normaResit, exResit)
//		}
//		if (w == 0) {
//			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
//			return 1;// , i, normaR1it, normaResit, exResit)
//		}
//		betta = (_ddotc(N, r1, _r) / r0_r)*(alfa / w);
//		_daxpbypcz(N, 1., r1, 1, betta, p0, -(betta*w), Ap, p0);
//		r0 = r1;
//		x0 = x1;
//
//	}
//	int n = normaResit.size() - 1;
//	cout << "сошлось при nitr=" << n << ",  ||res|| = " << normaResit[n] << " ,  ||R1|| = " << normaR1it[n] << endl;
//	return 1;// , niter, normaR1it, normaResit, exResit)
//}

vector<double> BiCGSTAB_Yacobi_CSR(const vector<double> &A, const vector<double> &b, const vector<double> &newx0, const vector<double> &RealSol,
	const vector<int> &iptr, const vector<int> &jptr, int &niter, double &eps, int &stopres,
	vector<double> &normaR1it, vector<double> &normaResit, vector<double> &exResit) {
	/*"""Программа принимает на вход систему Ax=b,
	причем А = квадратная матрица NxN элементов, хранимая в CSR формате.
	Программа возвращает вектор решений х;
	номер итерации на которой провесс расчет остановился;
	вектор содержащий норму невязок :*/
	cout << "BiCGSTAB + Yacobi preconditioning Solution start" << endl;
	int rank = newx0.size();
	vector<double> r0(rank);
	r0 = b - CSRmultiMV(A, newx0, iptr, jptr);

	vector<double> _r = r0,
		s0 = r0,
		p0(rank, 0.),
		v0(rank, 0.),
		_p0(rank),
		x1(rank),
		p1(rank),
		_p1(rank),
		v1(rank),
		s1(rank),
		q(rank),
		r1(rank),
		u(rank),
		res;
	// Complex number
	//vector<double> ..;
	// Real number
	double ro1 = 0., ro0 = 1., alfa0 = 1., w0 = 1.,
		betta = 0., alfa1 = 0., w1 = 0., v_r0 = 0.;
	double normab = innerProduct(b, b);
	//Yacobi_precondition(rank, A, iptr, jptr, p0, _p0);
	x1 = newx0 + r0;
	r1 = b - CSRmultiMV(A, x1, iptr, jptr);

	//int niter_stop = 1000;
	for (size_t i = 1; i < niter; i++)
	{
		ro1 = innerProduct(r1, _r);
		//	coutVector(r1, "r1", 'H');
		if (ro0 == 0.) {
			cout << "devision on ro1= 0)" << endl;
			cout << "сошлось при ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x1;// , i - 1, normaR1it, normaResit, exResit)
		}
		if (w0 == 0.) {
			cout << "devision on w0= 0)" << endl;
			cout << "сошлось при ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x1;// , i - 1, normaR1it, normaResit, exResit)}
		}
		if (alfa0 == 0.) {
			cout << "devision on alfa0= 0)" << endl;
			cout << "сошлось при ||res|| = " << normaResit[i] << " ,  ||R1|| = " << normaR1it[i] << endl;
			return x1;// , i - 1, normaR1it, normaResit, exResit)}
		}
		cout << "ro1=" << ro1 << " ro0=" << ro0 << " alfa0=" << alfa0 << "  w0 =" << w0 << endl;
		betta = (ro1 / ro0) / (alfa0 / w0);
		cout << "betta " << betta << endl;
		p1 = r1 - betta * (p0 - w0 * v0); //p1= r1, i = 0
									  //	coutVector(p1, "p1", 'H');
		Yacobi_precondition(rank, A, iptr, jptr, p1, _p1);
		//	coutVector(_p1, "_p1", 'H');
		v1 = CSRmultiMV(A, _p1, iptr, jptr);
		//	coutVector(v1, "v1", 'H');

		v_r0 = (innerProduct(v1, _r));
		if (v_r0 == 0.) {
			cout << "devision on v_r0= 0)" << endl;
			cout << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i - 1] << " ,  ||R1|| = " << normaR1it[i - 1] << endl;
			return x1;
		}
		alfa1 = ro1 * (1. / v_r0);//Почему  алфа так быстро и сильно растет????

		s1 = r1 - alfa1 * v1;
		Yacobi_precondition(rank, A, iptr, jptr, s1, q);
		//	coutVector(q, "q", 'H');
		u = CSRmultiMV(A, q, iptr, jptr);
		double u_sqr = innerProduct(u, u);
		if (u_sqr == 0.) {
			cout << "devision on u_sqr= 0)" << endl;
			cout << "сошлось при nitr=" << i << ", ||res|| = " << normaResit[i - 1] << " ,  ||R1|| = " << normaR1it[i - 1] << endl;
			return x1;
		}
		w1 = innerProduct(u, s1) *(1. / u_sqr);
		x1 = x1 + (alfa1*_p1 + w1 * q);
		//coutVector(x1, "solution", 'H');
		//coutVector(s1, "s1", 'H');
		//coutVector(u, "u", 'H');
		cout << "w1 " << w1 << endl;
		r1 = s1 - w1 * u;
		res = b - CSRmultiMV(A, x1, iptr, jptr);

		normaR1it.push_back(vectorNorma(r1));
		normaResit.push_back(vectorNorma(res));
		exResit.push_back(vectorNorma(RealSol - x1));
		cout << "normaR1[i-1] " << normaR1it[i - 1] << endl;
		if (normaR1it[i - 1] / normab <= eps) {
			cout << scientific << setprecision(2) << "сошлось при nitr=" << i << ",  ||res|| = " << normaResit[i - 1] << " ,  ||R1|| = "
				<< normaR1it[i - 1] << ", ||err|| = " << exResit[i - 1] << endl;
			return x1;
		}
		alfa0 = alfa1;
		p0 = p1;
		w0 = w1;
		v0 = v1;
		ro0 = ro1;
	}

	int n = normaResit.size() - 1;
	cout << scientific << setprecision(2) << "сошлось при ||res|| = " << normaResit[n] << " ,  ||R1|| = " << normaR1it[n] << ", ||err|| = " << exResit[n] << endl;
	return x1;// , niter, normaR1it, normaResit, exResit)
}


// ************************************************************************************************************************************************************************
// *********************************************************************** Other Solvers **********************************************************************************
// ************************************************************************************************************************************************************************

//метод  Простой Итерации
vector<vector<double>> SimpleIter(const vector<vector<double>> &A, const vector<double> &f, const  vector<double> &x0, double eps)
{//  в первой строке матрицы Х (со вскеми решениями на каждой итерации) и в первых ячейках res, error будут лежать нулевые прилижения, заданные вручную
 // поэтому цикл будет начинаться с 1, а не с 0
	int rank = A.size();
	vector<vector<double>> B, x;
	vector<double> res, normares, error;
	int count = 1;
	error.push_back(eps + 0.1);
	normares.push_back(eps + 0.1);
	x.push_back(x0);
	B = A;
	for (int i = 0; i < rank; i++)
	{
		B[i][i] -= 1.0;
	}
	B = B*(-1);
	while (error[count - 1] > eps)
	{
		x.push_back(innerProduct(B, x[count - 1]) + f);
		res = innerProduct(A,x[count]) + f*(-1);
		coutVector(res, "res", 'H');
		normares.push_back(vectorNorma(res));
		coutVector(x[count], "xk", 'H');
		coutVector(x[count - 1], "xk-1", 'H');
		error.push_back(vectorNorma(x[count] + x[count - 1] * (-1)));
		//	cout<<error.size() << endl;
		cout << count << "  error = " << error[count] << "  res = " << res[count + 1] << '\n';
		count++;
		if (count == 10)
			break;
	}
	x.push_back(normares);
	x.push_back(error);
	return(x);
}

//41метод гаусса прямой и обратный ход
vector<double> Gauss_straight(vector<vector<double>> &A,  vector<double> &b)
{
	vector<double> h(int(A.size()), 0);
	//vector<vector<double>> result = A;
	for (int i = 0; i < int(A.size()) - 1; i++)
	{
		//		cout << "i = " << i << endl;
		for (int k = i + 1; k < int(A[i].size()); k++)
		{

			double a = -A[k][i] / A[i][i];
			b[k] += b[i] * a;
			//cout << "i= " << i << "  k=" << k << "  a = " << a <<"   b = "<<b[i]<< endl;
			for (int j = i; j < int(A[i].size()); j++)
			{
				A[k][j] += A[i][j] * a;
			}
			//coutVector(result[k], 2);
		}
	}
	int rank = A.size();
	h[rank - 1] = b[rank - 1] / A[rank - 1][rank - 1];
	//cout << "hend=" << h[rank - 1] << endl;
	for (int i = rank - 2; i >= 0; i--)
	{
		double r = 0;
		for (int j = i + 1; j < int(A[i].size()); j++)
		{
			r += h[j] * A[i][j];
			//cout << "r=" << r << endl;
		}
		h[i] = (b[i] - r) / A[i][i];
		//cout << "i= " << i << " h = " << h[i] << endl;
	}
	return h;
}
//42
void GaussSolver(vector<vector<double>> &SLAU,  vector<double> &right, vector<double> &futureresult)
{
	//ifstream f;//Объявляем поток чтения из файла
	//vector<vector<double>> SLAU;// , SLAUt;
	//SLAU = inputf_vector2d(f, SLAUf);
	//vector<double> right_part;// (70, 0), b(70, 0);
	//right_part = inputf_vectord(f, rightf);
	futureresult = Gauss_straight(SLAU, right);
	//coutVector(solve1, "solve1", 'H');
	/*coutVector(right_part, 2);
	SLAU.push_back({ 1,2,3 });
	SLAU.push_back({ 4,5,6});
	SLAU.push_back({ 1,2,5 })*/
	//coutVector(SLAU, "SLAU");
	double norm = 5;
	double div = 0;
	double lumbda = 1;

	//dLU_decomposition(SLAU, Lmatrix, Umatrix);
	//cout << "L-matrix" << endl;
	//coutVector(Lmatrix,3);
	//cout << "U-matrix" << endl;
	//coutVector(Umatrix, 3);
	//cout << "L*U-matrix" << endl;
	//coutVector(multiMatrix(Lmatrix, Umatrix),2);
	//dLU_solver(Lmatrix, Umatrix, calc_result, right_part);
	//coutVector(calc_result, 2);
	//calc_of_all(x, y, calc_result, mass_array, n, end_period, st_period, q_period, start_point, end_point, lumbda, h);
	//outputf_vectors(futureresult, "solution_15_12.txt");
}


// ************************************************************************************************************************************************************************
// ********************************************************************* Preconditioners **********************************************************************************
// ************************************************************************************************************************************************************************
void diag_normalize(const int &rank, vector<double> &SLAU,
	const vector<int> &iptr, const  vector<int> &jptr, vector<double> &right) {
	//(int nb, int *idi, double *df_1, double *v, double *w)
	int di_n = 0;
	double di_elm = 1.;
	for (int i = 0; i < rank; i++)
	{
		int st_str = iptr[i],
			end_str = iptr[i + 1];	// начало строки
									//SLAU[p] = help;
									//cout << "start = " << st_str << "  end = " << end_str << endl;
		di_n = findElement(jptr, st_str, end_str, i);
		cout << "SLAU[" << di_n << "  " << SLAU[di_n] << endl;
		di_elm = SLAU[di_n];
		//	cout << "p = " << p << "  int(p) = " << int(p) << "  nelem = " << nelem << endl;

		right[i] = right[i] * (1. / abs(di_elm));
		for (int j = st_str; j < end_str; j++)
		{
			//		cout << "befor SLAU[j] " << SLAU[j] << endl;
			SLAU[j] = SLAU[j] * (1. / abs(di_elm));
			//	cout << "after SLAU[j] " << SLAU[j] << endl;
		}

	}

}

void Yacobi_precondition(const int &rank, const vector<double> &SLAU,
	const vector<int> &iptr, const  vector<int> &jptr, vector<double> &in, vector<double> &out) {
	//(int nb, int *idi, double *df_1, double *v, double *w)
	// один раз решаем систему SLAU*out = in, те просто делим in на диагональные элементы
	int di_n = 0;
	double di_elm = 1.;
	for (int i = 0; i < rank; i++)
	{
		int st_str = iptr[i],
			end_str = iptr[i + 1];	// начало строки
									//SLAU[p] = help;
									//cout << "start = " << st_str << "  end = " << end_str << endl;
		di_n = findElement(jptr, st_str, end_str, i);
		//	cout << "SLAU[" << di_n << "  " << SLAU[di_n] << endl;
		di_elm = SLAU[di_n];
		//	cout << "p = " << p << "  int(p) = " << int(p) << "  nelem = " << nelem << endl;
		out[i] = in[i] * (1. / abs(di_elm));
		//		coutVector(out, "out", 'H');

	}

}

void Ematrix(const int& rank, vector<vector<double>> &E) {
	E = vector<vector<double>>(rank, vector<double>(rank, 0.));
	for (size_t i = 0; i < rank; i++)
	{
		E[i][i] = 1.;
	}
}

double matrixDeterminante(const vector<vector<double>> &A) { return 0.; }
double matrixDiagDeterminante(const vector<vector<double>> &A) {
	double det = 1.;
	for (size_t i = 0; i < A.size(); i++)
	{
		det *= A[i][i];
	}
	return det;
}
void reverseMatrixByIterShultze(const vector<vector<double>> &A, vector<vector<double>> &revA) {
// https://calcsbox.com/post/iteracionnyj-metod-sulca-nahozdenia-obratnoj-matricy.html
	int rank = A.size(), k = 5;
	double eps = 0.0001;
	vector<vector<double>> E, Qi, Qi2, Adot_revA, QiSumm, // m = 2
		revAnext, revAprev = vector<vector<double>>(rank, vector<double>(rank, 0.1));
	Ematrix(rank, E);
	
	for (size_t i = 0; i < k; i++)
	{
		
		innerProduct(A, revAprev, Adot_revA);
		coutVector(Adot_revA, "Adot_revA");
		Qi = E - Adot_revA;
		coutVector(Qi, "Qi");
		innerProduct(Qi, Qi, Qi2);
		coutVector(Qi2, "Qi2");
		QiSumm = E + Qi + Qi2;
		coutVector(QiSumm, "Qisumm");
		innerProduct(revAprev, QiSumm, revAnext);
		coutVector(revAnext, "revAi");
		revAprev = revAnext;
	}
}