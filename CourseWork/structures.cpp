#include "pch.h"

void SLAE_CSR::coutSLAE(const string & name)
{
	// coutMatrixCSR()
	cout << name << endl;
	int rank = iptr.size() - 1;
	coutVector(iptr, "iptr", 'H');
	coutVector(jptr, "jptr", 'H');

	for (int i = 0; i < rank; i++)
	{
		cout << "i = " << i << ": ";
		for (int j = (int(iptr[i])); j<int(iptr[i + 1]); j++){

			cout << setprecision(2) << A[j] << ",  ";
		}
		cout << endl;
	}
	coutVector(b, "b:: ", 'H');
}

void Mesh::coutMesh()
{
	
	
	coutVector(X, "X", 'H');
	coutVector(Y, "Y", 'H');
	coutVector(elemNodes, "elemNodes");
	coutVector(elemBF, "elemBF");
	coutVector(nghbrs, "nghbrs");
	coutVector(bndrElms[0], "bndrElms[0]");
	coutVector(bndrElms[1], "bndrElms[1]");
	coutVector(bndrElms[2], "bndrElms[2]");
	coutVector(bndrElms[3], "bndrElms[3]");
	coutVector(elemBF, "elemBF");
}

void quadrElem::coutElem()
{
	cout << "x1 = " << x1 << ", x2 = " << x2 << ", y1 = " << y1 << ", y2 = " << y2 << endl;
	cout << "detJ = " << detJ << endl;
	coutVector(J, "Jacobian");

}
