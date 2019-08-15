#pragma once
class SolutionDiscreteAnalog
{
protected:
	// если struct - добавляем спецификатор const
	// если class - не нужно, тк все члены protected/private и доступ к ним можно получить только через getter
	AbstractModel* model;
	Mesh* mesh;
	SLAE_CSR slae;
	basisFunction* basis;
	SDAFactory* sdaFactory;
	double T0, Tn;
	double timeStep;
	int NtimeSteps;
	
public:
	SolutionDiscreteAnalog();
	// SolutionDiscreteAnalog(const string &item);
	SolutionDiscreteAnalog(AbstractModel *newModel, Mesh *newMesh, SDAFactory *newDAfitches, const vector<double>& timeParams);
	SolutionDiscreteAnalog(AbstractModel *newModel, Mesh *newMesh, SDAFactory*  newDAfitches);

	virtual ~SolutionDiscreteAnalog();
	
	// здесь будем определять какой базис использовать, пока только для квадратов.
	// 	virtual SLAE* createSLAE()=0;
	virtual void defineSolutionProcess(const string &basisItem, const int& rankBasis) = 0;
		// сборка SLAE
	virtual void assembleSlae(const int &nTimeStep) = 0;
	
	virtual void calcLocalMatrices(const currentParamsValue &curParams, const quadrElem &curElem,
		 vector<vector<double>> &resultA,  vector<double> &result) = 0;// какое уравнение - конве-диффузии, конв-дифф-реак или др

	virtual void BoundCondDerichlet(const int &nBoundary, const int &nTimeStep) = 0;

	virtual	void calcRightVectorNewTimeStep(const vector<double>& prevStepU);

	quadrElem calcLocalElem(int nElem);
	void calcPointSolution(const vector<double>&X, const vector<double>&Y, vector<double> &pointSol);
	double solutionInPoint(const double& x, const double& y);

	SLAE_CSR* getSLAE();
	Mesh* getMesh();
	AbstractModel* getModel();
	vector<string> getDescription();
	
	void coutSDA();
	double L2Norma();
	
	void addDescriptionToModel(const string &newLine);

	
protected:
	void calcLocalMassMatrix(const double &curGamma, const quadrElem &curElem,
		const vector<vector<double>> &FI, vector<vector<double>> &resultM);

	void calcLocalStiffnessMatrix(const double &curLyambda, const quadrElem &curElem, 
		const vector<vector<vector<double>>> &gradFI, vector<vector<double>> &resultG);
	
	void calcLocalAdvectiveMatrix(const vector<double> &curVelocity, const quadrElem &curElem,
		const vector<vector<double>> &FI, const vector<vector<vector<double>>> &gradFI,
		vector<vector<double>> &resultK);
	
	void calcLocalRightVector(const quadrElem &curElem, const vector<vector<double>> &FI, vector<double> &localb);
	
	virtual void calcLocalTimeMatrix( const quadrElem & curElem, const vector<vector<double>> &FI, vector<vector<double>> &resultM);
	
	virtual void calcLocalTimeVector(const quadrElem & curElem, const vector<double>& curF,
		const vector<vector<double>>& FI, vector<double>& resultB);

	int findPointInROw(const double &p, const vector<double> &row);

	vector<double> calcPointSolutionOnElmInGP(int nelem);
};

