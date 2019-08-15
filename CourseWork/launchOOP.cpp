#include"pch.h"
#include "launchOOP.h"

void createClassicSolutionforCG_MP1(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams)
{
	wstring path = L"ClassicSolutionforCG_MP1";// L"D:\source\repos\CourseWork2\CourseWork2\ClassicSolutionforCG_MP1";
	cout << "\n ClassicSolutionforCGplusSUPG_MP1()\n";
	int basisRank = 1;

	string basisItem = "lagrange", discrSol = "AdvDifSUPG";
	AbstractModel* model = new MP1biquadAS(indataGrid, indataParams);
	cout << sizeof(model) << endl;

	string br = to_string(basisRank);
	model->addDescription("Discrete model type::  " + discrSol);
	model->addDescription("basis::  " + basisItem + " " + br);
	model->coutModelParams();
	Mesh mesh = generateQuadElementMesh("CG", basisRank, indataGrid);
	//	mesh.coutMesh();

	SolutionCreater* creater = new SolutionCreaterCG();
	SolutionDiscreteAnalog* SDA = creater->createDiscreteAnalog(discrSol, basisItem, basisRank, model, &mesh);
	SDA->BoundCondDerichlet(1, -1);
	SDA->BoundCondDerichlet(2, -1);
	SDA->BoundCondDerichlet(3, -1);
	SDA->BoundCondDerichlet(4, -1);
	//	SDA->coutSDA();

	SLAE_CSR* slae1 = SDA->getSLAE();
	//	coutMatrixCSR(slae1, "slae");
	vector<double> x0(slae1->b.size(), 1.),
		normaR1It, normaResIt, exResIt;
	vector<double> realSol;
	model->getRealSol(mesh.X, mesh.Y, realSol);
	int niter = 1000, pivot = 50;
	double eps = 0.000000000001;
	// vector<double> GaussSolution = Gauss_straight(slae, right);
	SDA->addDescriptionToModel("Solver::  BiCGSTAB_CSR");
	SDA->addDescriptionToModel("input params: Niter:: " + to_string(niter));
	SDA->addDescriptionToModel("pivot:: " + to_string(pivot));
	SDA->addDescriptionToModel("eps:: " + to_string(eps));
	BiCGSTAB_CSR(slae1, x0, realSol, niter, eps, pivot, normaR1It, normaResIt, exResIt);
	coutVector(slae1->x, "sol", 'H');
	double L2normaError = SDA->L2Norma();
	cout << "L2 norma error = " << L2normaError << endl;
	SDA->addDescriptionToModel("Output params: Niter:: " + to_string(normaR1It.size()));
	SDA->addDescriptionToModel("L2 error norma:: " + to_string(L2normaError));
	SDA->addDescriptionToModel("R1 norma:: " + to_string(normaR1It.back()));
	SDA->addDescriptionToModel("Res norma:: " + to_string(normaResIt.back()));
	string newCatalogString = (discrSol+ " tau x120, eps=1e-12 " + basisItem + to_string(basisRank) + " Test_" + to_string(Ntest));
	wstring newCatalog(newCatalogString.begin(), newCatalogString.end());
	wcout << newCatalog << endl;
	outputfSolutionResults(path, newCatalog, realSol, slae1->x, indataGrid, mesh.X, mesh.Y, indataParams, SDA->getDescription(),
		normaR1It, normaResIt, exResIt);
	
}

void createDGSolutionforCG_MP1(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams)
{
	wstring path = L"DGSolutionforCG_MP1";// L"D:\source\repos\CourseWork2\CourseWork2\ClassicSolutionforCG_MP1";
	cout << "\n ClassicSolutionforCGplusSUPG_MP1()\n";
	int basisRank = 2;

	string basisItem = "lagrange", discrSol = "AdvDif";
	AbstractModel* model = new MP1biquadAS(indataGrid, indataParams);
	cout << sizeof(model) << endl;

	string br = to_string(basisRank);
	model->addDescription("Discrete model type::  " + discrSol);
	model->addDescription("basis::  " + basisItem + " " + br);
	model->coutModelParams();
	Mesh mesh = generateQuadElementMesh("DG", basisRank, indataGrid);
	mesh.coutMesh();

	SolutionCreater* creater = new SolutionCreaterDG();
	SolutionDiscreteAnalog* SDA = creater->createDiscreteAnalog(discrSol, basisItem, basisRank, model, &mesh);
	SDA->BoundCondDerichlet(1, -1);
	SDA->BoundCondDerichlet(2, -1);
	SDA->BoundCondDerichlet(3, -1);
	SDA->BoundCondDerichlet(4, -1);
	//	SDA->coutSDA();

	SLAE_CSR* slae1 = SDA->getSLAE();
	slae1->coutSLAE("slae");
	vector<double> x0(slae1->b.size(), 1.),
		normaR1It, normaResIt, exResIt;
	vector<double> realSol;
	model->getRealSol(mesh.X, mesh.Y, realSol);
	int niter = 1000, pivot = 50;
	double eps = 0.000000000001;
	// vector<double> GaussSolution = Gauss_straight(slae, right);
	SDA->addDescriptionToModel("Solver::  BiCGSTAB_CSR");
	SDA->addDescriptionToModel("input params: Niter:: " + to_string(niter));
	SDA->addDescriptionToModel("pivot:: " + to_string(pivot));
	SDA->addDescriptionToModel("eps:: " + to_string(eps));
	BiCGSTAB_CSR(slae1, x0, realSol, niter, eps, pivot, normaR1It, normaResIt, exResIt);
	coutVector(slae1->x, "sol", 'H');
	double L2normaError = SDA->L2Norma();
	cout << "L2 norma error = " << L2normaError << endl;
	SDA->addDescriptionToModel("Output params: Niter:: " + to_string(normaR1It.size()));
	SDA->addDescriptionToModel("L2 error norma:: " + to_string(L2normaError));
	SDA->addDescriptionToModel("R1 norma:: " + to_string(normaR1It.back()));
	SDA->addDescriptionToModel("Res norma:: " + to_string(normaResIt.back()));
	vector<double> pointSol;
	vector<double> newX, newY;
	vector<double>newIndatagrid = indataGrid;
	newIndatagrid[3] = indataGrid[3] * 3;
	newIndatagrid[8] = indataGrid[8] * 3;

	quadGridForQuadDomain(newIndatagrid, newX, newY);
	SDA->calcPointSolution(newX, newY, pointSol);
	
	string newCatalogString = (discrSol + " mu=4., eps=1e-12 " + basisItem + to_string(basisRank) + " Test_" + to_string(Ntest));
	wstring newCatalog(newCatalogString.begin(), newCatalogString.end());
	wcout << newCatalog << endl;
	outputfSolutionResults(path, newCatalog, realSol, pointSol, indataGrid, newX, newY, indataParams, SDA->getDescription(),
		normaR1It, normaResIt, exResIt);
}

void createClassicSolutionforMP2(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams)
{
	cout << "\ncreateClassicSLAEforMP1()\n";
	int basisRank = 1;
	string basisItem = "lagrange";
	AbstractModel* model = new MP2withBndryLayers(indataGrid, indataParams);
	// or: AbstractModel* model = new CG_MP1(indataGrid, indataParams);
	cout << sizeof(model) << endl;

	model->coutModelParams();
	Mesh mesh = generateQuadElementMesh("CG", basisRank, indataGrid);
	mesh.coutMesh();

	SolutionCreater* creater = new SolutionCreaterCG();
	SolutionDiscreteAnalog* SDA = creater->createDiscreteAnalog("AdvDif", basisItem, basisRank, model, &mesh);
	SDA->BoundCondDerichlet(1, -1);
	SDA->BoundCondDerichlet(2, -1);
	SDA->BoundCondDerichlet(3,-1);
	SDA->BoundCondDerichlet(4,-1);
	SDA->coutSDA();

	SLAE_CSR* slae1 = SDA->getSLAE();
//	coutMatrixCSR(slae1, "slae");
	vector<double> x0(slae1->b.size(), 1.),
		normaR1It, normaResIt, exResIt;
	vector<double> realSol;
	model->getRealSol(mesh.X, mesh.Y, realSol);
	int niter = 100, pivot = 11;
	double eps = 0.0000000000000001;
	// vector<double> GaussSolution = Gauss_straight(slae, right);
	BiCGSTAB_CSR(slae1, x0, realSol, niter, eps, pivot, normaR1It, normaResIt, exResIt);
	coutVector(slae1->x, "sol", 'H');
	double L2normaError = SDA->L2Norma();
	cout << "L2 norma error = " << L2normaError << endl;	
	
	
}

void createClassicSolutionforCG_MP4(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams)
{
	wstring path = L"ClassicSolutionforCG_MP4";// L"D:\source\repos\CourseWork2\CourseWork2\ClassicSolutionforCG_MP1";
	wcout << path << endl;
	int basisRank = 1;

	string basisItem = "lagrange", discrSol = "AdvDifSUPG";
	// AbstractModel* model = new MP1biquadAS();
	AbstractModel* model = new MP4withConstBC1(indataGrid, indataParams);
	string br = to_string(basisRank);
	model->addDescription("Discrete model type::  " + discrSol);
	model->addDescription("basis::  " + basisItem + " " + br);
	cout << sizeof(model) << endl;
	model->coutModelParams();
	Mesh mesh = generateQuadElementMesh("CG", basisRank, indataGrid);
	// mesh.coutMesh();

	SolutionCreater* creater = new SolutionCreaterCG();
	SolutionDiscreteAnalog* SDA = creater->createDiscreteAnalog(discrSol, basisItem, basisRank, model, &mesh);
	SDA->BoundCondDerichlet(1,-1);
	SDA->BoundCondDerichlet(2,-1);
	SDA->BoundCondDerichlet(3,-1);
	SDA->BoundCondDerichlet(4,-1);
//  SDA->coutSDA();

	SLAE_CSR* slae1 = SDA->getSLAE();
//	coutMatrixCSR(slae1, "slae");
	vector<double> x0(slae1->b.size(), 1.),
		normaR1It, normaResIt, exResIt;
	vector<double> realSol;
	model->getRealSol(mesh.X, mesh.Y, realSol);
	int niter = 100000, pivot = 100;
	double eps = 0.0000001;
	// vector<double> GaussSolution = Gauss_straight(slae, right);
	SDA->addDescriptionToModel("Solver::  BiCGSTAB_CSR");
	SDA->addDescriptionToModel("input params: Niter:: " + to_string(niter));
	SDA->addDescriptionToModel("pivot:: " + to_string(pivot));
	SDA->addDescriptionToModel("eps:: " + to_string(eps));
	BiCGSTAB_CSR(slae1, x0, realSol, niter, eps, pivot, normaR1It, normaResIt, exResIt);
	coutVector(slae1->x, "sol", 'H');
	
	SDA->addDescriptionToModel("Output params: Niter:: " + to_string(normaR1It.size()));
	// SDA->addDescriptionToModel("L2 error norma:: " + to_string(L2normaError));
	SDA->addDescriptionToModel("R1 norma:: " + to_string(normaR1It.back()));
	SDA->addDescriptionToModel("Res norma:: " + to_string(normaResIt.back()));
	// string newCatalogString = (discrSol + " tau x5 " + basisItem + to_string(basisRank) + " Test_" + to_string(Ntest));
	string newCatalogString = (discrSol + ", tau x6 " + basisItem + to_string(basisRank) + " Test_" + to_string(Ntest));
	wstring newCatalog(newCatalogString.begin(), newCatalogString.end());
	wcout << newCatalog << endl;
	
	outputfSolutionResults(path, newCatalog, realSol, slae1->x, indataGrid, mesh.X, mesh.Y, indataParams, SDA->getDescription(),
		normaR1It, normaResIt, exResIt);
}

void createCGSolutionforCG_MP7(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams)
{
	wstring path = L"ClassicSolutionforCG_MP7.2 adaptive mesh";// L"D:\source\repos\CourseWork2\CourseWork2\ClassicSolutionforCG_MP1";
	cout << "\n ClassicSolutionforCG_MP7()\n";
	int basisRank = 1;

	string basisItem = "lagrange", discrSol = "AdvDifSUPG";
	AbstractModel* model = new MP7withBndryLayer(indataGrid, indataParams);
	cout << sizeof(model) << endl;

	string br = to_string(basisRank);
	model->addDescription("Discrete model type::  " + discrSol);
	model->addDescription("basis::  " + basisItem + " " + br);
	model->coutModelParams();
	Mesh mesh = generateQuadElementMesh("CG", basisRank, indataGrid);
	mesh.coutMesh();	
	SolutionCreater* creater = new SolutionCreaterCG();
	SolutionDiscreteAnalog* SDA = creater->createDiscreteAnalog(discrSol, basisItem, basisRank, model, &mesh);
	SDA->BoundCondDerichlet(1,-1);
	SDA->BoundCondDerichlet(2,-1);
	SDA->BoundCondDerichlet(3,-1);
	SDA->BoundCondDerichlet(4,-1);
	//	SDA->coutSDA();

	SLAE_CSR* slae1 = SDA->getSLAE();
	//	coutMatrixCSR(slae1, "slae");
	vector<double> x0(slae1->b.size(), 1.),
		normaR1It, normaResIt, exResIt;
	vector<double> realSol;
	model->getRealSol(mesh.X, mesh.Y, realSol);
	int niter = 1000, pivot = 50;
	double eps = 0.000000001;
	// vector<double> GaussSolution = Gauss_straight(slae, right);
	SDA->addDescriptionToModel("Solver::  BiCGSTAB_CSR");
	SDA->addDescriptionToModel("input params: Niter:: " + to_string(niter));
	SDA->addDescriptionToModel("pivot:: " + to_string(pivot));
	SDA->addDescriptionToModel("eps:: " + to_string(eps));
	BiCGSTAB_CSR(slae1, x0, realSol, niter, eps, pivot, normaR1It, normaResIt, exResIt);
	coutVector(slae1->x, "sol", 'H');
	double L2normaError = SDA->L2Norma();
	cout << "L2 norma error = " << L2normaError << endl;
	SDA->addDescriptionToModel("Output params: Niter:: " + to_string(normaR1It.size()));
	SDA->addDescriptionToModel("L2 error norma:: " + to_string(L2normaError));
	SDA->addDescriptionToModel("R1 norma:: " + to_string(normaR1It.back()));
	SDA->addDescriptionToModel("Res norma:: " + to_string(normaResIt.back()));
	string newCatalogString = (discrSol + " tau x40 " + basisItem + to_string(basisRank) + " Test_" + to_string(Ntest));
	wstring newCatalog(newCatalogString.begin(), newCatalogString.end());
	wcout << newCatalog << endl;
	outputfSolutionResults(path, newCatalog, realSol, slae1->x, indataGrid, mesh.X, mesh.Y, indataParams, SDA->getDescription(),
		normaR1It, normaResIt, exResIt);
}

void createClassicSolutionforCG_MP8(const int Ntest, const vector<double>&indataGrid, const vector<double> &indataParams)
{

	wstring path = L"ClassicSolutionforCG_MP8";// L"D:\source\repos\CourseWork2\CourseWork2\ClassicSolutionforCG_MP1";
	wcout << path << endl;
	int basisRank = 1;

	string basisItem = "lagrange", discrSol = "AdvDifSUPG";
	AbstractModel* model = new MP8constSource(indataGrid, indataParams);
	cout << sizeof(model) << endl;

	string br = to_string(basisRank);
	model->addDescription(basisItem);
	model->addDescription(br);
	model->coutModelParams();
	Mesh mesh = generateQuadElementMesh("CG", basisRank, indataGrid);
	// mesh.coutMesh();

	SolutionCreater* creater = new SolutionCreaterCG();
	SolutionDiscreteAnalog* SDA = creater->createDiscreteAnalog(discrSol, basisItem, basisRank, model, &mesh);
	SDA->BoundCondDerichlet(1,-1);
	SDA->BoundCondDerichlet(2,-1);
	SDA->BoundCondDerichlet(3,-1);
	SDA->BoundCondDerichlet(4,-1);
	//	SDA->coutSDA();

	SLAE_CSR* slae1 = SDA->getSLAE();
	slae1->coutSLAE("full slae");
	vector<double> x0(slae1->b.size(), 0.),
		normaR1It, normaResIt, exResIt;
	vector<double> realSol;
	model->getRealSol(mesh.X, mesh.Y, realSol);
	int niter = 1000, pivot = 50;
	double eps = 0.0000000000000001;
	// vector<double> GaussSolution = Gauss_straight(slae, right);
	BiCGSTAB_CSR(slae1, x0, realSol, niter, eps, pivot, normaR1It, normaResIt, exResIt);
	coutVector(slae1->x, "sol", 'H');
	//double L2normaError = SDA->L2Norma();
	//cout << "L2 norma error = " << L2normaError << endl;

	string newCatalogString = (discrSol+basisItem + to_string(basisRank) + "Test_" + to_string(Ntest));
	wstring newCatalog(newCatalogString.begin(), newCatalogString.end());
	wcout << newCatalog << endl;
	outputfSolutionResults(path, newCatalog, realSol, slae1->x, indataGrid, mesh.X, mesh.Y, indataParams, SDA->getDescription(),
		normaR1It, normaResIt, exResIt);
}


void createClassicSolutionforCG_RP1(const int Ntest, const vector<double>&indataGrid,
	const vector<double> &indataParams, const vector<double> &timeParams)
{
	wstring path = L"ClassicSolutionforCG_RP1+source";// L"D:\source\repos\CourseWork2\CourseWork2\ClassicSolutionforCG_MP1";
	int basisRank = 1;

	string basisItem = "lagrange", discrSol = "AdvDif";
	AbstractModel* model = new RP1_filterInPorousMedia(indataGrid, indataParams);
	cout << sizeof(model) << endl;

	string br = to_string(basisRank);
	model->addDescription("Discrete model type::  " + discrSol);
	model->addDescription("basis::  " + basisItem + " " + br);
	model->coutModelParams();
	Mesh mesh = generateQuadElementMesh("CG", basisRank, indataGrid);
	mesh.coutMesh();

	TimeDepSolutionCreaterCG* creater = new TimeDepSolutionCreaterCG();
	SolutionDiscreteAnalog* SDA = creater->createDiscreteAnalog(discrSol, basisItem, basisRank, model, &mesh, timeParams);
	SDA->BoundCondDerichlet(3, 0);
	SDA->BoundCondDerichlet(4, 0);
	//	SDA->coutSDA();

	SLAE_CSR* slae1 = SDA->getSLAE();
	//	coutMatrixCSR(slae1, "slae");
	vector<double> x0(slae1->b.size(), 1.),
		normaR1It, normaResIt, exResIt;
	vector<double> realSol;
	model->getRealSol(mesh.X, mesh.Y, realSol);
	int niter = 100, pivot = 50;
	double eps = 0.000000000001;
	// vector<double> GaussSolution = Gauss_straight(slae, right);
	SDA->addDescriptionToModel("Solver::  BiCGSTAB_CSR");
	SDA->addDescriptionToModel("input params: Niter:: " + to_string(niter));
	SDA->addDescriptionToModel("pivot:: " + to_string(pivot));
	SDA->addDescriptionToModel("eps:: " + to_string(eps));
	BiCGSTAB_CSR(slae1, x0, realSol, niter, eps, pivot, normaR1It, normaResIt, exResIt);
	coutVector(slae1->x, "sol", 'H');
		SDA->addDescriptionToModel("Output params: Niter:: " + to_string(normaR1It.size()));
	SDA->addDescriptionToModel("L2 error norma::  Null ");
	SDA->addDescriptionToModel("R1 norma:: " + to_string(normaR1It.back()));
	SDA->addDescriptionToModel("Res norma:: " + to_string(normaResIt.back()));
	vector<double> prevU = slae1->b;

	string addDiscript = " taux15, eps=1e-12", timeCat = "\\T=" + to_string(0);
	string newCatalogString = (discrSol + addDiscript + basisItem + to_string(basisRank)
		+ " Test_" + to_string(Ntest));
	wstring mainCatalog(newCatalogString.begin(), newCatalogString.end()),
		timeStepCatalog (timeCat.begin(), timeCat.end());
	
	wcout << mainCatalog << endl;
	double timeStep = (timeParams[1] - timeParams[0]) / timeParams[2];
	outputfSolutionResults(path, mainCatalog+timeStepCatalog, realSol, slae1->x, indataGrid, mesh.X, mesh.Y, indataParams, SDA->getDescription(),
		normaR1It, normaResIt, exResIt);
	for ( int i = 1; i < timeParams[2]; i++)
	{
		cout << "time step " << i << endl;
		SDA->calcRightVectorNewTimeStep(prevU);
		SDA->BoundCondDerichlet(3, i);
		SDA->BoundCondDerichlet(4, i);
		timeCat = "\\T=" + to_string(i);
		timeStepCatalog = wstring(timeCat.begin(), timeCat.end());
		slae1 = SDA->getSLAE();
		BiCGSTAB_CSR(slae1, x0, realSol, niter, eps, pivot, normaR1It, normaResIt, exResIt);
		outputfSolutionResults(path, mainCatalog+timeStepCatalog, realSol, slae1->x, indataGrid, mesh.X, mesh.Y, indataParams, SDA->getDescription(),
			normaR1It, normaResIt, exResIt);
	}
}


void outputfSolutionResults(const wstring path, const wstring newCatalogName,
	const vector<double> &analytSol, const vector<double> &solution, 
	const vector<double> &indataMesh, const vector<double> &meshX, const vector<double> &meshY,
	const vector<double> modelParams, const vector<string>& discription,
	const vector<double> &normaR1, const vector<double> &normaRes, const vector<double> &normaExRes ){
	
	wstring filenameSol = L"Solution",
		filename_paramsFEM = L"Parameters",
		filename_titleFEM = L"Title_FEM",
		filename_normaR1 = L"NormaR1",
		filename_normaRes = L"NormaRes",
		filename_normaERROR = L"NormaERROR",
		filename_RealSol = L"RealSol",
		filename_X = L"X",
		filename_Y = L"Y",
		filename_discrptn = L"SolDiscription",
		filename_mesh = L"indata";

	int EndIter = normaR1.size();
	//EndIterYacobi = normaR1itYacobi.size(), EndIterSUPGYacobi = normaR1itSUPGYacobi.size();
	wofstream fout(filename_titleFEM);
	//{
	//	fout << left << newCatalogName <<
	//		/*"\n Mesh " << meshX.size() << "x" << meshY.size() << " elms \n";
	//	fout << left << scientific << setprecision(1) << "s=" << sigma << "; " << "g=" << fixed << setprecision(0) << gamma << "; " <<
	//		"U={" << Utest[0] << ", " << Utest[1] << "}, Pe=" << Pecle << "\n";
	//	fout << left << "it=" << EndIter << ", Res=" << scientific << setprecision(1) << normaResit[EndIter - 1]
	//		<< ", R1=" << normaR1it[EndIter - 1] << ", Err=" << exResit[EndIter] <<*/ "\n";
	//}
	//fout.close(); // закрываем файл

	//fout.open(filename_paramsFEM);
	///*{
	//	fout << left << ntest << endl <<
	//		charsize << endl <<
	//		setprecision(1) << sigma << endl
	//		<< gamma << endl
	//		<< Utest[0] << endl
	//		<< Utest[1] << endl
	//		<< Pecle << "\n"
	//		<< EndIter << endl;

	//	if (normaRes[EndIter - 1] != -NAN)
	//	{
	//		fout << normaResit[EndIter - 1] << endl;
	//	}
	//	else { fout << -1 << endl; }
	//	if (normaR1it[EndIter - 1] != -NAN)
	//	{
	//		fout << normaR1it[EndIter - 1] << endl;
	//	}
	//	else { fout << -1 << endl; }
	//	if (exResit[EndIter - 1] != -NAN)
	//	{
	//		fout << exResit[EndIter - 1] << endl;
	//	}
	//	else { fout << -1 << endl; }

	//}*/
	//fout.close(); // закрываем файл

	outputf_vectors(modelParams, filename_paramsFEM, path, newCatalogName);
	outputf_vectors(solution, filenameSol, path, newCatalogName);
	outputf_vectors(analytSol, filename_RealSol, path, newCatalogName);
	outputf_vectors(normaR1, filename_normaR1, path, newCatalogName);
	outputf_vectors(normaRes, filename_normaRes, path, newCatalogName);
	outputf_vectors(normaExRes, filename_normaERROR, path, newCatalogName);
	outputf_vectors(meshX, filename_X, path, newCatalogName);
	outputf_vectors(meshY, filename_Y, path, newCatalogName);
	outputf_vectors(discription, filename_discrptn, path, newCatalogName);
	outputf_vectors(indataMesh, filename_mesh, path, newCatalogName);
	}