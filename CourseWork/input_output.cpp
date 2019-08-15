#include "pch.h"
// ФУНКЦИИ ВЫВОДА И ВВОДА ДАННЫХ В ФАЙЛ

#include"input_output.h"


void outputf_vectors(const vector<vector<double>> &A, const char namefile[100])
{
	ofstream fout1(namefile); // создаём объект класса ofstream для записи и связываем его с файлом cppstudio.txt
//	fout1 << A.size() << "  " << A[0].size()<<"\n";
	for (int i = 0; i < int(A.size()); i++)
	{
		//cout << "A[" << i << "].size() = " << A[i].size() << endl;
		//fout1 << i << "  ";
		for (int j = 0; j < int(A[i].size()); j++)
		{
			fout1 << setw(15) << A[i][j] << " ";
		}
		fout1 << "\n";
	}
	fout1.close(); // закрываем файл
}
void outputf_vectors(const vector<double> &A,  const char namefile[100])
{
	ofstream fout1(namefile); // создаём объект класса ofstream для записи и связываем его с файлом cppstudio.txt
//	fout1 << A.size() <<"\n";
	for (int i = 0; i < int(A.size()); i++)
	{
		fout1 << setw(15) << A[i] << endl;
	}
	fout1.close(); // закрываем файл
}
void outputf_vectors(const vector<vector<int>> &A, const char namefile[100])
{
	ofstream fout1(namefile); // создаём объект класса ofstream для записи и связываем его с файлом cppstudio.txt
	fout1 << A.size() << "  " << A[0].size()<<"\n";

	for (int i = 0; i < int(A.size()); i++)
	{
		
		//cout << "A[" << i << "].size() = " << A[i].size() << endl;
		for (int j = 0; j < int(A[i].size()); j++)
		{
			fout1 << setw(8) << A[i][j] << " ";
		}
		fout1 << "\n";
	}
	fout1.close(); // закрываем файл
}
void outputf_vectors(const vector<int> &A, const char namefile[100])
{
	ofstream fout1(namefile); // создаём объект класса ofstream для записи и связываем его с файлом cppstudio.txt
	fout1 << A.size() << endl;
	for (int i = 0; i < int(A.size()); i++)
	{
		fout1 << setw(15) << A[i] << endl;
	}
	fout1.close(); // закрываем файл
}
void outputf_sets(const set<int> &A, const char namefile[100])
{
	ofstream fout1(namefile); // создаём объект класса ofstream для записи и связываем его с файлом cppstudio.txt
	fout1 << A.size() << endl;
	for (set<int>::iterator it = A.begin(); it != A.end(); ++it)
	{
		fout1 << setw(5) << *it << endl;
	}
	fout1.close(); // закрываем файл
}
void outputf(ofstream &f, const double &a, const char s[100])//Функция для записи													 //числа в файл
{
	f.open(s);//Открываем файл
			  //Проверка успешности открытия файла:
	if (f.fail()) {
		cout << "\n Ошибка открытия файла";
		exit(1);
	}
	f << a;//записываем переменную в файл
	f.close();//Закрываем файл
}
vector<vector<vector<double>>> inputf_vector3d(ifstream &f,  char s[100])
{//Объявление указателей и переменных
	int N, M, L;
	//Выделение ОП под строки и столбцы матрицы
	vector<vector<vector<double>>> result;
	//Открытие файла для чтения
	f.open(s);
	//Проверка успешности открытия файла для чтения
	if (f.fail()) {
		cout << "\n Ошибка открытия файла";
		exit(1);
	}
	f >> N;
	f >> M;
	f >> L;
	for (int i = 0; i< N; i++)
	{
		vector<double> b(L, 0);
		vector<vector<double>> c((N), b);
		for (int j = 0; j < N; j++)
		{
			for (int l = 0; l< L; l++)
			{
				f >> c[j][l];
			}
		}
		result.push_back(c);
	}
	//Закрытие файла
	f.close();
	return result;
}
vector<vector<double>> inputf_vector2d(ifstream &f,  const char s[100])
{
	f.open(s);
	if (f.fail()) {
		cout << "\n Ошибка открытия файла";
		exit(1);
	}
	int N, M;
	f >> N;
	f >> M;
	vector<double> b(N, 0);
	vector<vector<double>> result((M), b);

	for (int i = 0; i< N; i++)
	{

		for (int j = 0; j < M; j++)
		{
			f >> result[i][j];
		}
	}
	//Закрытие файла
	f.close();
	return result;
}
vector<vector<int>> inputf_vector2i(ifstream &f,  const char s[100])
{

	f.open(s);
	if (f.fail()) {
		cout << "\n Ошибка открытия файла";
		exit(1);
	}
	int N, M;
	f >> N;
	f >> M;
	vector<int> b(M, 0);
	vector<vector<int>> result(N, b);

	for (int i = 0; i< N; i++)
	{

		for (int j = 0; j < M; j++)
		{
			f >> result[i][j];
		}
	}
	f.close();
	return result;
}
vector<double> inputf_vectord(ifstream &f, const char s[100])
{
	f.open(s);
	int N;
	//N = inputf(f,s);
	f >> N;
	cout << N << endl;
	vector<double> result(N, 0);
	if (f.fail()) {
		cout << "\n Ошибка открытия файла";
		exit(1);
	}
	for (int i = 0; i<N; i++)
	{
		f >> result[i];
	}
	f.close();
	return result;
}
vector<int> inputf_vectori(ifstream &f, const char s[100])
{
	f.open(s);
	int N;
	//N = inputf(f,s);
	f >> N;
	cout << N << endl;
	vector<int> result(N, 0);
	if (f.fail()) {
		cout << "\n Ошибка открытия файла";
		exit(1);
	}
	for (int i = 0; i<N; i++)
	{
		f >> result[i];
	}
	f.close();
	return result;
}
double inputf(ifstream &f, const char s[100])//Функция для чтения числа из файла
{
	double a;
	f.open(s);//Открываем файл
			  //Проверка успешности открытия файла:
	if (f.fail()) {
		cout << "\n Ошибка открытия файла";
		exit(1);
	}
	f >> a;//Читаем переменную из файла
	f.close();//Закрываем файл
	return a;
}

void outputf_vectors(const vector<vector<double>> &A, const wstring namefile, const wstring path)
{
	std::experimental::filesystem::create_directories(path);
	// outputf(, path, L"test.txt");
	std::ofstream f;

	f.open(path + L"\\" + namefile);//Открываем файл
			  //Проверка успешности открытия файла:
	if (f.fail()) {
		std::cout << "\n Ошибка открытия файла";
		exit(1);
	}

//	fout1 << A.size() << "  " << A[0].size()<<"\n";
	for (int i = 0; i < int(A.size()); i++)
	{
		//cout << "A[" << i << "].size() = " << A[i].size() << endl;
		//fout1 << i << "  ";
		for (int j = 0; j < int(A[i].size()); j++)
		{
			f << setw(8) << A[i][j] << " ";
		}
		f << "\n";
	}
	f.close(); // закрываем файл
}
void outputf_vectors(const vector<double> &A, const wstring &namefile, const wstring &path, const wstring newCatalog)
{
	std::experimental::filesystem::create_directories(path+ L"\\"+newCatalog);
	// outputf(, path, L"test.txt");
	std::ofstream f;

	f.open(path +L"\\" + newCatalog + L"\\" + namefile+L".txt");//Открываем файл
			  //Проверка успешности открытия файла:
	if (f.fail()) {
		std::cout << "\n Ошибка открытия файла";
		exit(1);
	}
	//	fout1 << A.size() <<"\n";
	for (int i = 0; i < int(A.size()); i++)
	{
		f << setw(15)<< setprecision(15)<< A[i] << endl;
	}
	f.close(); // закрываем файл
}

void outputf_vectors(const vector<string> &A, const wstring &namefile, const wstring &path, const wstring newCatalog)
{
	std::experimental::filesystem::create_directories(path + L"\\" + newCatalog);
	// outputf(, path, L"test.txt");
	std::ofstream f;

	f.open(path + L"\\" + newCatalog + L"\\" + namefile+L".txt");//Открываем файл
			  //Проверка успешности открытия файла:
	if (f.fail()) {
		std::cout << "\n Ошибка открытия файла";
		exit(1);
	}
	//	fout1 << A.size() <<"\n";
	for (int i = 0; i < int(A.size()); i++)
	{
		f << A[i] << endl;
	}
	f.close(); // закрываем файл
}

void outputf_vectors(const vector<vector<int>> &A, const wstring namefile, const wstring path)
{
	std::experimental::filesystem::create_directories(path);
	// outputf(, path, L"test.txt");
	std::ofstream f;

	f.open(path + L"\\" + namefile);//Открываем файл
			  //Проверка успешности открытия файла:
	if (f.fail()) {
		std::cout << "\n Ошибка открытия файла";
		exit(1);
	}

	for (int i = 0; i < int(A.size()); i++)
	{

		//cout << "A[" << i << "].size() = " << A[i].size() << endl;
		for (int j = 0; j < int(A[i].size()); j++)
		{
			f << setw(8) << A[i][j] << " ";
		}
		f << "\n";
	}
	f.close(); // закрываем файл
}
void outputf_vectors(const vector<int> &A, const wstring namefile, const wstring path)
{
	std::experimental::filesystem::create_directories(path);
	// outputf(, path, L"test.txt");
	std::ofstream f;

	f.open(path + L"\\" + namefile);//Открываем файл
			  //Проверка успешности открытия файла:
	if (f.fail()) {
		std::cout << "\n Ошибка открытия файла";
		exit(1);
	}
	f << A.size() << endl;
	for (int i = 0; i < int(A.size()); i++)
	{
		f << setw(15) << A[i] << endl;
	}
	f.close(); // закрываем файл
}
void outputf_sets(const set<int> &A, const wstring namefile, const  wstring path)
{
	std::experimental::filesystem::create_directories(path);
	// outputf(, path, L"test.txt");
	std::ofstream f;

	f.open(path + L"\\" + namefile);//Открываем файл
			  //Проверка успешности открытия файла:
	if (f.fail()) {
		std::cout << "\n Ошибка открытия файла";
		exit(1);
	}

	f << A.size() << endl;
	for (set<int>::iterator it = A.begin(); it != A.end(); ++it)
	{
		f << setw(5) << *it << endl;
	}
	f.close(); // закрываем файл
}
void outputf(ofstream &f, const double &a, const wstring s, const wstring path)//Функция для записи													 //числа в файл
{
	
	// outputf(, path, L"test.txt");
	
	f.open(path + L"\\" + s);//Открываем файл
			  //Проверка успешности открытия файла:
	if (f.fail()) {
		std::cout << "\n Ошибка открытия файла";
		exit(1);
	}

	f.open(s);//Открываем файл
			  //Проверка успешности открытия файла:
	if (f.fail()) {
		cout << "\n Ошибка открытия файла";
		exit(1);
	}
	f << a;//записываем переменную в файл
	f.close();//Закрываем файл
}

int outputf(double a, char s[100])//Функция для записи													 //числа в файл
{
	std::ofstream f;
	f.open(s);//Открываем файл
			  //Проверка успешности открытия файла:
	if (f.fail()) {
		std::cout << "\n Ошибка открытия файла";
		exit(1);
	}
	f << a;//записываем переменную в файл
	f.close();//Закрываем файл
	return 4;
}

void create_folder(double a, const wstring folder, const wstring file) {
	/* если дериктории нет, то создаем с заданным названием (в конце всегда есть _N,
	если папка с таким N уже есть, то в конце приписываем N+1, и тд*/
}
void outputf(double a, const wstring folder, const wstring file)//Функция для записи													 //числа в файл
{
	std::ofstream f;

	f.open(folder + L"\\" + file);//Открываем файл
			  //Проверка успешности открытия файла:
	if (f.fail()) {
		std::cout << "\n Ошибка открытия файла";
		exit(1);
	}
	f << a;//записываем переменную в файл
	f.close();//Закрываем файл
	//return 4;
}

void launch_out(){

		int result = 0;
		wstring path(L"new_folder");
		/* L"" - wide-char литерал, имеет тип const wchar_t[n], где n кол-во символов + '\0'.
		CreateDirectory является макросом, который в зависимости от конфигурации проекта развернется либо в CreateDirectoryA,
		если юникод выключен и в CreateDirectoryW, если юникод включен. Соотстветственно LPCTSTR разворачивается либо в char*, либо в wchar_t*.
		При явном вызове CreateDirectoryW соотстветственно туда должен передаваться wchar_t*, string просто так в wchar_t* не конвертируешь,
		однако можно использовать wstring.*/

	//	create_directories(path);
	//	outputf(6, path, L"test.txt");	
}


void read_gmesh_data(string gmsh_filename) {
	void gmsh_data_read(string gmsh_filename, int node_dim, int node_num,
		double node_x[], int element_order, int element_num, int element_node[]);

	void gmsh_size_read(string gmsh_filename, int &node_num, int &node_dim,
		int &element_num, int &element_order);

}
//bool CreatePath(std::wstring &wsPath)
//{
//	DWORD attr;
//	int pos;
//	bool result = true;
//
//	// Check for trailing slash:
//	pos = wsPath.find_last_of(SLASH);
//	if (wsPath.length() == pos + 1)  // last character is "\"
//	{
//		wsPath.resize(pos);
//	}
//
//	// Look for existing object:
//	attr = GetFileAttributesW(wsPath.c_str());
//	if (0xFFFFFFFF == attr)  // doesn't exist yet - create it!
//	{
//		pos = wsPath.find_last_of(SLASH);
//		if (0 < pos)
//		{
//			// Create parent dirs:
//			result = CreatePath(wsPath.substr(0, pos));
//		}
//		// Create node:
//		result = result && CreateDirectoryW(wsPath.c_str(), NULL);
//	}
//	else if (FILE_ATTRIBUTE_DIRECTORY != attr)
//	{  // object already exists, but is not a dir
//		SetLastError(ERROR_FILE_EXISTS);
//		result = false;
//	}
//
//	return result;
//}
// TODO: Установите ссылки на любые требующиеся дополнительные заголовки в файле STDAFX.H
// , а не в данном файле
