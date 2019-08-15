#include "pch.h"
#include "mesh_trg.h"


mesh_trg::mesh_trg()
{

}


mesh_trg::~mesh_trg()
{

}

void bndry_elms(const vector<vector<int>>& trgs, const vector<int> &bndry_edges,
	vector<vector<int>>& nghbr_list, vector<short> &fool) {
	
	int Nedges = bndry_edges.size();
	vector<vector<bool>> checked_pair;
	int Nedge_i = 0;
	int Nelem = trgs.size();

	for (int i = 0; i < Nelem; i++){
		for (int j = i + 1; j < Nedges && fool[i]<= 2; j++){
			checked_pair = check_edges(trgs[i], trgs[j]);
			Nedge_i = decoder(checked_pair[0]);
			if (Nedge_i == 0) nghbr_list[i][0] = 1; // 1 point lies on boundary
			if (Nedge_i != -1) {
				nghbr_list[j][Nedge_i] = -1; // нет соседа для этой грани элемента, тк она принадлежит границе
				nghbr_list[i][0] = 2;// 2 points lie on boundary
				fool[j]++;
			}
		}
	}
}

void find_nghbours(const vector<vector<int>>& trgs, const vector<int> &bndry_edges, vector<vector<int>>& nghbr_list){
	// в nghbr_list в первой ячейке каждого элемента будет храниться инф-ия о принадлжености/непринадлжености его к границе.
	int Ncur_trg = 0, Nelem = trgs.size();
	vector<short> fool(Nelem, 0);
	nghbr_list = vector<vector<int>>(Nelem, vector<int>(4, 0));
	bndry_elms(trgs, bndry_edges, nghbr_list, fool);
	vector<vector<bool>> checked_pair;
	int Nedge_i = 0, Nedge_j = 0;
	for (int i = 0; i < Nelem && fool[i]!=3; i++)
	{
		for (int j = i+1; j < Nelem && fool[j] != 3; j++)
		{
			checked_pair = check_edges(trgs[i], trgs[j]);
			Nedge_i = decoder(checked_pair[0]);
			if (Nedge_i != -1) {
				Nedge_j = decoder(checked_pair[1]);
				nghbr_list[i][Nedge_i] = j;
				nghbr_list[j][Nedge_j] = i;
				fool[j]++;
				fool[i]++;
			}
			
		}
	}
	
}

int decoder(vector<bool> input) {
	if (vector<bool>{1, 1, 0} == input) return 1;
	if (vector<bool>{0, 1, 1} == input) return 2;
	if (vector<bool>{1, 0, 1} == input) return 3;
	if (vector<bool>{0, 0, 0} == input) return -1;// нет общих ребер
	else return 0; // нет общих ребер, но есть общая вершина

}

vector<vector<bool>> check_edges(const vector<int>& cur, vector<int> other)
{
	vector< vector<bool>> result(2, vector<bool>(3,0));
	// для треугольников с лин БФ
	for (int i = 1; i < cur.size(); i++)	{
		for (int j = 1; j < other.size(); j++){
			result[0][i-1] = result[1][j-1] = (cur[i] == other[j]);
		}
	}
	return result;
}
