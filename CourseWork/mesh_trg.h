#pragma once
class mesh_trg
{
public:
	mesh_trg();
	~mesh_trg();
	
	/*string gmsh_filename, int node_dim, int node_num,
	double node_x[], int element_order,
	int element_num, int element_node[]);
*/
		// Main entieties
//	Mainpoints; //nodes - формирующие область
//	Maincurves; // 
//	mainsurfaces;

//	nodes;
//	triangles;
//	boundaries;
//	curves;





}; 

void bndry_elms(const vector<vector<int>>& trgs, const vector<int> &bndry_edges,
	vector<vector<int>>& nghbr_list, vector<short> &fool);

void find_nghbours(const vector<vector<int>>& trgs, const vector<int> &bndry_edges, vector<vector<int>>& nghbr_list);

int decoder(vector<bool> input);

vector<vector<bool>> check_edges(const vector<int>& cur, vector<int> other);

