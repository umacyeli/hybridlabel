#include<string>
#include<vector>
#include<set>
#include<map>
#include<unordered_set>
#include<unordered_map>
#include<queue>

using namespace std;

class Graph{
private:
	vector<vector<int> > v;//in label if the graph is undirected, only in label is used
	//vector<int> > u;//out label
	vector<vector<int> > vw;//weight for in label if the graph is unweighted, all weight are set to 1
	vector<vector<int> > uw;//weight for out label

	vector<vector<int> > adj;
	vector<vector<int> > adjw;

	vector<pair<int, int> > degree;
	vector<pair<int, int> > edges;


	int n;//the number of vertices
	int weighted;//0: unweighted 1: weighted
	int directed;//0: undirected 1: directed

	vector<set<int> > hneighbors;
	vector<set<int> > hreach;
	vector<int> drank;//drank[v] = rank. get rank of a vertice v
	vector<int> dorder;//dorder[r] = v. get vertice by its rank.
	vector<vector<int> > label;
	vector<vector<int> > labeld;
	vector<vector<int> > nlabel;//links+shortcuts
	vector<vector<int> > nlabeld;
	vector<vector<int> > hidden_label;
	vector<vector<int> > hidden_labeld;
	vector<vector<int> > glabel;//guidences
	vector<vector<int> > glabeld;
	vector<vector<int> > olabel;
	vector<vector<int> > olabeld;


	vector<int> com;
	vector<vector<pair<int, int> > > plabel;
	vector<set<int> > pset;
	int olabel_amount;	

	string graphFile;//path to the given graph file
	
	vector<unordered_set<int> > P;
	vector<unordered_set<int> > border_nodes;
	vector<unordered_set<string> > border_edges;
	string hash_edge(pair<int, int> se);

public:
	int rv;
	int V;
	Graph(string graphFile, int weighted, int directed);
	int read();
	int loadLabel(const char* filename);
	int loadRankOrder(const char* filename);
	int compressLabel(const char* filename);
	int decomposeLabel();
	int loadNewlabel(const char* filename);
	int QueryNewLabel(int u, int v);
	int lightenLabel(int v, int num);
	int query(int u, int v);
	int pllquery(int u, int v);
	/*
	int triade(double ep);
	int init(unordered_set<int> &dirty);	
	unordered_set<int> intersect(int v, int u);
	int clean(unordered_set<int> &dirty, double ep);
	double jaccard(pair<int, int> se);
	int deleteEdge(pair<int, int> se);
	int extract(unordered_set<int> &R, unordered_set<int> &bn, unordered_set<string> &be, unordered_set<int> &dirty);
	*/
};
