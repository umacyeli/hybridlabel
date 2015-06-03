#include "Graph.h"
#include<fstream>
#include<algorithm>
#include<sstream>
#include<iostream>

using namespace std;

Graph::Graph(string graphFile, int weighted, int directed){
	this->graphFile = graphFile;
	this->weighted = weighted;
	this->directed = directed;
}

int Graph::loadRankOrder(const char* filename){
	ifstream ifs(filename);
	if(!ifs.good()) return 0;
	ifs >> V;
	drank.resize(V);
	dorder.resize(V);
	int r,v;
	while(ifs >> r >> v){
		drank[v] = r;
		dorder[r] = v;
	}
	ifs.close();	
	return 1;
}

int Graph::loadLabel(const char* filename){
	label.clear();
	labeld.clear();
	label.resize(V);
	labeld.resize(V);
	hneighbors.resize(V);
	hreach.resize(V);	

	ifstream ifs(filename);
	ifs >> V;
	if(!ifs.good()) return 0;
	int v;
	int ln;
	while(ifs >> v >> ln){
		vector<int>& label_v = label[v];
		vector<int>& labeld_v = labeld[v];
		set<int>& hneighbors_v = hneighbors[v];
		int u,d;
		for(int i = 0; i < ln; i++){
			ifs >> u >> d;	
			label_v.push_back(u);
			labeld_v.push_back(d);
			if( d == 1)
				hneighbors_v.insert(u);
		}
		if( v == V-1){
			ifs >> olabel_amount;
			break;
		}
	}

	ifs.close();
	return 1;	
}

int Graph::loadNewlabel(const char* filename){
	nlabel.clear();
	nlabeld.clear();
	hidden_label.clear();
	hidden_labeld.clear();
	nlabel.resize(V);
	nlabeld.resize(V);
	hidden_label.resize(V);
	hidden_labeld.resize(V);	

	ifstream ifs(filename);
	ifs >> V;
	if(!ifs.good()) return 0;
	int v;
	int ln;
	while(ifs >> v >> ln){
		int u,d;
	       for(int i = 0; i < ln; i++){
			ifs >> u >> d;
			if(d == 1){
				nlabel[v].push_back(u);
				nlabeld[v].push_back(d);
			}else if(d > 1){
				hidden_label[v].push_back(u);
				hidden_labeld[v].push_back(d);
			}
	       }		       
	       if( v == V-1)
		       break;
	}		

	ifs.close();
	return 1;

}

int Graph::QueryNewLabel(int u, int v){
	queue<int> qu, qv;	
	queue<int> dqu, dqv;
	int du = 0, dv = 0;
	map<int, int> md, mi;
	int maxu, maxv;
	int lb = INT_MAX, up = -1;

	qu.push(u);
	qv.push(v);
	dqu.push(du);
	dqv.push(dv);
	md[u] = du;
	md[v] = dv;
	mi[u] = 1;
	mi[v] = 1;
	maxu = drank[u];
	maxv = drank[v];

	for(int i = 0; i < hidden_label[u].size(); i++){
		int w = hidden_label[u][i];
		int wd = hidden_labeld[u][i];
		if(mi[w] != 0){
			md[w] += wd;
			mi[w] = 3;
			int mdw = md[w];
			if(mdw < lb)
				lb = mdw;
		}else{
			md[w] = wd;
			mi[w] = 3;
		}
	}
	for(int i = 0; i < hidden_label[v].size(); i++){
		int w = hidden_label[v][i];
		int wd = hidden_labeld[v][i];
		if(mi[w] != 0){
			md[w] += wd;
			mi[w] = 3;
			int mdw = md[w];
			if(mdw < lb)
				lb = mdw;
		}else{
			md[w] = wd;
			mi[w] = 3;
		}
	}
	
	while(qu.size() !=0 || qv.size() != 0){
		queue<int>& fq,sq,dfq,fsq;
		int& maxf, maxs;
		if(maxu < maxv){
			fq = qv;
			sq = qu;
			dfq = dqv;
			dsq = dqu;
			maxf = maxv;
			maxs = maxu;
		}else{
			fq = qu;
			sq = qv;
			dfq = dqu;
			dsq = dqv;
			maxf = maxu;
			maxs = maxv;
		}
		int w = fq.front();
		int dw = dfq.front() + 1; 
		fq.pop();
		dfq.pop();
		for(int i = 0; i < nlabel[w].size(); i++){
			int ww = nalebl[w][i];
			if(mi[ww]==1){
				md[ww] += dw;
				return md[ww];
			}else if(mi[ww] >2){
				md[ww] += dw;
				if(md[ww] < lb)
					lb = md[ww];	
			}else{
				mi[ww] = 1;
				md[ww] = dw;
			}
			fq.push(ww);
			dfq.push(dw);
			if(drank[ww] < maxf)
				maxf = drank[ww];	
		}	
		w = sq.front();
		dw = dsq.front() + 1;
		sq.pop();
		dsq.pop();

	}	

}

int Graph::compresslabel(const char* filename){
	ofstream ofs(filename);
	ofs << V << endl;
	cout << "Start compressing!" << endl;
	int reduce_count = 0;
	for(int i = 0; i < V; i++){
		int v = dorder[i];//get the ith rank vertice
		vector<int> new_label;
		vector<int> new_labeld;
		for(auto& u : hneighbors[v]){
			hreach[v].insert(u);//can reach its higher rank neighbors
			for(auto& q : hreach[u]){
				hreach[v].insert(q);//can reach 
			}
		}
		set<int>::iterator it;
		for(int j = 0; j < label[v].size(); j++){
			if(labeld[v][j] > 1){
				it = hreach[v].find(label[v][j]);		
				if( it != hreach[v].end())
					reduce_count++;
				else{
					new_label.push_back(label[v][j]);
					new_labeld.push_back(labeld[v][j]);
				}
			}else{
				new_label.push_back(label[v][j]);
				new_labeld.push_back(labeld[v][j]);
			}
		}

		ofs << v << " " << new_label.size() << " "; 
		for(int j = 0; j < new_label.size(); j++){
			ofs << new_label[j] << " " << new_labeld[j] << " ";
		}
		ofs << endl;

		if(i%5000 == 0){
			cout << "progress: " << i << "/" << V << endl;
			cout << label[v].size() << " : " << hreach[v].size() << endl;
		}
		/*
		if(i == 5000){
			cout << drank[v] << endl;
			for(auto& u : hreach[v])
				cout << drank[u] << " ";
			cout << endl;
			for(int i = 0; i < label[v].size(); i++)
				cout << drank[label[v][i]] << "," << labeld[v][i] << " ";
			cout <<endl;
		}*/
	}
	cout << reduce_count << "/" << olabel_amount << endl;
	return 1;
	ofs.close();
}

//int query(int v, int u){
//	if( v >= V || u >=V ) return v == u ? 0 : INT_MAX;
//}

int Graph::read(){
	//map<pair<int, int>, int> weights;

	ifstream ifs(graphFile);	
	if(!ifs.good()) return 0;
	int v,u,w;
	long s;
	n = -1;
	
	//while(ifs >> v >> u >> w >> s){	
	while(ifs >> v >> u){	
		if(v > n)
			n = v;
		if(u > n)
			n = u;
		--v;--u;

		if(u != v)
			edges.push_back(make_pair(u, v));
		//valid_edges[se] = 1;
		//weights[se] = w;	
	}
	
	adj.resize(n);
	
	for(int i = 0; i < edges.size(); i++){
		adj[edges[i].first].push_back(edges[i].second);
		if(directed == 0)
			adj[edges[i].second].push_back(edges[i].first);
	}

	degree.resize(n);
	int max_degree = -1;
	for(int i = 0; i < n; i++){//O(V*d*logd)
		degree[i] = make_pair( adj[i].size(), i );
	}
	
	sort(degree.rbegin(), degree.rend());	

	vector<int> order(n);
	for(int i = 0; i < n; i++){
		order[degree[i].second] = i;
	}
	
	set<int> top_1_triangle;
	set<int> top_1_neighbors;
	map<int, int> top_1_triangle_count;
	set<int>::iterator it;
	int topv = n * 0.02;
	
	
	set<int> top_v_set;
	for(int i = 0; i < topv; i++){
		top_v_set.insert(degree[i].second);	
	}

	set<int> top_1_2hop;
	for(int i = 0; i < adj[degree[0].second].size(); i++){
		u = adj[degree[0].second][i];
		top_1_2hop.insert(u);
		for(int j = 0; j < adj[u].size(); j++){
			top_1_2hop.insert(adj[u][j]);
		}
	}
		
	/*
	v = degree[0].second;
	top_1_triangle.insert(v);
	for(int i = 0; i < adj[v].size(); i++){
		int u1 = adj[v][i];			
		top_1_neighbors.insert(u1);
		
		//2-hop neighbors
		//for(int j = 0; j < adj[u1].size(); j++)
		//	top_1_neighbors.insert(adj[u1][j]);

		
		set<int> v_n;
		set<int> u1_n;
		for(int j = 0; j < adj[v].size(); j++) v_n.insert(adj[v][j]);
		for(int j = 0; j < adj[u1].size(); j++) u1_n.insert(adj[u1][j]);
		for(auto& u2 : u1_n){
			it = v_n.find(u2);
			if(it != v_n.end()){//the intersection
				set<int> u2_n;
				for(int j = 0; j < adj[u2].size(); j++) u2_n.insert(adj[u2][j]);
				for(auto& uj : u2_n){
					it = u1_n.find(uj);
					if(it != u2_n.end()){
						top_1_triangle_count[uj]++;
					}
				}
			
			}
		}
	}
	
	vector<pair<int, int> > sort_t;
	for(auto& x : top_1_triangle_count){
		sort_t.push_back( make_pair(x.second, x.first) );
	}
	sort(sort_t.rbegin(), sort_t.rend());
	
	for(int i = 0; i < sort_t.size() && i < n* 0.02; i++){
		top_1_triangle.insert(sort_t[i].second);
	}
	


	int hit_count = 0;
	for(int i = 0; i < topv; i++){
		u = degree[i].second;
		
		it = top_1_triangle.find(u);
		if(it != top_1_triangle.end() )
			hit_count++;
		
		//it = top_1_neighbors.find(u);
		//if(it != top_1_neighbors.end() )
		//	hit_count++;
	}	

	cout << hit_count << "/" << topv << "," << n << "," << top_1_triangle.size() << "," << top_1_neighbors.size() << endl;
	*/
	return 1;
}
/*
int Graph::triade(double ep){
	unordered_set<int> dirty;
	init(dirty);	
	clean(dirty, ep);
	while(edges.size()!=0){
		unordered_set<int> R;
		unordered_set<int> bn;

		unordered_set<string> be;
		
		extract(R, bn, be, dirty);	
		
		P.push_back(R);
		border_nodes.push_back(bn);
		border_edges.push_back(be);

		clean(dirty, ep);	
	}	
}

int Graph::init(unordered_set<int> &dirty){
	for(auto& e: edges){//O(E*d) or O(|E| + w) according to the paper
		int v, u;
		istringstream iss(e);
		iss >> v >> u;

		unordered_set<int> iset = intersect(v, u);//O(d)

		T[e] = iset.size();
		//J[e] = (double)T[e] / (double)(degree[v] + degree[i] - T[e]);	
		dirty.insert(v);
		dirty.insert(u);
	}	
	return 1;
}

unordered_set<int> Graph::intersect(int v, int u){//O(d)
	unordered_set<int> iset;
	vector<int> nv = adj[v];
	vector<int> nu = adj[u];	
	
	//merge join
	for(int i = 0, j = 0; i < nv.size()&& j < nu.size() ; ){
		int vv = nv[i];
		int uu = nu[j];
		if(vv == uu){
			i++;
			j++;
			
			pair<int, int> ve, ue;
			if(vv > v) ve = make_pair(v, vv);
			else ve = make_pair(vv, v);	
			if(uu > u) ue = make_pair(u, uu);
			else ue = make_pair(uu, u);	
		
			unordered_set<string>::const_iterator gotve = edges.find(hash_edge(ve));
			unordered_set<string>::const_iterator gotue = edges.find(hash_edge(ue));
			if( gotve == edges.end() || gotue == edges.end() )
				continue;	

			iset.insert(vv);
		}else if(vv > uu)
			j++;
		else//vv < uu
			i++;
	}

	return iset;
}

int Graph::clean(unordered_set<int> &dirty, double ep){
	while(dirty.size()!=0){
		for(unordered_set<int>::iterator it = dirty.begin(); it != dirty.end(); ++it){
			int v = *it;
			for(int i = 0; i < adj[v].size(); i++){
				int u = adj[v][i];	
				pair<int, int> se;
				if(v > u) se = make_pair(u, v);
				else se = make_pair(v, u);	
				unordered_set<string>::const_iterator got = edges.find(hash_edge(se));
				if( got != edges.end() ){
					if(jaccard(se) < ep){
						deleteEdge(se);//O(1)
						dirty.insert(v);//O(1)
						dirty.insert(u);//O(1)
						break;
					}
				}
			}	
			dirty.erase(it);//will it cause a problem?
		}	
	}	
	return 1;
}

double Graph::jaccard(pair<int, int> se){//O(1)
	int v = se.first;
	int u = se.second;
	return (double)T[hash_edge(se)]/(double)(degree[v] + degree[u] - T[hash_edge(se)]);//O(1)
}

int Graph::deleteEdge(pair<int, int> se){//O(1)
	int v = se.first;
	int u = se.second;

   unordered_set<int> iset = intersect(v, u);//O(d)
       for(unordered_set<int>::iterator it = iset.begin(); it != iset.end(); it++){
	               int k = *it;
		               
		               pair<int, int> kv;
			               if(v > k) kv = make_pair(k, v);
				               else kv = make_pair(v, k);
					               
					               T[hash_edge(kv)]--;
						               
						               pair<int, int> ku;
							               if(u > k) ku = make_pair(k, u);
								               else ku = make_pair(u, k);
									           
									               T[hash_edge(ku)]--;
										           }
           

	edges.erase(hash_edge(se));	//O(1)
	A[degree[v]].erase(v);//O(1)
	A[degree[u]].erase(u);//O(1)

	degree[v]--;
	degree[u]--;
	A[degree[v]].insert(v);
	A[degree[u]].insert(u);

	return 1;
}

int Graph::extract(unordered_set<int> &R, unordered_set<int> &bn, unordered_set<string> &be, unordered_set<int> &dirty){
	int v = -1;
	unordered_map<int, int> theta;
	for(int i = A.size() - 1; i !=0 ; i--){ 	//O(d)?
		if(A[i].size() != 0)
			v = *(A[i].begin());//v has the largest degree
	}
	for(int i = 0; i < adj[v].size(); i++){//O(d^3) ?!!!!!
		int u1 = adj[v][i];
		unordered_set<int> iset = intersect(u1, v);//O(d)
		for(unordered_set<int>::iterator it = iset.begin(); it != iset.end(); it++){
			int u2 = *it;
			unordered_set<int> jset = intersect(u1, u2);
			for(unordered_set<int>::iterator jit = jset.begin(); jit != jset.end(); jit++){
				theta[*jit]++;
			}
		}
	}
	
	vector<pair<int, int> > sorted_theta;
       for(unordered_map<int, int>::iterator it = theta.begin(); it != theta.end(); it++)
       sorted_theta.push_back(make_pair(it->second, it->first));	       
	sort(sorted_theta.rbegin(), sorted_theta.rend());

	for(int i = 0; i < sorted_theta.size() && i < degree[v]; i++){
		R.insert(sorted_theta[i].second);
	}

	for(int i = 0; i < adj[v].size(); i++){
		int u = adj[v][i];

		pair<int, int> se;
		if(v > u) se = make_pair(u, v);
		else se = make_pair(v, u);	
		
		unordered_set<string>::const_iterator got = edges.find(hash_edge(se));
		if( got != edges.end() )
			R.insert(u);
	}

	for(unordered_set<int>::iterator it = R.begin(); it != R.end(); it++){//line 11
		int vv = *it;
		for(int i = 0; i < adj[vv].size(); i++){
			int uu = adj[vv][i];		
			pair<int, int> se;
			if(vv > uu) se = make_pair(uu, vv);
			else se = make_pair(vv, uu);	
			unordered_set<string>::const_iterator got = edges.find(hash_edge(se));
			if( got != edges.end() ){
				dirty.insert(uu);
				
				unordered_set<int>::const_iterator gotuu = R.find(uu);	
				if(gotuu == R.end()){//vv is a border node && se is a border edge
					bn.insert(vv);
					be.insert(hash_edge(se));	
				}

				//should delete edges here
				deleteEdge(se);
			}
		}
	}
	
	for(unordered_set<int>::iterator it = R.begin(); it != R.end(); it++){//line 11
		int vv = *it;
		dirty.erase(vv);//will it cause problem?
	}

	return 1;
}

string Graph::hash_edge(pair<int, int> se){
	string l = to_string(se.first);
	string r = to_string(se.second);
	string key = l + "\t" + r;
	return key;	
}
*/
