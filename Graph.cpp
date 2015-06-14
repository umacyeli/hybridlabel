#include "Graph.h"
#include<fstream>
#include<algorithm>
#include<sstream>
#include<iostream>
#include<queue>

using namespace std;

Graph::Graph(string graphFile, int weighted, int directed){
	this->graphFile = graphFile;
	this->weighted = weighted;
	this->directed = directed;
	this->rv = 0;
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
	olabel.resize(V);
	olabeld.resize(V);
	hneighbors.resize(V);
	hreach.resize(V);	
	com.resize(V, 1);
	glabel.resize(V);
	glabeld.resize(V);
	nlabel.resize(V);
	nlabeld.resize(V);
	plabel.resize(V);
	pset.resize(V);

	ifstream ifs(filename);
	ifs >> V;
	if(!ifs.good()) return 0;
	int v;
	int ln;
	while(ifs >> v >> ln){
		vector<int>& label_v = label[v];
		vector<int>& labeld_v = labeld[v];
		
		vector<int>& glabel_v = glabel[v];
		vector<int>& glabeld_v = glabeld[v];

		vector<int>& nlabel_v = nlabel[v];
		vector<int>& nlabeld_v = nlabeld[v];

		set<int>& hneighbors_v = hneighbors[v];
		int u,d,s;
		int cnt = 0;
		for(int i = 0; i < ln; i++){
			ifs >> u >> d >> s;	
			label_v.push_back(u);
			labeld_v.push_back(d);
			olabel[v].push_back(u);
			olabeld[v].push_back(d);

			if(d == 1 || s == 0){
				nlabel_v.push_back(u);
				nlabeld_v.push_back(d);
			}else{
				glabel_v.push_back(u);
				glabeld_v.push_back(d);
			}

			if( d == 1)
				hneighbors_v.insert(u);
		}
		if( v == V-1){
			ifs >> olabel_amount;
			break;
		}
	}
	int sumg = 0;
	for(v = 0; v < V; v++){
		sumg += glabel[v].size();
	}

	cout << sumg << "/" << olabel_amount << endl;

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

//literally decompose PLL labels into three classes: direct links, shortcuts, guidences(which can be dropped)
int Graph::decomposeLabel(){
	cout << "Start decomposing!" << endl;
	int reduce_count = 0;

	com.resize(V);
	glabel.resize(V);
	glabeld.resize(V);
	nlabel.resize(V);
	nlabeld.resize(V);
	plabel.resize(V);
	pset.resize(V);

	for(int i = 0; i < V; i++){
		int v = dorder[i];//get the ith rank vertice
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
				if( it != hreach[v].end()){
					for(auto& e : hneighbors[v]){
						if(pllquery(e,label[v][j]) == (labeld[v][j] - 1)){
							reduce_count++;
							glabel[v].push_back(label[v][j]);
							glabeld[v].push_back(labeld[v][j]);
							break;
						}else{
							nlabel[v].push_back(label[v][j]);
							nlabeld[v].push_back(labeld[v][j]);
						}
					}
				}
				else{
					nlabel[v].push_back(label[v][j]);
					nlabeld[v].push_back(labeld[v][j]);
				}
			}else{
					nlabel[v].push_back(label[v][j]);
					nlabeld[v].push_back(labeld[v][j]);
			}
		}
		if(i%5000 == 0){
			cout << "progress: " << i << "/" << V << "," << reduce_count << endl;
			cout << label[v].size() << " : " << hreach[v].size() << endl;
		}
		com[v] = 1;
	}
	cout << reduce_count << "/" << olabel_amount << endl;
	return 1;
}

int Graph::compressLabel(const char* filename){
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

int Graph::lightenLabel(int v, int num){
	if(num <= 0) return 0;

	vector<int> neighbor_v;
	vector<int> neighbord_v;

	for(int i = 0; i < nlabel[v].size(); i++){
		if(nlabeld[v][i] == 1){
			neighbor_v.push_back(nlabel[v][i]);
			neighbord_v.push_back(1);
		}	
	}

	vector<int> &glabel_v = glabel[v];
	vector<int> &glabeld_v = glabeld[v];
	vector<pair<int, int> > gld;


	for(int i = 0; i < glabel_v.size(); i++)	gld.push_back(make_pair(glabeld_v[i], glabel_v[i]));	
	sort(gld.begin(), gld.end());
	
	vector<pair<int, int> > &plabel_v = plabel[v];

	set<int>::iterator it;
	for(int i = 0; i < num; i++){
		if(gld.size() == 0) break;

		int hv = gld[gld.size()-1].second;
		int hvd = gld[gld.size()-1].first;	
		int pv = -1;
		int pvd = -1;
		for(int j = gld.size() - 2; j >= 0; j--){
			int hu = gld[j].second;
			int hud = gld[j].first;
			if(hud == (hvd - 1)){
				if(pllquery(hu, hv) == 1 && drank[hu] > drank[hv]){
					pv = hu;
					pvd = hud;
					break;
				}
			}	
			/*
			if(j < 0) break;
			int hu = gld[j].second;
			int hud = gld[j].first;
			it = hreach[hu].find(hv);	
			if( it != hreach[hu].end() ){
				if(pllquery(hu, hv) == (hud - hvd)){
					pv = hu;
					pvd = hud;
					break;
				}
			}
			if( j == 0 && pv == -1){
				for(int k = 0; k < neighbor_v.size(); k++){
					hu = neighbor_v[k];
					hud = 1;
					it = hreach[hu].find(hv);
					if(it != hreach[hu].end() ){
						if(pllquery(hu, hv) == (hud - hvd)){
							pv = hu;
							pvd = hud;
							break;
						}
					}
				}
			}
			*/
		}	
		if( pv == -1){
			for(int k = 0; k < neighbor_v.size(); k++){
				int hu = neighbor_v[k];
				int hud = 1;
				if(hud == (hvd - 1)){
					if(pllquery(hu, hv) == 1 && drank[hu] > drank[hv]){
						pv = hu;
						pvd = hud;
						break;
					}
				}	
			}
		}
		if(pv != -1){
			it = pset[v].find(pv);
			if( it == pset[v].end() ){
				plabel_v.push_back(make_pair(-pvd, pv));	
				pset[v].insert(pv);
				sort(plabel_v.rbegin(), plabel_v.rend());
				com[v] = 0;
			}
			glabel_v.pop_back();
			glabeld_v.pop_back();
			gld.pop_back();
			rv++;
		}
	}
	
	//re organized the labels for PLL search
	vector<pair<int, int> > slabelv;
	vector<int> &label_v = label[v];
	vector<int> &labeld_v = labeld[v];

	for(int i = 0; i < nlabel[v].size(); i++) slabelv.push_back(make_pair(nlabel[v][i], nlabeld[v][i]));
	for(int i = 0; i < glabel[v].size(); i++) slabelv.push_back(make_pair(glabel[v][i], glabeld[v][i]));

	sort(slabelv.begin(), slabelv.end());

	label_v.resize(slabelv.size());
	labeld_v.resize(slabelv.size());
	for(int i = 0; i < slabelv.size(); i++){
		label_v[i] = slabelv[i].first;
		labeld_v[i] = slabelv[i].second;
	}
	slabelv.clear();

	return 1;
}

int Graph::pllquery(int u, int v){
	if(u==v)
		return 0;	

	int mind = INT_MAX;
	int i = 0;
	int j = 0;
	vector<int>& label_v = olabel[v];
	vector<int>& labeld_v = olabeld[v];
	vector<int>& label_u = olabel[u];
	vector<int>& labeld_u = olabeld[u];

	for( ; i < label_v.size() && j < label_u.size(); ){
		
		int lv = label_v[i];
		int lvd = labeld_v[i];
		int lu = label_u[j];
		int lud = labeld_u[j];
	
		if(lv == lu){
			if(lvd + lud < mind)
				mind = lvd + lud;
			i++; j++;
		}else{
			if(lv > lu) j++;
			if(lu > lv) i++;
		}


	}

	return mind;
}

int Graph::query(int u, int v){

	if(u==v)
		return 0;	

	int mind = INT_MAX;

	int i = 0;
	int j = 0;
	vector<int>& label_v = label[v];
	vector<int>& labeld_v = labeld[v];
	vector<int>& label_u = label[u];
	vector<int>& labeld_u = labeld[u];

	for( ; i < label_v.size() && j < label_u.size(); ){
		
		int lv = label_v[i];
		int lvd = labeld_v[i];
		int lu = label_u[j];
		int lud = labeld_u[j];
	
		if(lv == lu){
			if(lvd + lud < mind)
				mind = lvd + lud;
			i++; j++;
		}else{
			if(lv > lu) j++;
			if(lu > lv) i++;
		}


	}
	
	if(com[u] == 1 && com[v] == 1){//merge join search just as pll

		/*	
		vector<int> vvisited(V, 0);
		vector<int> uvisited(V, 0);
		vector<int> dist(V, 0);
		priority_queue<pair<int, int> > vheap, uheap;
		int vmin = INT_MAX;
		int umin = INT_MAX;
		
		for(int i = 0; i < nlabel[v].size(); i++){
			vvisited[nlabel[v][i]] = 1;
			dist[nlabel[v][i]] += nlabeld[v][i];
		}
		for(int i = 0; i < glabel[v].size(); i++){
			vvisited[glabel[v][i]] = 1;
			dist[glabel[v][i]] += glabeld[v][i];
		}
		for(int i = 0; i < nlabel[u].size(); i++){
			uvisited[nlabel[u][i]] = 1;
			dist[nlabel[u][i]] += nlabeld[u][i];
			if(vvisited[nlabel[u][i]] == 1 && uvisited[nlabel[u][i]] == 1)
				if(dist[nlabel[u][i]] < mind) mind = dist[nlabel[u][i]];
		}
		for(int i = 0; i < glabel[u].size(); i++){
			uvisited[glabel[u][i]] = 1;
			dist[glabel[u][i]] += glabeld[u][i];
			if(vvisited[glabel[u][i]] == 1 && uvisited[glabel[u][i]] == 1)
				if(dist[glabel[u][i]] < mind) mind = dist[glabel[u][i]];	
		}*/
		
		return mind;
	}	
	if( plabel[v].size() != 0 && plabel[u].size() == 0 && mind <= -plabel[v][0].first)
		return mind;
	
	if( plabel[v].size() == 0 && plabel[u].size() != 0 && mind <= -plabel[u][0].first)
	       	return mind;
	if( plabel[v].size() != 0 && plabel[u].size() != 0 && mind <= -plabel[v][0].first && mind <= -plabel[u][0].first )
	       	return mind;

	vector<int> vvisited(V, 0);
	vector<int> uvisited(V, 0);
	vector<int> dist(V, 0);
	vector<int> vdist(V,INT_MAX);
	vector<int> udist(V,INT_MAX);
	priority_queue<pair<int, int> > vheap, uheap;
	int vmin = INT_MAX;
	int umin = INT_MAX;
	
	for(i = 0; i < label_v.size(); i++){
	       	vvisited[label_v[i]] = 1;
		dist[label_v[i]] += labeld_v[i];
	}
	for(i = 0; i < label_u.size(); i++){
	       	uvisited[label_u[i]] = 1;
		dist[label_u[i]] += labeld_u[i];	
	}
/*	
	for(int i = 0; i < nlabel[v].size(); i++){
		vvisited[nlabel[v][i]] = 1;
		dist[nlabel[v][i]] += nlabeld[v][i];
	}
	for(int i = 0; i < glabel[v].size(); i++){
		vvisited[glabel[v][i]] = 1;
		dist[glabel[v][i]] += glabeld[v][i];
	}
	for(int i = 0; i < nlabel[u].size(); i++){
		uvisited[nlabel[u][i]] = 1;
		dist[nlabel[u][i]] += nlabeld[u][i];
		if(vvisited[nlabel[u][i]] == 1 && uvisited[nlabel[u][i]] == 1)
			if(dist[nlabel[u][i]] < mind) mind = dist[nlabel[u][i]];
	}
	for(int i = 0; i < glabel[u].size(); i++){
		uvisited[glabel[u][i]] = 1;
		dist[glabel[u][i]] += glabeld[u][i];
		if(vvisited[glabel[u][i]] == 1 && uvisited[glabel[u][i]] == 1)
			if(dist[glabel[u][i]] < mind) mind = dist[glabel[u][i]];	
	}
*/
	for(i = 0; i < plabel[v].size(); i++){
		pair<int, int> e = plabel[v][i];
		if(vvisited[e.second] == 1)
			dist[e.second] += e.first;//this dist should be subtracted
		vvisited[e.second] = 0;
		vdist[e.second] = -e.first;
		//dist[e.second] = -e.first;
		vheap.push(make_pair(e.first, e.second));		
	}
	for(i = 0; i < plabel[u].size(); i++){
		pair<int, int> e = plabel[u][i];
		if(uvisited[e.second] == 1)
			dist[e.second] += e.first;
		uvisited[e.second] = 0;
		udist[e.second] = -e.first;
		//dist[e.second] = -e.first;
		uheap.push(make_pair(e.first, e.second));		
	}
	while(vheap.size() != 0 || uheap.size() != 0){	
		if( vheap.size() != 0 && mind <= -vheap.top().first && uheap.size() == 0)//lower bound can be stronger by replacing with vheap + label_u and uheap + label_v
			return mind;
		if( uheap.size()!= 0 && mind <= -uheap.top().first && vheap.size() == 0)
			return mind;

		if( vheap.size() != 0 && uheap.size() != 0 && mind <= -uheap.top().first && mind <= -vheap.top().first)
			return mind;
		
		if(vheap.size() != 0){
			int vtop = vheap.top().second;
			int vtopd = -vheap.top().first;	
			vheap.pop();
			
			if(vvisited[vtop] == 0){//or else it has been fixed in nlabel
				vvisited[vtop] = 1;//has problem here
				dist[vtop] += vtopd;	
			}else continue;
			if(uvisited[vtop] == 1){
				if(dist[vtop] < mind) mind = dist[vtop];
				return mind;
			}
			for(i = 0; i < label[vtop].size(); i++){
				int vn = label[vtop][i];
				int vnd = labeld[vtop][i] + vtopd;
				if(vvisited[vn] == 1) continue;//vn has been fixed in nlabel
				if(vnd < vdist[vn]){
					vdist[vn] = vnd;
					vheap.push(make_pair(-vnd, vn));		
				}
			}		
		}	
		
		if(uheap.size() != 0){
			int utop = uheap.top().second;
			int utopd = -uheap.top().first;	
			uheap.pop();
			if(uvisited[utop] == 0){//or else it has been fixed in nlabel
				uvisited[utop] = 1;
				dist[utop] += utopd;	
			}else continue;
			if(vvisited[utop] == 1){
				if(dist[utop] < mind) mind = dist[utop];
				return mind;
			}
			for(i = 0; i < label[utop].size(); i++){
				int un = label[utop][i];
				int und = labeld[utop][i] + utopd;
				if(uvisited[un] == 1) continue;//vn has been fixed in nlabel
				if(und < udist[un]){
					udist[un] = und;
					uheap.push(make_pair(-und, un));		
				}
			}		
		}	
	}
	return mind;
}



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
