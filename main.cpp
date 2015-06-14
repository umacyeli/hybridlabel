#include<iostream>
#include<string>
#include "Graph.h"
#include <sys/time.h>

using namespace std;

int main(){
	//string graphFile = "/Users/jagsly/Downloads/graph_data/slashdot-threads/out.slashdot-threads";
	//string graphFile = "/Users/jagsly/Downloads/graph_data/facebook-wosn-links/out.facebook-wosn-links";
	//string graphFile = "/Users/jagsly/Downloads/graph_data/soc-Epinions1/out.soc-Epinions1";//u,v
	//string graphFile = "/Users/jagsly/Downloads/graph_data/youtube-u-growth/out.youtube-u-growth";//u,v,w,s
	//string graphFile = "/Users/jagsly/Downloads/graph_data/flickr-links/out.flickr-links";//u,v
	//string graphFile = "/Users/jagsly/Downloads/graph_data/ego-facebook/out.ego-facebook";
	//string graphFile = "/Users/jagsly/Downloads/graph_data/WikiTalk.txt";
	//string graphFile = "/Users/jagsly/Downloads/graph_data/loc-gowalla_edges/out.loc-gowalla_edges";
	string graphFile = "/Users/jagsly/Downloads/graph_data/loc-brightkite_edges/out.loc-brightkite_edges";//u,v
	//string graphFile = "/Users/jagsly/Downloads/graph_data/petster-friendships-dog/out.petster-friendships-dog-uniq";
	//string graphFile = "/Users/jagsly/Downloads/graph_data/petster-friendships-cat/out.petster-friendships-cat-uniq";
	Graph g(graphFile, 0, 0);	
	char* label_file = "/Users/jagsly/Downloads/graph_data/loc-brightkite_edges/out.loc-brightkite_edges.label.gpll";
	char* degree_file = "/Users/jagsly/Downloads/graph_data/loc-brightkite_edges/out.loc-brightkite_edges.label.gpll.drank";
	//char* label_file = "/Users/jagsly/Downloads/graph_data/soc-Epinions1/out.soc-Epinions1.label.gpll";
	//char* degree_file = "/Users/jagsly/Downloads/graph_data/soc-Epinions1/out.soc-Epinions1.label.gpll.drank";
	//char* new_label_file = "/Users/jagsly/Downloads/graph_data/loc-brightkite_edges/out.loc-brightkite_edges.newlabel";
	//char* label_file = "/Users/jagsly/Downloads/graph_data/slashdot-threads/out.slashdot-threads.label.pll";
	//char* degree_file = "/Users/jagsly/Downloads/graph_data/slashdot-threads/out.slashdot-threads.label.pll.drank";
	//char* new_label_file = "/Users/jagsly/Downloads/graph_data/slashdot-threads/out.slashdot-threads.newlabel";
	//char* label_file = "/Users/jagsly/Downloads/graph_data/example_graph/hub_edges.label.pll";
	//char* degree_file = "/Users/jagsly/Downloads/graph_data/example_graph/hub_edges.label.pll.drank";
	//char* new_label_file = "/Users/jagsly/Downloads/graph_data/example_graph/hub_edges.newlabel";
	
	srand (time(NULL));
	double useconds;
	
	g.loadRankOrder(degree_file);
	g.loadLabel(label_file);
	//g.decomposeLabel();
	for(int i = 0; i < 100; i ++){
		int times = 10000;	
		vector<int> v(times);
		vector<int> u(times);
		for(int j = 0; j < times; j++){
			v[j] = rand() % g.V;
			u[j] = rand() % g.V;
		}

		struct timeval start;
		gettimeofday(&start, NULL);
		int icount = 0;
		for(int j = 0; j < times; j++){
			int qr =  g.query(v[j], u[j]);
			
			int pr = g.pllquery(v[j],u[j]);
			
			if(pr < 999 && pr != qr){
				cout << pr << " ? " << qr << endl;
			       	icount ++;
			}
		}
		struct timeval stop;
		gettimeofday(&stop, NULL);
		useconds = (stop.tv_sec - start.tv_sec)*1000000 + stop.tv_usec - start.tv_usec;
		//cout << "query time: " << useconds/times << endl;
		//cout << "rg: " << g.rv << endl;
		cout << useconds/times << "," << g.rv << "," << icount << endl;
		for(int j = 0; j < g.V; j++)
			g.lightenLabel(j,1);
	}

}
