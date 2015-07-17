//
//  HLabel.h
//  hybridLabeling
//
//  Created by Jagsly on 7/12/15.
//  Copyright (c) 2015 Jagsly. All rights reserved.
//

#ifndef hybridLabeling_HLabel_h
#define hybridLabeling_HLabel_h

#include<vector>
#include<string>
#include<queue>
#include <boost/heap/d_ary_heap.hpp>

using namespace boost::heap;


using namespace std;

typedef struct label{//compact label
    int hid;//hub id
    int d;//hub distance
} label;

typedef struct flabel{//full label
    int hid;
    int d;
    int hp;//hub parent
} flabel;

typedef struct parent{
    int pid;//parent vertex id;
    int hop;//hop limit for this parent
    int pd;//parent distance
} parent;

struct compare_pairpair
{
    bool operator () (const pair<int, pair<int, int> >& p1, const pair<int, pair<int, int> >& p2) const{
        
        if( p1.first != p2.first )
            return p1.first < p2.first;
        
        if( p1.second.first != p2.second.first )
            return p1.second.first > p2.second.first;
        
        return p1.second.second > p2.second.second;
        
        
    }
};


typedef boost::heap::d_ary_heap<pair<int, pair<int, int> >, boost::heap::arity<2>, boost::heap::mutable_<true>, boost::heap::compare<compare_pairpair> >  BinaryHeap;

class HLabel{
    
public:
    //loading
    int load(string indexFile, vector<vector<flabel> > &fullindex);//load full label into memory
    int loadRank(string rankFile);//load rank order file

    //preprocessing
    int full2invert(vector<vector<flabel> > fullindex, vector<vector<label> > &invertedindex);//creating inverted index list
    int full2compact(vector<vector<flabel> > fullindex);//importing full label into index_
    int full2pll(vector<vector<flabel> > fullindex);//importing full label into pllindex_
  
    //removing label
    int removelabel(int v, vector<vector<flabel> > &fullindex, vector<vector<label> > invertedindex);//remove v from every label it exists(except itself)
    
    //query
    int pllquery(int s, int t);
    
    void initialization();// call it before using query
    int query(int s, int t);
    void clean();
    
    HLabel();
    ~HLabel();
    
    //stat
    double hittest(int v, vector<vector<label> > invertedindex, double samRate);//testing the hit ratio of v(monte carlo with sampling ratio samRate)
    int fanouttest(string foStat);// x-axis: from high rank to low rank vertices, y-axis: number of labels
    
    int parents_distribution(int n); //the rank removed to
    int parents_fanout_distribution(int n); //the rank removed to
    
    
    int V;//total number of Vertex
    int L;//total number of original Label
    int RL;//total number of Reduced Label

    vector<int> rank;//get rank of v by rank[v]
    vector<int> inv;//get vertex of i-th rank by inv[i]
    
private:
        vector<vector<label> > index_;//final labels used to answer the queries
    vector<vector<parent> > parents_;//final parents used to answer the queries
    vector<vector<label> > pllindex_;//pll index

    //essential structure for distance query
    vector<int> dlist;//a hash structure for distance merging
    queue<int> vque;// visit queue for each
    vector<int> visited;// 1 - visited by s, 2 - visisted by t, 3 - visited by both sides
    vector<BinaryHeap::handle_type> handles[2];
    vector<bool> reach[2];
    
    queue<int> pvque;//mark the parents
    vector<bool> pmark[2];//mark the parents
};


#endif
