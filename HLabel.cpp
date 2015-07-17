//
//  HLabel.cpp
//  hybridLabeling
//
//  Created by Jagsly on 7/12/15.
//  Copyright (c) 2015 Jagsly. All rights reserved.
//

#include "HLabel.h"
#include<fstream>
#include<iostream>
#include<set>

using namespace std;


HLabel::HLabel(){
    V = -1;
    L = -1;
    RL = 0;
}

HLabel::~HLabel(){
    index_.clear();
    parents_.clear();
    pllindex_.clear();
    rank.clear();
    inv.clear();
}

bool parentcompare(parent lhs, parent rhs) { return lhs.pd < rhs.pd; }

int HLabel::load( string indexFile, vector<vector<flabel> > &fullindex){
    
    ifstream ifs(indexFile);
    if( ifs.bad() == true)
        return 0;
    
    ifs >> V;
    
    //Initialization
    fullindex.clear();
    fullindex.resize(V);
    parents_.clear();
    parents_.resize(V);
    
    for(int i = 0; i < V; i++){
        int v, vlabel_num;
        ifs >> v >> vlabel_num;
        for(int j = 0; j < vlabel_num; j++){
            int uid, ud, up;
            ifs >> uid >> ud >> up;
            flabel ul = {uid, ud, up};
            fullindex[i].push_back(ul);
        }
    }
    
    ifs >> L;

    if (ifs.bad()) return 0;
    ifs.close();
    return 1;
}

int HLabel::loadRank(string rankFile){
    
    ifstream ifs(rankFile);
    if( ifs.bad() == true)
        return 0;
    
    ifs >> V;
    
    //Initialization
    rank.clear();
    rank.resize(V);
    inv.clear();
    inv.resize(V);
    
    int v, r;
    while( ifs >> r >> v){
        rank[v] = r;
        inv[r] = v;
    }
    
    if (ifs.bad()) return 0;
    ifs.close();
    
    return 1;
}

int HLabel::full2invert(vector<vector<flabel> > fullindex, vector<vector<label> > &invertedindex){
    
    invertedindex.clear();
    invertedindex.resize(V);
    
    for(int i = 0; i < V; i++){
        for(int j = 0; j < fullindex[i].size(); j++){
            flabel ufl = fullindex[i][j];
            label ul = {i, ufl.d};
            invertedindex[ufl.hid].push_back(ul);
        }
    }
    
    return 1;
}

int HLabel::full2pll(vector<vector<flabel> > fullindex){
    
    pllindex_.clear();
    pllindex_.resize(V);
    
    for(int i = 0; i < V; i++){
        for(int j = 0; j < fullindex[i].size(); j++){
            flabel ufl = fullindex[i][j];
            label ul = {ufl.hid, ufl.d};
            pllindex_[i].push_back(ul);
        }
    }
    return 1;
}

int HLabel::full2compact(vector<vector<flabel> > fullindex){
    
    index_.clear();
    index_.resize(V);
    
    for(int i = 0; i < V; i++){
        for(int j = 0; j < fullindex[i].size(); j++){
            flabel ufl = fullindex[i][j];
            label ul = {ufl.hid, ufl.d};
            index_[i].push_back(ul);
        }
    }
    return 1;
}

//the vertex v can be removed from the labels of vertex u only under certain condition:
//(1) v != u
//(2) the parent of v in the label of u is not u(pv != u)
//after removing v, we need to merge p(v) and the parent list parents_[u] based on:
//(1) if p(v) exists in parents_[u], we do nothing or remove v in parents_[u] (RL++)
//(2) if v exists in parents_[u], we replace (v, v_hop, vd) by (p(v), v_hop + 1, p(v)d) or we do nothing when p(v) has already been in the parent list(RL++)
//(3) if nothing match c1 and c2 in parents_[u] add (p(v), v_hop=1, p(v)d)
int HLabel::removelabel(int v, vector<vector<flabel> > &fullindex, vector<vector<label> > invertedindex){
    
    for(int i = 0; i < invertedindex[v].size(); i++){//traversing all vertices(u) containing v as a label
        int u = invertedindex[v][i].hid;
        int dvu = invertedindex[v][i].d;
        
        if( v == u) continue;// c1
        
        for(int j = 0; j < fullindex[u].size(); j++){
            int w = fullindex[u][j].hid;
            if( w == v){//found the v in the label of u
                int pv = fullindex[u][j].hp;//the parent of v in the label of u
                if( pv != u){//c2, starting to merge parents
                    
                    for(int k = 0; k < parents_[u].size(); k++){
                        int tp = parents_[u][k].pid;
                        if(tp == pv){//c1, we do nothing
                           
                            if(tp != v )
                               RL++;
                            
                            pv = -1;
                            
                            for(int q = 0; q < parents_[u].size(); q++){
                                int tmpv = parents_[u][q].pid;
                                if ( tmpv == v ){//if v is in the parent list and p(v) has been added, then remove v in the parent list
                                    if( parents_[u][k].hop < ( parents_[u][q].hop + 1 ) ){
                                        parents_[u][k].hop = parents_[u][q].hop + 1;
                                    }
                                    parents_[u].erase(parents_[u].begin() + q);
                                    break;
                                }
                            }

                        }
                        
                        if(tp == v){//c2, we replacing old parents
                            
                            bool oldexisted_flag = false;
                            //first we check whether pv has already existed
                            for(int q = 0; q < parents_[u].size(); q++){
                                int tmpv = parents_[u][q].pid;
                                if( tmpv == pv ){//the old parent has already existed
                                    oldexisted_flag = true;
                                    //update the pv in parent list
                                    if(parents_[u][q].hop < (parents_[u][k].hop+1) )
                                        parents_[u][q].hop = parents_[u][k].hop + 1;
                                    //and erase v from parent list
                                    parents_[u].erase(parents_[u].begin() + k );
                                }
                            }
                            
                            if(oldexisted_flag == false ){
                                parents_[u][k].pid = pv;
                                parents_[u][k].hop++;
                            
                                for(int q = 0; q < fullindex[u].size(); q++){
                                    int tmpv = fullindex[u][q].hid;
                                    if ( tmpv == pv ){
                                        parents_[u][k].pd = fullindex[u][q].d;
                                        break;
                                    }
                                }
                            }
                            
                            RL++;
                            pv = -1;
                            break;
                        }
                    }
                    
                    if(pv != -1){
                        int npd = -1;
                        for(int q = 0; q < fullindex[u].size(); q++){
                            int tmpv = fullindex[u][q].hid;
                            if ( tmpv == pv ){
                                npd = fullindex[u][q].d;
                                break;
                            }
                        }
                        parent np = {pv, 1, npd};
                        parents_[u].push_back(np);
                    }
                   
                    //removing the label
                    fullindex[u].erase(fullindex[u].begin() + j);
                    
                }
            }
        }
        
        sort(parents_[u].begin(), parents_[u].end(), parentcompare);
        
    }
    return 1;
}

int HLabel::pllquery(int s, int t){
    if( s == t ) return 0;
    
    int distance = INT32_MAX;
    int i = 0, j = 0;
    vector<label> index_s = pllindex_[s], index_t = pllindex_[t];
    
    for( ; i < index_s.size() && j < index_t.size() ; ){
        int sv = index_s[i].hid;
        int sd = index_s[i].d;
        int tv = index_t[j].hid;
        int td = index_t[j].d;
        
        if( sv == tv ){
            int tmp_distance = sd + td;
            
            if( tmp_distance < distance )
                distance = tmp_distance;
            
            i++; j++;
            
        }else{
            if (sv > tv) j++;
            if (sv < tv) i++;
        }
        
    }
    
    return distance;
}

//call it before using query(int s, int)
void HLabel::initialization(){
    
    //need to re initialize
    dlist.resize(V);
    visited.resize(V);
    reach[0].resize(V);
    reach[1].resize(V);
    
    //no need to re initialize
    handles[0].resize(V);
    handles[1].resize(V);
    
    pmark[0].resize(V);
    pmark[1].resize(V);
    
    return;
}

//3 phases for distance query:
//(1) merge current label using dlist to get the current min distance
//(2) check the first parents in both lists(sorted by distances), if mind <= both of first parents, then return;
//(3) push the parents whose distance is smaller than the mind for both lists to both heap
//(4) round-and-robin fashion deal with the heap whose current top is smaller than current mind

int HLabel::query(int s, int t){

    if( s == t ) return 0;
    
    int distance = INT32_MAX;
    int i = 0, j = 0;
    vector<label> index_s = index_[s], index_t = index_[t];
    
    for( ; i < index_s.size() && j < index_t.size() ; ){
        int sv = index_s[i].hid;
        int sd = index_s[i].d;
        int tv = index_t[j].hid;
        int td = index_t[j].d;
        
        if(reach[0][sv] == false){
            vque.push(sv);
            dlist[sv] += sd;
            visited[sv] += 1;
            reach[0][sv] = true;
        }
        
        if(reach[1][tv] == false){
            vque.push(tv);
            dlist[tv] += td;
            visited[tv] += 2;
            reach[1][tv] = true;
        }
        
        if( sv == tv ){
            int tmp_distance = sd + td;
            if( tmp_distance < distance )
                distance = tmp_distance;
            i++; j++;
            
        }else{
            if (sv > tv) j++;
            if (sv < tv) i++;
        }
        
    }
 //   if( s == 6121 && t == 4849){
   //     int kkkk  = 0;
    //}
    vector<parent> parent_s = parents_[s], parent_t = parents_[t];
    
    if( parent_s.size()!=0 && distance <= parent_s[0].pd && parent_t.size() !=0 && distance <= parent_t[0].pd) { clean(); return distance; }
    
    
    if( parent_s.size() == 0  && parent_t.size() !=0 && distance <= parent_t[0].pd) { clean(); return distance; }
    
    if( parent_s.size() != 0  && distance <= parent_s[0].pd && parent_t.size() == 0) { clean(); return distance; }
    
    
    BinaryHeap heap[2];//heap[0] for sHeap, heap[1] for tHeap
    
    /*
    for(int i = 0; i < parent_s.size(); i++)
        if(parent_s[i].pd < distance){
       //     reach[0][parent_s[i].pid] = true;
         //   vque.push(parent_s[i].pid);
            handles[0][parent_s[i].pid] = heap[0].push(make_pair(-parent_s[i].pd, make_pair(parent_s[i].pid, parent_s[i].hop)));
            
            pmark[0][parent_s[i].pid] = true;
            pvque.push(parent_s[i].pid);
        }
    
    for(int i = 0; i < parent_t.size(); i++)
        if(parent_t[i].pd < distance){
           // reach[1][parent_t[i].pid] = true;
           // vque.push(parent_t[i].pid);
            handles[1][parent_t[i].pid] = heap[1].push(make_pair(-parent_t[i].pd, make_pair(parent_t[i].pid, parent_t[i].hop)));
            
            pmark[1][parent_t[i].pid] = true;
            pvque.push(parent_t[i].pid);
        }
    
    int roller = 1;//0 for heap[0] (s heap) 1 for heap[1] (t heap)
    
    
    while( (heap[0].size() != 0 && -heap[0].top().first < distance )|| ( heap[1].size() != 0 && -heap[1].top().first < distance )){// TODO: max test
        if(roller == 1)
            roller = 0;
        else//roller = 0
            roller = 1;
        
        BinaryHeap &tmpheap = heap[roller];
        vector<bool> &tmpreach = reach[roller];
        vector<BinaryHeap::handle_type> &tmphandles = handles[roller];

        if(tmpheap.size() == 0 || distance <= -tmpheap.top().first) continue;//forget about the invalid heap
        

        
        int cv = tmpheap.top().second.first;
        int chop = tmpheap.top().second.second;
        int cd = -tmpheap.top().first;
        
        tmpheap.pop();
        
        if( (visited[cv] == (roller + 1) || visited[cv] == 3 ) && pmark[roller][cv] == false ) continue;//it has been fixed
    
        if( pmark[roller][cv] == false){
            visited[cv] += (roller + 1);// 0+1 - visited by s, 1+1 - visited by t
            dlist[cv] += cd;
            tmpreach[cv] = true;
            vque.push(cv);
        }
        
        //if (s == 4250 && t == 3893) {
    //    if ( s == 2356 && t == 1468 ){
        //if( s == 6121 && t == 4849){
          //  cout << "heap[" << roller << "] is popping vertex " << cv << " d=" << cd << endl;
        //}
        
        if(visited[cv] == 3)//scanned by both side
            if( dlist[cv] < distance )
                distance = dlist[cv];
        
        if(chop == 0) continue; //we do not extend the search from cv
        
    
        //TODO: push update...
        vector<label> index_top = index_[cv];
        for( int i = 0; i < index_top.size(); i++ ){
            int nv = index_top[i].hid;//new vertex
            if( nv == cv ) continue;
            int nd = index_top[i].d;// new distance
            int thed = nd + cd;
            
     
            
            if( (visited[nv] == 1 || visited[nv] == 3) && roller == 0)//nv has been fixied in heap[0]
                continue;
            if( (visited[nv] == 2 || visited[nv] == 3) && roller == 1)//nv has been fixied in heap[1]
                continue;
            
            if( tmpreach[nv] == false){//directly push into heap

                tmphandles[nv] = tmpheap.push(make_pair(-thed, make_pair(nv, chop-1)));//push new distance
                
    
                tmpreach[nv] = true;
                
                
                vque.push(nv);
                
            }else{//update the heap
                if(visited[nv] != (roller+1) && visited[nv] != 3){
                    if( -(*tmphandles[nv]).first > thed){//need to update
                  //      if(nv == 1585)
                    //        cout << -(*tmphandles[nv]).first << "," << thed << endl;
                 //       cout << "s=" << s << ",t=" << t << " nv=" << nv << " ind " << -(*tmphandles[nv]).first << "," << "thed " << thed << " heapsize=" << tmpheap.size()  << endl;
                        
                        tmpheap.increase(tmphandles[nv], make_pair(-thed, make_pair(nv, chop-1)));//chop-1 is correct?
                    }
                }
                
            }
            
            
        }
        
    }
    */
    
    heap[0].clear();
    heap[1].clear();
    
    clean();
    
    return distance;
}

int HLabel::fanouttest(string foStat){
    ofstream ofs(foStat);
    
    for(int i = 0; i < V; i++){
        int v = inv[i];
        int vs = pllindex_[v].size();
        ofs << v << " " << vs << endl;
    }
    ofs.close();
    return 1;
}

void HLabel::clean(){

    while(!vque.empty()){
        dlist[vque.front()] = 0;
        visited[vque.front()] = 0;
        reach[0][vque.front()] = false;
        reach[1][vque.front()] = false;
        vque.pop();
    }
    
    while(!pvque.empty()){
        pmark[0][pvque.front()] = false;
        pmark[1][pvque.front()] = false;
        pvque.pop();
    }
    
    for(int i = 0; i < dlist.size(); i++){
        if(visited[i] !=0 || reach[0][i] != false || reach[1][i] != false || dlist[i] != 0)
            cout << "noooooooo!" << endl;
    }
}

double HLabel::hittest(int v, vector<vector<label> > invertedindex, double samRate){
    
    int vn = invertedindex[v].size();//
    int pn = vn * vn;//the amount of pairs
    
    int sn = samRate * pn;//sampling ratio
    
    bool hit = false;
    int hn = 0;//hit number
    
    srand( 256 );//static sampling
    //srand( time( 0 ) );
    
    while(sn != 0 ){
        sn--;
        
        int s = rand()%vn;
        int t = rand()%vn;
        
        if( s == t && s == v ){ hn++; continue;}
        if( s == t && s != v ) continue;
        
        int distance = INT32_MAX;
        int hv = -1;//hit vertices
        int i = 0, j = 0;
        vector<label> index_s = pllindex_[s], index_t = pllindex_[t];
        
        for( ; i < index_s.size() && j < index_t.size() ; ){
            int sv = index_s[i].hid;
            int sd = index_s[i].d;
            int tv = index_t[j].hid;
            int td = index_t[j].d;
            
            if( sv == tv ){
                int tmp_distance = sd + td;
                if( tmp_distance <= distance ){
                    if(distance == tmp_distance){//we assign the sv(tv) != v to hv
                        if( hv == v )//change it
                            hv = sv;
                        //else we do nothing
                    }else{//distance < tmp_distance
                        hv = sv;
                        distance = tmp_distance;
                    }
                }
                i++; j++;
                
            }else{
                if (sv > tv) j++;
                if (sv < tv) i++;
            }
            
        }
        
        if( hv == v)
            hn++;
    
    }
    sn = samRate * pn;
    cout << "hn=" << hn << ", sn=" << sn << endl;
    double hitrate = (double)hn/(double)sn;
    return hitrate;
}

int HLabel::parents_distribution(int n){
    set<int> tnr;
    for(int i = 0; i < n; i ++ ){
        int topv = inv[i];
        tnr.insert(topv);
    }
    for(int i = 0; i < V; i++){
        for(int j = 0; j < parents_[i].size(); j++){
            int pid = parents_[i][j].pid;
            tnr.insert(pid);
        }
    }
    return tnr.size();
}

int HLabel::parents_fanout_distribution(int n){
    set<int> tnr;
    set<int> non_tnr;
    int top_fanout = 0;
    int nontop_fanout = 0;
    for(int i = 0; i < n; i ++ ){
        int topv = inv[i];
        top_fanout += index_[topv].size();
        tnr.insert(topv);
    }
    
    set<int>::iterator finder;
    for(int i = 0; i < V; i++){
        for(int j = 0; j < parents_[i].size(); j++){
            int pid = parents_[i][j].pid;
            finder = tnr.find(pid);
            if(finder == tnr.end()){
                finder = non_tnr.find(pid);
                if(finder == non_tnr.end()){
                    non_tnr.insert(pid);
                    nontop_fanout += index_[parents_[i][j].pid].size();
                }
            }
        }
    }
    
    cout << top_fanout/tnr.size() << " vs " << nontop_fanout/non_tnr.size() << endl;
    
    
    return 1;
}