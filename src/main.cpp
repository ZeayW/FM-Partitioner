#include "net.h"
#include "cell.h"
#include <ctime>
#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <ctype.h>
#include <set>
#include <map>
#include <cmath>
#include <list>
#include <vector>


using namespace std;

vector <int> cellstack;
vector <Net*> nets, *bestnets;
vector <Cell*> cells, *bestcells;

int k = 0, bestk;
bool *bestset;
vector <int> bestA, bestB;
map <string, int> cell2id, net2id; 
int ccnt = 0, ncnt = 0;

double error; 
int totalCellsize = 0, aCellsize = 0, bCellsize = 0, cutsz = 0;
int bestacnnt = 0, bestbcnnt = 0, bestaCellsize = 0, bestbCellsize = 0;
int aCellCnt = 0, bCellCnt = 0, aGain = 0, bestg = 0;
int Pmax = 0;

int cutSize = 0;
ifstream ifsCell, ifsNet;
ofstream of;

double tstart, tend;
map <int, Node*> bucketlist[2];

void parseCells(istream & in){
    string str;
    int size;
    while (in >> str >> size){
        cell2id[str] = ccnt;
        if (aCellsize <= bCellsize){
            Cell *c = new Cell(str, size, 0, ccnt);
            cells.push_back(c);
            aCellsize += size;
            aCellCnt++;
        }
        else {
            Cell *c = new Cell(str, size, 1, ccnt);
            cells.push_back(c);
            bCellsize += size;
            bCellCnt++;
        }
        totalCellsize += size; 
        ccnt++;
    }
}

void calAB(){
    for (int i = 0; i < ncnt; i++){
        vector <int> & vl = nets[i]->cellList;
        nets[i]->B = 0;
        nets[i]->A = 0;
        for (int j = 0; j < vl.size(); j++){
            Cell * cell = cells[vl[j]];
            if (cell->set) nets[i]->B++;
            else nets[i]->A++;
        }
    }
}

void parseNets(istream & in){
    string str, tmp;
    while (in >> tmp){ // NET
        in >> str;  // nxxx
        net2id[str] = ncnt;
        in >> tmp;  // {
        Net *n = new Net(str);
        nets.push_back(n);
        while (in >> tmp && tmp[0] != '}'){
            vector <int> & l = cells[cell2id[tmp]]->netList;
            if (!l.size() || l[l.size()-1] != ncnt) {
                l.push_back(ncnt);
                cells[cell2id[tmp]]->pins++;
                nets[ncnt]->cellList.push_back(cell2id[tmp]);
                if (cells[cell2id[tmp]]->set) nets[ncnt]->B++;
                else nets[ncnt]->A++;
            }
        }
        ncnt++;
    }
}

void outputFile(ostream & out){
    out << "cut_size " << cutsz << endl;
    out << "A " << aCellCnt << endl;
    for (int i = 0; i < ccnt; i++)
        if (!cells[i]->set)
            out << cells[i]->name << endl;
    out << "B " << bCellCnt << endl;
    for (int i = 0; i < ccnt; i++)
        if (cells[i]->set)
            out << cells[i]->name << endl;
}


void parse(int argc, char ** argv){
    char* cellFile = argv[1];
    char* netFile = argv[2];
    char* outputFile = argv[3];
    ifsCell.open(cellFile, ios::in);
    if (!ifsCell.is_open()){
        cerr << "Cannot open the cells file: " << cellFile << endl;		
        exit(-1);
    }
    else{
        parseCells(ifsCell);
    }
    ifsCell.close();

    ifsNet.open(netFile, ios::in);
    if (!ifsNet.is_open()){
        cerr << "Cannot open the nets file: " << netFile << endl;	
        exit(-1);
    }
    else{
        parseNets(ifsNet);
    }
    ifsNet.close();

    of.open(outputFile, ios::out);
    if (!of.is_open())
        cerr << "Cannot open the output file: " << outputFile << endl;		
    
}

void remove(Cell * c){
    Node *p = c->to;
    p->prev->next = p->next;
    if (p->next != NULL) p->next->prev = p->prev;
}

// insert cell c with gain g to the front of bucket_list[g]
void insert_front(Cell * c){
    int gain = c->gain;
    bool set = c->set;
    Node *p = c->to;
    Node *pre = bucketlist[set][gain];
    Node* cur = bucketlist[set][gain]->next;
    //while(cur!=NULL && cells[cur->id]->size<c->size){
    //    cout<<p->id<<" "<<c->size<<" "<< cells[cur->id]->size<<endl;
    //    pre = pre->next;
    //    cur = cur->next;
    //}
    //cout<<endl;
    p->prev = pre;
    p->next = cur;
    pre->next = p;
    
    //p->prev = bucketlist[set][gain];
    //p->next = bucketlist[set][gain]->next;
    //bucketlist[set][gain]->next = p;
    if (p->next != NULL) p->next->prev = p;
}

void move(Cell * c){
    remove(c);
    insert_front(c);
}

void reverse(){
    int i = cellstack.size()-1;
    for (; i > k; i--)
        cells[cellstack[i]]->set = !cells[cellstack[i]]->set;
}

// build the bucket list
void buildbucketlist(){
    bucketlist[0].clear();
    bucketlist[1].clear();
    // init the bucket list( head node for all p)
    for (int i = -Pmax; i <= Pmax; i++) {
        if (bucketlist[0][i] == NULL) bucketlist[0][i] = new Node(-1);
        if (bucketlist[1][i] == NULL) bucketlist[1][i] = new Node(-1);
    }
    // insert all the cells to the bucket list
    for (int i = 0; i < ccnt; i++)
        insert_front(cells[i]);
}

void initGain(){
    for (int i = 0; i < ccnt; i++){
        cells[i]->gain = 0;
        cells[i]->lock = 0;
    }
    
    aGain = 0;
    bestg = aGain;
    bestacnnt = aCellCnt;
    bestbcnnt = bCellCnt;
    bestaCellsize = aCellsize;
    bestbCellsize = bCellsize;
    bestk = k;

    // for each cell, init its gain
    for (int i = 0; i < ccnt; i++){
        // for each net n on cell i, check if n is critical
        for (int j = 0 ; j < cells[i]->netList.size(); j++){
            int n = cells[i]->netList[j];
            // if cell i belongs to set A 
            if (cells[i]->set == 0) {
                if (nets[n]->A == 1) cells[i]->gain++;
                if (nets[n]->B == 0) cells[i]->gain--;
            }
            else {
                if (nets[n]->B == 1) cells[i]->gain++;
                if (nets[n]->A == 0) cells[i]->gain--;
            }
        }
    }
    buildbucketlist();
}

// find the cell with maximum gain
Cell * findMaxGain(bool set){
    int p = Pmax;
    // find the max gain (find the first list that is not empty)
    while (p >= -Pmax && bucketlist[set][p]->next == NULL){p--;}
    // find the first cell with maximum gain
    if (p<-Pmax){
        return NULL;
    }
    Cell * ans = cells[bucketlist[set][p]->next->id];
    return ans;
}

//update all the gain of all the cells that are related to the base cell c
void updateGain(Cell * c){
    aGain += c->gain;

    c->lock = true;
    int num = c->to->id;
    cellstack.push_back(num);
    // if c belong to set A
    if (!c->set) {
        int szn = c->netList.size();
        for(int i = 0; i < szn; i++){
            int id = c->netList[i];
            Net * net = nets[id];
            int szc = net->cellList.size();
            //check critical nets before the move
            // if T(n)=0 then increment gains of all free cells on n
            if (net->B == 0){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!cells[idc]->lock) {
                        cells[idc]->gain++;
                        move(cells[idc]);
                    }
                }
            }
            // else if T(n) = 1 then decrement gain of the only T cell on n,
            // if it is free
            else if (net->B == 1){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!cells[idc]->lock && cells[idc]->set) {
                        cells[idc]->gain--;
                        move(cells[idc]);
                    }
                }
            }
            // change F(n), T(n)
            net->A--;
            net->B++;
            c->set = true;
            // check for critical nets after the move
            // if F(n) = 0 then decrement gains of all free cells on n
            if (net->A == 0){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!cells[idc]->lock) {
                        cells[idc]->gain--;
                        move(cells[idc]);
                    }
                }
            }
            // else if F(n) = 1 then increment gain of the only F cell on n,
            // if it is free
            else if (net->A == 1){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!cells[idc]->lock && !cells[idc]->set) {
                        cells[idc]->gain++;
                        move(cells[idc]);
                    }
                }
            }
        }

        // remove the base cell
        remove(c);
        aCellsize -= c->size;
        bCellsize += c->size;
        aCellCnt--;
        bCellCnt++;
    }
    else {
        int szn = c->netList.size();
        for(int i = 0; i < szn; i++){
            int id = c->netList[i];
            Net * net = nets[id];
            int szc = net->cellList.size();
            if (net->A == 0){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!cells[idc]->lock) {
                        cells[idc]->gain++;
                        move(cells[idc]);
                    }
                }
            }
            else if (net->A == 1){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!cells[idc]->lock && !cells[idc]->set) {
                        cells[idc]->gain--;
                        move(cells[idc]);
                    }
                }
            }
            net->B--;
            net->A++;
            c->set = false;
            if (net->B == 0){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!cells[idc]->lock) {
                        cells[idc]->gain--;
                        move(cells[idc]);
                    }
                }
            }
            else if (net->B == 1){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!cells[idc]->lock && cells[idc]->set) {
                        cells[idc]->gain++;
                        move(cells[idc]);
                    }
                }
            }
        }
        remove(c);
        bCellsize -= c->size;
        aCellsize += c->size;
        bCellCnt--;
        aCellCnt++;
    }
    if (aGain > bestg){
        bestg = aGain;
        bestacnnt = aCellCnt;
        bestbcnnt = bCellCnt;
        bestaCellsize = aCellsize;
        bestbCellsize = bCellsize;
        bestk = k;
    }
        
    return;
}

int iter = 0;

// judge if the movement is valid or not
bool isValid(bool set,int sz){
    if (set){
        return abs(aCellsize-bCellsize-2*sz) < error;
    }
    else{
        return abs(bCellsize-aCellsize-2*sz) < error;
    }
}

void FMAlgorithm(){
    bool flag = false;
    initGain();
    int count = 0;
    k = 0;
    bestk = 0;
    cellstack.clear();
    while (!flag && count++ < ccnt){ 
        Cell * a = findMaxGain(0), * b = findMaxGain(1);
        int thred_n = (iter<=5)?3 : 10;
        
        if (a!=NULL && b!=NULL){
            if (a->gain >= b->gain) {
                int n = 0;
                
                while( !isValid(true,a->size) && a->to->next!=NULL && n++<=thred_n){
                    a = cells[a->to->next->id];
                }
                if (isValid(true,a->size)) updateGain(a);
                else if (isValid(false,b->size)) updateGain(b);
                else flag = true;
            }
            else {
                int n = 0;
                while( !isValid(false,b->size) && b->to->next!=NULL && n++<=thred_n){
                    b = cells[b->to->next->id];
                }
                if (isValid(false,b->size)) updateGain(b);
                else if (isValid(true,a->size)) updateGain(a);
                else flag = true;
            }
        }
        else if (a==NULL){
            int n = 0;
            while( !isValid(false,b->size) && b->to->next!=NULL ){
                b = cells[b->to->next->id];
            }
            if (isValid(false,b->size)) updateGain(b);
            else flag = true;
        }
        else {
            // check if balance
            int n = 0;
            while( !isValid(true,a->size) && a->to->next!=NULL){
                a = cells[a->to->next->id];
            }
            if (isValid(true,a->size)) updateGain(a);
            else flag = true;
        }
        k++;
    
    }
    
    if (bestg > 0 ) {
        iter++;
        
        k = bestk;
        aCellCnt = bestacnnt;
        bCellCnt = bestbcnnt;
        aCellsize = bestaCellsize;
        bCellsize = bestbCellsize;
        reverse();
        calAB();

        cout << "Iter " << iter << ", Gains: " << bestg << endl;
        if (iter>=40){
            return;
        }
        FMAlgorithm();
        
    }
    else { bestk = -1; k = -1; return;}
}


int main(int argc, char *argv[]){
   
    parse(argc, argv);

    // calculate intial cut size
    cutsz = 0;
    for (int i = 0; i < ncnt; i++)
        if (nets[i]->A && nets[i]->B) cutsz++;
    cout << "Initial Cut Size = " << cutsz << endl;
    cout << endl;

    error = (double) totalCellsize/10;
    for (int i = 0; i < ccnt; i++)
        if (cells[i]->pins > Pmax) Pmax = cells[i]->pins;

    bestaCellsize = aCellsize;
    bestbCellsize = bCellsize;
    tstart = clock();
    FMAlgorithm();
	tend = clock();

    k = bestk;
    aCellCnt = bestacnnt;
    bCellCnt = bestbcnnt;
    aCellsize = bestaCellsize;
    bCellsize = bestbCellsize;
    reverse();
    calAB();
    cutsz = 0;
    for (int i = 0; i < ncnt; i++)
        if (nets[i]->A && nets[i]->B) cutsz++;
    cout << "Final Cut Size = " << cutsz << endl;
    outputFile(of);
    of.close();
    cout << endl;
    cout << "FM Algorithm Run Time: " << (double)(tend-tstart)/CLOCKS_PER_SEC << " sec\n";
    cout << "Total Run Time: " << (double)clock()/CLOCKS_PER_SEC << " sec\n";
    
}
