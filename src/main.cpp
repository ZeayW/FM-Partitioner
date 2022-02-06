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
int totalCellsize = 0, aCellsize = 0, bCellsize = 0, cs = 0;
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

void countCutSize(){
    cs = 0;
    for (int i = 0; i < ncnt; i++)
        if (nets[i]->A && nets[i]->B) cs++;
}


//calculate Pmax to intialize the bucket list
void countPmax(){
    for (int i = 0; i < ccnt; i++)
        if (cells[i]->pins > Pmax) Pmax = cells[i]->pins;
}

// error = n/10
void countError(){
    error = (double) totalCellsize/10;
}

void outputFile(ostream & out){
    out << "cut_size " << cs << endl;
    out << "A " << aCellCnt << endl;
    for (int i = 0; i < ccnt; i++)
        if (!cells[i]->set)
            out << cells[i]->name << endl;
    out << "B " << bCellCnt << endl;
    for (int i = 0; i < ccnt; i++)
        if (cells[i]->set)
            out << cells[i]->name << endl;
}


void parseInput(int argc, char ** argv){

    char opt = 0;
    while ((opt = getopt(argc, argv, "c:n:o:h?")) != -1){
        switch (opt){
            case 'c':
                ifsCell.open(optarg, ios::in);
                if (!ifsCell.is_open())
                    cout << "Cannot open the cells file at [-" << opt << ' ' << optarg << ']' << endl;		
                break;
            case 'n':
                ifsNet.open(optarg, ios::in);
                if (!ifsNet.is_open())
                    cout << "Cannot open the nets file at [-" << opt << ' ' << optarg << ']' << endl;		
                break;
            case 'o':
                of.open(optarg, ios::out);
                if (!of.is_open()) {
                    cout << "Cannot open the output file at [-" << opt << ' ' << optarg << ']' << endl;		
                }
                break;
            case 'h' :
			case '?' :
            default:
                cerr << "Usage: " << argv[0] << " -c <cells file name> -n <nets file name> -o <output file name>\n";
                exit(EXIT_FAILURE);
        }
    }

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


void reverse(){
    int i = cellstack.size()-1;
    for (; i > k; i--)
        cells[cellstack[i]]->set = !cells[cellstack[i]]->set;
    
}

void store(){
    bestg = aGain;
    bestacnnt = aCellCnt;
    bestbcnnt = bCellCnt;
    bestaCellsize = aCellsize;
    bestbCellsize = bCellsize;
    bestk = k;

}



void restore(){

    k = bestk;
    aCellCnt = bestacnnt;
    bCellCnt = bestbcnnt;
    aCellsize = bestaCellsize;
    bCellsize = bestbCellsize;

    reverse();
    calAB();

}


void initGain(){
    for (int i = 0; i < ccnt; i++){
        cells[i]->gain = 0;
        cells[i]->lock = 0;
    }
    
    aGain = 0;
    store();

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
    if (aGain > bestg)
        store();
    return;
}

int pass = 0;

bool isValid(int set,int sz){
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
        // 
        Cell * a = findMaxGain(0), * b = findMaxGain(1);
        int thred_n = (pass<=5)?3 : 10;
        
        if (a!=NULL && b!=NULL){
            if (a->gain >= b->gain) {
                int n = 0;
                
                while(abs(aCellsize-bCellsize-2*a->size) >= error && a->to->next!=NULL && n++<=thred_n){
                    a = cells[a->to->next->id];
                }
                if (isValid(0,a->size)) updateGain(a);
                else if (isValid(1,b->size)) updateGain(b);
                else flag = true;
            }
            else {
                int n = 0;
                while(abs(bCellsize-aCellsize-2*b->size) >= error && b->to->next!=NULL && n++<=thred_n){
                    b = cells[b->to->next->id];
                }
                if (abs(bCellsize-aCellsize-2*b->size) < error) updateGain(b);
                else if (abs(aCellsize-bCellsize-2*a->size) < error) updateGain(a);
                else flag = true;
            }
        }
        else if (a==NULL){
            int n = 0;
            while(abs(bCellsize-aCellsize-2*b->size) >= error && b->to->next!=NULL ){
                b = cells[b->to->next->id];
            }
            if (abs(bCellsize-aCellsize-2*b->size) < error) updateGain(b);
            else flag = true;
        }
        else {
            // check if balance
            int n = 0;
            while(abs(aCellsize-bCellsize-2*a->size) >= error && a->to->next!=NULL){
                a = cells[a->to->next->id];
            }
            if (abs(aCellsize-bCellsize-2*a->size) < error) updateGain(a);
            else flag = true;
        }
        k++;
    
    }
    
    if (bestg > 0 ) {
        pass++;
        
        restore();
        cout << "Pass " << pass << endl;
        cout << "Best Partial Sum of Gains: " << bestg << endl;
        cout << "Total Sum of Gains (Should be 0): " << aGain << endl;
        cout << endl;
        if (pass>=40){
            return;
        }
        FMAlgorithm();
        
    }
    else { bestk = -1; k = -1; return;}
}

// adjust the sets to be balanced
void adjust(){
    if (abs(aCellsize-bCellsize) < error) return;
    else {
        cout << "...Need balancing.\n";
        int i;
        for (i = 0; i < ccnt && abs(aCellsize-bCellsize) >= error; i++){
            Cell * c = cells[i];
            // if |A|>|B| and c belong to A , then move c to B
            if (aCellsize > bCellsize && !c->set){
                aCellsize -= c->size;
                bCellsize += c->size;
                c->set = true;
            }
            else if (aCellsize < bCellsize && c->set){
                aCellsize += c->size;
                bCellsize -= c->size;
                c->set = false;
            }
        }
        if (i == ccnt && abs(aCellsize-bCellsize) >= error) {
            cerr << "(ERROR)...This testcase can never be balanced!\n";
            //exit(EXIT_FAILURE);
        }
    }
}


int main(int argc, char *argv[]){
    cout<<"hello"<<endl;
    ios_base::sync_with_stdio(false);
	
    parseInput(argc, argv);
    if (ifsCell.is_open()) parseCells(ifsCell);
    else parseCells(cin);
    ifsCell.close();

    if (ifsNet.is_open()) parseNets(ifsNet);
    else parseNets(cin);
    ifsNet.close();
    //bestset = new bool[ccnt]();
    countCutSize();
    cout << "Initial Cut Size = " << cs << endl;
    cout << endl;
    countError();
    countPmax();
    adjust();
    
    bestaCellsize = aCellsize;
    bestbCellsize = bCellsize;
    tstart = clock();
    FMAlgorithm();
	tend = clock();
    restore();
    countCutSize();
    cout << "Final Cut Size = " << cs << endl;
    if (of.is_open()) outputFile(of);
    else outputFile(cout);
    of.close();
    cout << endl;
    cout << "FM Algorithm Run Time: " << (double)(tend-tstart)/CLOCKS_PER_SEC << " sec\n";
    cout << "Total Run Time: " << (double)clock()/CLOCKS_PER_SEC << " sec\n";
    
}
