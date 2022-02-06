#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <ctype.h>
#include <set>
#include <map>
#include <list>
#include <vector>
#include <cmath>
#include "net.h"
#include "cell.h"
#include <ctime>

using namespace std;

typedef vector <int> vi;

vector <Net*> vn, *bestvn;
vector <Cell*> vc, *bestvc;
vector <int> cellstack;
int k = 0, bestk;
bool *bestset;
vector <int> bestA, bestB;
map <string, int> mc, mn; 
int ccnt = 0, ncnt = 0;

int cutSize = 0;
ifstream ifc, ifn;
ofstream of;
double error; 
int totalCellsize = 0, aCellsize = 0, bCellsize = 0, cs = 0;
int aCellCnt = 0, bCellCnt = 0, aGain = 0, bestg = 0;

int bestacnnt = 0, bestbcnnt = 0, bestaCellsize = 0, bestbCellsize = 0;
int Pmax = 0;
double tstart, tend;
map <int, Node*> blist[2];

void parseCells(istream & in){
    string str;
    int size;
    while (in >> str >> size){
        mc[str] = ccnt;
        if (aCellsize <= bCellsize){
            Cell *c = new Cell(str, size, 0, ccnt);
            vc.push_back(c);
            aCellsize += size;
            aCellCnt++;
        }
        else {
            Cell *c = new Cell(str, size, 1, ccnt);
            vc.push_back(c);
            bCellsize += size;
            bCellCnt++;
        }
        totalCellsize += size; 
        ccnt++;
    }
}

void calAB(){
    for (int i = 0; i < ncnt; i++){
        vi & vl = vn[i]->cellList;
        vn[i]->B = 0;
        vn[i]->A = 0;
        for (int j = 0; j < vl.size(); j++){
            Cell * cell = vc[vl[j]];
            if (cell->set) vn[i]->B++;
            else vn[i]->A++;
        }
    }
}

void parseNets(istream & in){
    string str, tmp;
    while (in >> tmp){ // NET
        in >> str;  // nxxx
        mn[str] = ncnt;
        in >> tmp;  // {
        Net *n = new Net(str);
        vn.push_back(n);
        while (in >> tmp && tmp[0] != '}'){
            vi & l = vc[mc[tmp]]->netList;
            if (!l.size() || l[l.size()-1] != ncnt) {
                l.push_back(ncnt);
                vc[mc[tmp]]->pins++;
                vn[ncnt]->cellList.push_back(mc[tmp]);
                if (vc[mc[tmp]]->set) vn[ncnt]->B++;
                else vn[ncnt]->A++;
            }
        }
        ncnt++;
    }
}

void countCutSize(){
    cs = 0;
    for (int i = 0; i < ncnt; i++)
        if (vn[i]->A && vn[i]->B) cs++;
}


void test(){
    for (int i = 0; i < ccnt; i++){
        cout << vc[i]->name << " s"
	     << vc[i]->size << " g"
             << vc[i]->gain << " ab"
	     << vc[i]->set  << " p"
	     << vc[i]->pins << " l"
	     << vc[i]->lock << ' ';
	for (int j = 0; j < vc[i]->netList.size(); j++){
	    int id = vc[i]->netList[j];
	    cout << vn[id]->name << ' ';
        }
	cout << endl;
    }
    cout << "...\n";
    cout << "# of cells = " << ccnt << endl;
    cout << "Total Cell Size = " << totalCellsize << endl;
    cout << "Set A Size = " << aCellsize << endl;
    cout << "Set B Size = " << bCellsize << endl;
    cout << "|A - B| = " << abs(aCellsize-bCellsize) << endl;
    cout << "...\n";
    for (int i = 0; i < ncnt; i++){
        cout << vn[i]->name << ' ';
        for (int j = 0; j < vn[i]->cellList.size(); j++){
            int id = vn[i]->cellList[j];
            cout << vc[id]->name << ' ';
        }
        cout << endl;
    }
    cout << "...\n";
    return;
}

//calculate Pmax to intialize the bucket list
void countPmax(){
    for (int i = 0; i < ccnt; i++)
        if (vc[i]->pins > Pmax) Pmax = vc[i]->pins;
}

// error = n/10
void countError(){
    error = (double) totalCellsize/10;
}

void outputFile(ostream & out){
    out << "cut_size " << cs << endl;
    out << "A " << aCellCnt << endl;
    for (int i = 0; i < ccnt; i++)
        if (!vc[i]->set)
            out << vc[i]->name << endl;
    out << "B " << bCellCnt << endl;
    for (int i = 0; i < ccnt; i++)
        if (vc[i]->set)
            out << vc[i]->name << endl;
}


void parseInput(int argc, char ** argv){

    char opt = 0;
    while ((opt = getopt(argc, argv, "c:n:o:h?")) != -1){
        switch (opt){
            case 'c':
                ifc.open(optarg, ios::in);
                if (!ifc.is_open())
                    cout << "Cannot open the cells file at [-" << opt << ' ' << optarg << ']' << endl;		
                break;
            case 'n':
                ifn.open(optarg, ios::in);
                if (!ifn.is_open())
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


void traverse(){
    
    for (int k = 0; k < 2; k++){
        cout << "---- " << ((!k) ? "A" : "B") << " ----\n";
        for (int i = Pmax ; i >= -Pmax; i--){
            cout << '[' << i << ']' << ' ';
            Node *trav = blist[k][i]->next;
            while (trav != NULL){
                cout << vc[trav->id]->name << "->";
                trav = trav->next;
            }
            cout << endl;
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
    Node *pre = blist[set][gain];
    Node* cur = blist[set][gain]->next;
    //while(cur!=NULL && vc[cur->id]->size<c->size){
    //    cout<<p->id<<" "<<c->size<<" "<< vc[cur->id]->size<<endl;
    //    pre = pre->next;
    //    cur = cur->next;
    //}
    //cout<<endl;
    p->prev = pre;
    p->next = cur;
    pre->next = p;
    
    //p->prev = blist[set][gain];
    //p->next = blist[set][gain]->next;
    //blist[set][gain]->next = p;
    if (p->next != NULL) p->next->prev = p;
}

void move(Cell * c){
    remove(c);
    insert_front(c);
}

// build the bucket list
void buildBlist(){
    blist[0].clear();
    blist[1].clear();
    // init the bucket list( head node for all p)
    for (int i = -Pmax; i <= Pmax; i++) {
        if (blist[0][i] == NULL) blist[0][i] = new Node(-1);
        if (blist[1][i] == NULL) blist[1][i] = new Node(-1);
    }
    // insert all the cells to the bucket list
    for (int i = 0; i < ccnt; i++)
        insert_front(vc[i]);
}

// find the cell with maximum gain
Cell * findMaxGain(bool set){
    int p = Pmax;
    // find the max gain (find the first list that is not empty)
    while (p >= -Pmax && blist[set][p]->next == NULL){p--;}
    // find the first cell with maximum gain
    if (p<-Pmax){
        return NULL;
    }
    Cell * ans = vc[blist[set][p]->next->id];
    return ans;
}


void reverse(){
    int i = cellstack.size()-1;
    for (; i > k; i--)
        //cout << cellstack[i] << " ";
        vc[cellstack[i]]->set = !vc[cellstack[i]]->set;
    
}

void store(){
    bestg = aGain;
    //bestvc = &vc;
    //bestvn = &vn;
    bestacnnt = aCellCnt;
    bestbcnnt = bCellCnt;
    bestaCellsize = aCellsize;
    bestbCellsize = bCellsize;
    //bestset.clear();
    //bestA.clear();
    //bestB.clear();
    bestk = k;
    //for (int i = 0 ; i < ccnt; i++)
    //    bestset[i] = vc[i]->set;
}



void restore(){
    //vc = *bestvc;
    //vn = *bestvn;
    k = bestk;
    //cout << k;
    aCellCnt = bestacnnt;
    bCellCnt = bestbcnnt;
    aCellsize = bestaCellsize;
    bCellsize = bestbCellsize;
    //for (int i = 0 ; i < ccnt; i++)
    //    vc[i]->set = bestset[i];
    reverse();
    calAB();
    //cout << "???\n";
}


void initGain(){
    for (int i = 0; i < ccnt; i++){
        vc[i]->gain = 0;
        vc[i]->lock = 0;
    }
    
    aGain = 0;
    store();

    // for each cell, init its gain
    for (int i = 0; i < ccnt; i++){
        // for each net n on cell i, check if n is critical
        for (int j = 0 ; j < vc[i]->netList.size(); j++){
            int n = vc[i]->netList[j];
            // if cell i belongs to set A 
            if (vc[i]->set == 0) {
                if (vn[n]->A == 1) vc[i]->gain++;
                if (vn[n]->B == 0) vc[i]->gain--;
            }
            else {
                if (vn[n]->B == 1) vc[i]->gain++;
                if (vn[n]->A == 0) vc[i]->gain--;
            }
        }
    }
    buildBlist();
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
            Net * net = vn[id];
            int szc = net->cellList.size();
            //check critical nets before the move
            // if T(n)=0 then increment gains of all free cells on n
            if (net->B == 0){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock) {
                        vc[idc]->gain++;
                        move(vc[idc]);
                    }
                }
            }
            // else if T(n) = 1 then decrement gain of the only T cell on n,
            // if it is free
            else if (net->B == 1){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock && vc[idc]->set) {
                        vc[idc]->gain--;
                        move(vc[idc]);
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
                    if (!vc[idc]->lock) {
                        vc[idc]->gain--;
                        move(vc[idc]);
                    }
                }
            }
            // else if F(n) = 1 then increment gain of the only F cell on n,
            // if it is free
            else if (net->A == 1){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock && !vc[idc]->set) {
                        vc[idc]->gain++;
                        move(vc[idc]);
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
            Net * net = vn[id];
            int szc = net->cellList.size();
            if (net->A == 0){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock) {
                        vc[idc]->gain++;
                        move(vc[idc]);
                    }
                }
            }
            else if (net->A == 1){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock && !vc[idc]->set) {
                        vc[idc]->gain--;
                        move(vc[idc]);
                    }
                }
            }
            net->B--;
            net->A++;
            c->set = false;
            if (net->B == 0){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock) {
                        vc[idc]->gain--;
                        move(vc[idc]);
                    }
                }
            }
            else if (net->B == 1){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock && vc[idc]->set) {
                        vc[idc]->gain++;
                        move(vc[idc]);
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
                    a = vc[a->to->next->id];
                }
                if (abs(aCellsize-bCellsize-2*a->size) < error) updateGain(a);
                else if (abs(bCellsize-aCellsize-2*b->size) < error) updateGain(b);
                else flag = true;
            }
            else {
                int n = 0;
                while(abs(bCellsize-aCellsize-2*b->size) >= error && b->to->next!=NULL && n++<=thred_n){
                    b = vc[b->to->next->id];
                }
                if (abs(bCellsize-aCellsize-2*b->size) < error) updateGain(b);
                else if (abs(aCellsize-bCellsize-2*a->size) < error) updateGain(a);
                else flag = true;
            }
        }
        else if (a==NULL){
            int n = 0;
            while(abs(bCellsize-aCellsize-2*b->size) >= error && b->to->next!=NULL ){
                b = vc[b->to->next->id];
            }
            if (abs(bCellsize-aCellsize-2*b->size) < error) updateGain(b);
            else flag = true;
        }
        else {
            // check if balance
            int n = 0;
            while(abs(aCellsize-bCellsize-2*a->size) >= error && a->to->next!=NULL){
                a = vc[a->to->next->id];
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
            Cell * c = vc[i];
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
    if (ifc.is_open()) parseCells(ifc);
    else parseCells(cin);
    ifc.close();

    if (ifn.is_open()) parseNets(ifn);
    else parseNets(cin);
    ifn.close();
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
