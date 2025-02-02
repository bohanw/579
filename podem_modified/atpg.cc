/* stuck-at fault ATPG for combinational circuit
 * Last update: 2006/12/09 */
#include <iostream>
#include "circuit.h"
#include "GetLongOpt.h"
#include <cmath>
using namespace std;

vector<vector<char> > PIvector;
vector<vector<char> > Redun_PIvector;
int testsNum = 0;
int testsNumX = 0;
int TotalReorderedswitchActivity;
int BestReorderedSwitchActivit;
int NumOfX;
int TotalNumOfX = 0;
int NonRedundencyComb = 0;
vector<vector<GATE*> > TestPI;
vector<GATE*> PreviousNetlist;
int BeginCompare = 0;
int totalSwitchActivity = 0;
extern GetLongOpt option;

int test = 0;
//generate all stuck-at fault list
void CIRCUIT::GenerateAllFaultList()
{
    cout << "Generate stuck-at fault list" << endl;
    register unsigned i, j;
    GATEFUNC fun;
    GATEPTR gptr, fanout;
    FAULT *fptr;
    for (i = 0;i<No_Gate();++i) {
        gptr = Netlist[i]; fun = gptr->GetFunction();
        if (fun == G_PO) { continue; } //skip PO
        //add stem stuck-at 0 fault to Flist
        fptr = new FAULT(gptr, gptr, S0);
        Flist.push_front(fptr);
        //add stem stuck-at 1 fault to Flist
        fptr = new FAULT(gptr, gptr, S1);
        Flist.push_front(fptr);

        if (gptr->No_Fanout() == 1) { continue; } //no branch faults

        //add branch fault
        for (j = 0;j< gptr->No_Fanout();++j) {
            fanout = gptr->Fanout(j);
            fptr = new FAULT(gptr, fanout, S0);
            fptr->SetBranch(true);
            Flist.push_front(fptr);
            fptr = new FAULT(gptr, fanout, S1);
            fptr->SetBranch(true);
            Flist.push_front(fptr);
        } //end all fanouts
    } //end all gates
    //copy Flist to undetected Flist (for fault simulation)
    UFlist = Flist;
    return;
}

//stuck-at fualt PODEM ATPG (fault dropping)
void CIRCUIT::Atpg()
{
    cout << "Run stuck-at fault ATPG" << endl;
    unsigned i, total_backtrack_num(0), pattern_num(0);
    ATPG_STATUS status;
    FAULT* fptr;
    list<FAULT*>::iterator fite;
    
    //Prepare the output files
    ofstream OutputStrm;
    if (option.retrieve("output")){
        OutputStrm.open((char*)option.retrieve("output"),ios::out);
        if(!OutputStrm){
              cout << "Unable to open output file: "
                   << option.retrieve("output") << endl;
              cout << "Unsaved output!\n";
              exit(-1);
        }
    }

    if (option.retrieve("output")){
	    for (i = 0;i<PIlist.size();++i) {
		OutputStrm << "PI " << PIlist[i]->GetName() << " ";
	    }
	    OutputStrm << endl;
    }
    for (fite = Flist.begin(); fite != Flist.end();++fite) {
        fptr = *fite;
        if (fptr->GetStatus() == DETECTED) { continue; }
        //run podem algorithm
        status = Podem(fptr, total_backtrack_num);
        switch (status) {
            case TRUE:
                fptr->SetStatus(DETECTED);
                ++pattern_num;
                //run fault simulation for fault dropping
                for (i = 0;i < PIlist.size();++i) { 
			    ScheduleFanout(PIlist[i]); 
                if (option.retrieve("output")){ OutputStrm <<PIlist[i]->GetValue();}
		    }
                if (option.retrieve("output")){ OutputStrm << endl;}

                // save prvious netlist state
                //PreviousNetlist = Netlist;

                for (i = PIlist.size();i<Netlist.size();++i) {
                    Netlist[i]->SetValue(X);}

                LogicSim();
                FaultSim();
                break;
            case CONFLICT:
                fptr->SetStatus(REDUNDANT);
                break;
            case FALSE:
                fptr->SetStatus(ABORT);
                break;
        }

        //BeginCompare = 1;


    } //end all faults

    //compute fault coverage
    unsigned total_num(0);
    unsigned abort_num(0), redundant_num(0), detected_num(0);
    unsigned eqv_abort_num(0), eqv_redundant_num(0), eqv_detected_num(0);
    for (fite = Flist.begin();fite!=Flist.end();++fite) {
        fptr = *fite;
        switch (fptr->GetStatus()) {
            case DETECTED:
                ++eqv_detected_num;
                detected_num += 1;//fptr->GetEqvFaultNum();
                break;
            case REDUNDANT:
                ++eqv_redundant_num;
                redundant_num += 1;//fptr->GetEqvFaultNum();
                break;
            case ABORT:
                ++eqv_abort_num;
                abort_num += 1;//fptr->GetEqvFaultNum();
                break;
            default:
                cerr << "Unknown fault type exists" << endl;
                break;
        }
    }
    total_num = detected_num + abort_num + redundant_num;

    cout.setf(ios::fixed);
    cout.precision(2);
    cout << "---------------------------------------" << endl;
    cout << "Test pattern number = " << pattern_num << endl;
    cout << "Total backtrack number = " << total_backtrack_num << endl;
    cout << "---------------------------------------" << endl;
    cout << "Total fault number = " << total_num << endl;
    cout << "Detected fault number = " << detected_num << endl;
    cout << "Undetected fault number = " << abort_num + redundant_num << endl;
    cout << "Abort fault number = " << abort_num << endl;
    cout << "Redundant fault number = " << redundant_num << endl;
    cout << "---------------------------------------" << endl;
    cout << "Total equivalent fault number = " << Flist.size() << endl;
    cout << "Equivalent detected fault number = " << eqv_detected_num << endl;
    cout << "Equivalent undetected fault number = " << eqv_abort_num + eqv_redundant_num << endl;
    cout << "Equivalent abort fault number = " << eqv_abort_num << endl;
    cout << "Equivalent redundant fault number = " << eqv_redundant_num << endl;
    cout << "---------------------------------------" << endl;
    cout << "Fault Coverge = " << 100*detected_num/double(total_num) << "%" << endl;
    cout << "Equivalent FC = " << 100*eqv_detected_num/double(Flist.size()) << "%" << endl;
    cout << "Fault Efficiency = " << 100*detected_num/double(total_num - redundant_num) << "%" << endl;
    cout << "---------------------------------------" << endl;
    return;
}

//run PODEM for target fault
//TRUE: test pattern found
//CONFLICT: no test pattern
//FALSE: abort
ATPG_STATUS CIRCUIT::Podem(FAULT* fptr, unsigned &total_backtrack_num)
{
    unsigned i, backtrack_num(0);
    GATEPTR pi_gptr(0), decision_gptr(0);
    ATPG_STATUS status;

    //set all values as unknown
    for (i = 0;i<Netlist.size();++i) { Netlist[i]->SetValue(X); }
    //mark propagate paths
    MarkPropagateTree(fptr->GetOutputGate());
    //propagate fault free value
    status = SetUniqueImpliedValue(fptr);
    switch (status) {
        case TRUE:
            LogicSim();
            //inject faulty value
            if (FaultEvaluate(fptr)) {
                //forward implication
                ScheduleFanout(fptr->GetOutputGate());
                LogicSim();
            }
            //check if the fault has propagated to PO
            if (!CheckTest()) { status = FALSE; }
            break;
        case CONFLICT:
            status = CONFLICT;
            break;
        case FALSE: break;
    }

    while(backtrack_num < BackTrackLimit && status == FALSE) {
        //search possible PI decision
        pi_gptr = TestPossible(fptr);
        if (pi_gptr) { //decision found
            ScheduleFanout(pi_gptr);
            //push to decision tree
            GateStack.push_back(pi_gptr);
            decision_gptr = pi_gptr;
        }
        else { //backtrack previous decision
            while (!GateStack.empty() && !pi_gptr) {
                //all decision tried (1 and 0)
                if (decision_gptr->GetFlag(ALL_ASSIGNED)) {
                    decision_gptr->ResetFlag(ALL_ASSIGNED);
                    decision_gptr->SetValue(X);
                    ScheduleFanout(decision_gptr);
                    //remove decision from decision tree
                    GateStack.pop_back();
                    decision_gptr = GateStack.back();
                }
                //inverse current decision value
                else {
                    decision_gptr->InverseValue();
                    ScheduleFanout(decision_gptr);
                    decision_gptr->SetFlag(ALL_ASSIGNED);
                    ++backtrack_num;
                    pi_gptr = decision_gptr;
                }
            }
            //no other decision
            if (!pi_gptr) { status = CONFLICT; }
        }
        if (pi_gptr) {
            LogicSim();
            //fault injection
            if(FaultEvaluate(fptr)) {
                //forward implication
                ScheduleFanout(fptr->GetOutputGate());
                LogicSim();
            }
            if (CheckTest()) { status = TRUE; }
        }
    } //end while loop

    //clean ALL_ASSIGNED and MARKED flags
    list<GATEPTR>::iterator gite;
    for (gite = GateStack.begin();gite != GateStack.end();++gite) {
        (*gite)->ResetFlag(ALL_ASSIGNED);
    }
    for (gite = PropagateTree.begin();gite != PropagateTree.end();++gite) {
        (*gite)->ResetFlag(MARKED);
    }

    //clear stacks
    GateStack.clear(); PropagateTree.clear();
    
    //assign true values to PIs
    vector<int> Xindex;
    if (status ==  TRUE) {
        vector<char> temp;
        vector<char> temp2;
        NumOfX = 0;
        PIvector.push_back(temp);
        Redun_PIvector.push_back(temp2);
		for (i = 0;i<PIlist.size();++i) {
		    switch (PIlist[i]->GetValue()) {
			case S1: {PIvector[testsNum].push_back('1'); Redun_PIvector[testsNumX].push_back('1');break;}
			case S0: {PIvector[testsNum].push_back('0'); Redun_PIvector[testsNumX].push_back('0');break;}
			case D: {PIvector[testsNum].push_back('1'); Redun_PIvector[testsNumX].push_back('1'); PIlist[i]->SetValue(S1); break;}
			case B: {PIvector[testsNum].push_back('0'); Redun_PIvector[testsNumX].push_back('0'); PIlist[i]->SetValue(S0); break;}
            // set unknown value for experiment
			case X: PIvector[testsNum].push_back('0'); Redun_PIvector[testsNumX].push_back('0'); NumOfX = NumOfX + 1;Xindex.push_back(i); TotalNumOfX = TotalNumOfX + 1; PIlist[i]->SetValue(S0);break;// PIlist[i]->SetValue(VALUE(2.0 * rand()/(RAND_MAX + 1.0))); break;
			default: cerr << "Illigal value" << endl; break;
		    } 
		}//end for all PI
        testsNum = testsNum + 1;
        testsNumX = testsNumX + 1;
        NonRedundencyComb = NonRedundencyComb + 1;
    } //end status == TRUE



// X-filling to 0 and 1

// if too many X in a vector - random filling
if (status ==  TRUE) {
cout << "NumOfX =" << NumOfX << ',';
test = test + 1;
cout << test <<endl;
int NumOfBits = PIlist.size();
if((NumOfX > NumOfBits/2) || (NumOfX > 10)){

    int ParentTestIndex_1 = testsNumX - 1;
    Redun_PIvector.push_back(Redun_PIvector[ParentTestIndex_1]);
        for(int k = 0; k < Xindex.size(); k++){
            Redun_PIvector[testsNumX][Xindex[k]] = '1';
        }
        
        testsNumX = testsNumX + 1;
        
    // for(int k = 0; k < Xindex.size(); k++){
    //     Redun_PIvector[testsNumX - 1][Xindex[k]] = (char)(VALUE(2.0 * rand()/(RAND_MAX + 1.0)) + 48);
    // }
} 
else {
//if not too many X in a vector - exhasusted all combanation
    int comb = pow(2,NumOfX);
    int ParentTestIndex = testsNumX - 1;
    for(int j = 1; j < comb; j++){
        Redun_PIvector.push_back(Redun_PIvector[ParentTestIndex]);
        for(int k = 0; k < Xindex.size(); k++){
            Redun_PIvector[testsNumX][Xindex[k]] = (char)(((j >> k) & 1) + 48);
        }
        testsNumX = testsNumX + 1;
 //       cout << "testsNumX"<< testsNumX <<endl;
    }
}

}
/////////////////////////////////////////////////////////////////////

    total_backtrack_num += backtrack_num;
    return status;
}


vector<vector<char> >& CIRCUIT::getRedunPI(){
    //cout << "redsize = " << Redun_PIvector.size();
    return Redun_PIvector;
}
void CIRCUIT::printParameters(){
    cout << "NumOfX = " << TotalNumOfX << endl;
    cout << "NonRedundencyComb = " << NonRedundencyComb << endl;
    cout << "testsNum = "<< testsNum << endl;
}

// print the assigned input in termimal

void CIRCUIT::printPI(vector<vector<char> >& PIprint){
    for(int i = 0; i < PIprint.size(); i++){
        cout<<"test"<< i << '=';
        for(int j = 0; j < PIprint[i].size(); j++){
            cout << PIprint[i][j]<<',';
        }

        cout << endl;
    }
}

vector<vector<char> >& CIRCUIT::getPIvector(){
    return PIvector;
}
// calculate switch activity

int CIRCUIT::CalSwitchActivity(vector<vector<char> >& testvector){
    int count = 0;
    int size1 = testvector.size();
    //cout << "size1 = "<<size1<<endl;
    for(int i = 1; i < size1; i++){
        for(int j = 0; j < testvector[i].size() && j < testvector[i-1].size(); j++){
           // cout << "i = " <<  testvector[i].size() <<',' << " j=" <<testvector[i-1].size()<<endl;
            if(testvector[i][j] != testvector[i-1][j]) count = count + 1;
        }
    }

    return count;
}


vector<vector<char> > CIRCUIT::greedy(vector<vector<char> >& vectorin){
    int NumOfTest = vectorin.size();
    TotalReorderedswitchActivity = 0;
    int switchActivity;
    for(int i = 0; i < NumOfTest - 1; i++){
        switchActivity = NumOfTest;
        //cout<<"process test" << i << endl;
        int nextvector = 0;
        for(int j = i+1; j < NumOfTest; j++){
            int count = 0;
            for(int k = 0; k < vectorin[i].size(); k++){
                if(vectorin[i][k] != vectorin[j][k]) count = count + 1;
            }

            if(count <= switchActivity) {
                switchActivity = count;
                nextvector = j;
            }
        }

        TotalReorderedswitchActivity = TotalReorderedswitchActivity + switchActivity;
        vector<char> temp;
        temp = vectorin[i+1];
        vectorin[i+1] = vectorin[nextvector];
        vectorin[nextvector] = temp;
    }
    // for(int m = 0; m < vectorin.size(); m++){
    //     for(int n = 0; n < vectorin[m].size(); n++){
    //         cout <<vectorin[m][n];
    //     }
    //     cout <<' '<<endl;
    // }
    //cout << "reordered switchActivity=" << TotalReorderedswitchActivity << endl;
    return vectorin;
}

vector<vector<char> > CIRCUIT::reorder(vector<vector<char> >& vectorin){
    int NumOfTest = vectorin.size();

    vector<vector<char> > best = greedy(vectorin);

    // for(int m = 0; m < best.size(); m++){
    //     for(int n = 0; n < best[m].size(); n++){
    //         cout <<best[m][n];
    //     }
    //     cout <<' '<<endl;
    // }

    //cout << "TotalReorderedswitchActivity = "<<TotalReorderedswitchActivity << endl;
    BestReorderedSwitchActivit = TotalReorderedswitchActivity;
    //cout<<"BestReorderedSwitchActivit = "<< BestReorderedSwitchActivit<< endl;

    for(int i = 1; i < NumOfTest; i++){
        //cout<<"process test" << i << endl;
        vector<char> temp_v = vectorin[i];
        vectorin[i] = vectorin[0];
        vectorin[0] = temp_v;
        vector<vector<char> > temp2 = greedy(vectorin);
        //cout<<"BestReorderedSwitchActivit = "<< BestReorderedSwitchActivit<< endl;
        if(TotalReorderedswitchActivity < BestReorderedSwitchActivit){
            BestReorderedSwitchActivit = TotalReorderedswitchActivity;
            best = temp2;
        }

    }


    for(int m = 0; m < best.size(); m++){
        for(int n = 0; n < best[m].size(); n++){
            cout <<best[m][n];
        }
        cout <<' '<<endl;
    }
    cout<<"BestReorderedSwitchActivit = "<< BestReorderedSwitchActivit<< endl;
    return best;
    // choose the start vector

}

int CIRCUIT::ComputeTotalsw(vector<vector<char> >& vectorin){
    int testNum = vectorin.size();
    if(testNum == 0) return 0;

    for (int j = 0;j<vectorin[0].size();++j){
        switch(vectorin[0][j]){
            case '0': PIlist[j]->SetValue(S0);break;
            case '1': PIlist[j]->SetValue(S1);break;
            default: cout<<"what??"<<endl; break;}
        ScheduleFanout(PIlist[j]);}

    for (int k = PIlist.size();k<Netlist.size();++k) {Netlist[k]->SetValue(X);}
        LogicSim();


    vector<GATE*> OldNetlist = Netlist;

    int Total_sw = 0;
    for(int i = 1; i < testNum; i++){
    for (int j = 0;j<vectorin[i].size();++j){
        switch(vectorin[i][j]){
            case '0': PIlist[j]->SetValue(S0);break;
            case '1': PIlist[j]->SetValue(S1);break;
            default: cout<<"what??"<<endl; break;}
        ScheduleFanout(PIlist[j]);}
    for (int k = PIlist.size();k<Netlist.size();++k) {Netlist[k]->SetValue(X);}
        LogicSim();

    for (int k = 0;k<Netlist.size();++k) {

        if(OldNetlist[k]->GetValue() == Netlist[k]->GetValue()) Total_sw++;
    }

    OldNetlist = Netlist;

    }

    return Total_sw;
}

double CIRCUIT::ComputeFaultCoverage(vector<vector<char> >& vectorin){
    int testNum = vectorin.size();
    if(testNum == 0) return 0.0;
    //cout << "testNum" << testNum<<endl;
    list<FAULT*>::iterator fite;
    FAULT* fptr;
    for (fite = Flist.begin();fite!=Flist.end();++fite) {
        fptr = *fite;
        fptr->SetStatus(UNKNOWN);
    }

    UFlist = Flist;


    for(int i = 0; i < testNum; i++){
        for (int j = 0;j<vectorin[i].size();++j){
            switch(vectorin[i][j]){
                case '0': PIlist[j]->SetValue(S0);break;
                case '1': PIlist[j]->SetValue(S1);break;
                default: cout<<"what??"<<endl; break;
            }

            ScheduleFanout(PIlist[j]); 
        }
        for (int k = PIlist.size();k<Netlist.size();++k) {Netlist[k]->SetValue(X);}
            LogicSim();
            FaultSim();
    }

    int remain_fault = 0;
    for (fite = UFlist.begin();fite!=UFlist.end();++fite) {
        remain_fault ++;
    }
    //compute fault coverage
    unsigned total_num(0);
    unsigned abort_num(0), redundant_num(0), detected_num(0), unknown_fault(0);
    unsigned eqv_abort_num(0), eqv_redundant_num(0), eqv_detected_num(0);

    int debug = 0;

    for (fite = Flist.begin();fite!=Flist.end();++fite) {
        fptr = *fite;
        debug = debug + 1;
        switch (fptr->GetStatus()) {
            case DETECTED:
                ++eqv_detected_num;
                detected_num += 1;//fptr->GetEqvFaultNum();
                break;
            case REDUNDANT:
                ++eqv_redundant_num;
                redundant_num += 1;//fptr->GetEqvFaultNum();
                break;
            case ABORT:
                ++eqv_abort_num;
                abort_num += 1;//fptr->GetEqvFaultNum();
                break;
            default:
                //cerr << "Unknown fault type exists" << endl;
                unknown_fault++;
                break;
        }
    }
    total_num = detected_num + abort_num + redundant_num + unknown_fault;
    double FaultCoverge = detected_num/double(total_num);

    // cout << "**************************************" << endl;
    // cout << "Total fault number = " << total_num << endl;
    // cout << "Detected fault number = " << detected_num << endl;
    // cout << "Undetected fault number = " << abort_num + redundant_num << endl;
    // cout << "Abort fault number = " << abort_num << endl;
    // cout << "Redundant fault number = " << redundant_num << endl;
    // cout << "**************************************" << endl;
    // cout << "Fault Coverge = " << 100*detected_num/double(total_num) << "%" << endl;
    // cout << "Equivalent FC = " << 100*eqv_detected_num/double(Flist.size()) << "%" << endl;
    // cout << "Fault Efficiency = " << 100*detected_num/double(total_num - redundant_num) << "%" << endl;
    // cout << "**************************************" << endl;

    // cout << "debug = " << debug<<endl;
    // cout << "remain_fault = " << remain_fault << endl;

    return FaultCoverge;

}

//inject fault-free value and do simple backward implication
//TRUE: fault can be backward propagated to PI
//CONFLICT: conflict assignment
//FALSE: fault can not be propagated to PI
ATPG_STATUS CIRCUIT::SetUniqueImpliedValue(FAULT* fptr)
{
    register ATPG_STATUS status(FALSE);
    GATEPTR igptr(fptr->GetInputGate());
    //backward implication fault-free value
    switch (BackwardImply(igptr, NotTable[fptr->GetValue()])) {
        case TRUE: status = TRUE; break;
        case CONFLICT: return CONFLICT; break;
        case FALSE: break;
    }
    if (!fptr->Is_Branch()) { return status; }
    //if branch, need to check side inputs of the output gate
    GATEPTR ogptr(fptr->GetOutputGate());
    VALUE ncv(NCV[ogptr->GetFunction()]);
    //set side inputs as non-controlling value
    for (unsigned i = 0;i < ogptr->No_Fanin();++i) {
        if (igptr == ogptr->Fanin(i)) { continue; }
        switch (BackwardImply(ogptr->Fanin(i), ncv)) {
            case TRUE: status = TRUE; break;
            case CONFLICT: return CONFLICT; break;
            case FALSE: break;
        }
    }
    return status;
}

//apply the input values of gate according to its output value
//TRUE: fault can be backward propagated to PI
//CONFLICT: conflict assignment
//FALSE: fault can not be propagated to PI
ATPG_STATUS CIRCUIT::BackwardImply(GATEPTR gptr, VALUE value)
{
    register unsigned i;
    register ATPG_STATUS status(FALSE);
    GATEFUNC fun(gptr->GetFunction());

    if (fun == G_PI) { //reach PI
        //conflict assignment
        if (gptr->GetValue() != X && gptr->GetValue() != value) {
            return CONFLICT;
        }
        gptr->SetValue(value); ScheduleFanout(gptr);
        return TRUE;
    }
    //not PI gate
    switch (fun) {
        case G_NOT:
            switch (BackwardImply(gptr->Fanin(0), NotTable[value])) {
                case TRUE: status = TRUE; break;
                case FALSE: break;
                case CONFLICT: return CONFLICT; break;
            }
            break;
        case G_BUF:
            switch (BackwardImply(gptr->Fanin(0), value)) {
                case TRUE: status = TRUE; break;
                case FALSE: break;
                case CONFLICT: return CONFLICT; break;
            }
            break;
        case G_AND:
        case G_OR:
            if (value != NCV[fun]) { break; }
            for (i = 0;i<gptr->No_Fanin();++i) {
                switch (BackwardImply(gptr->Fanin(i), NCV[fun])) {
                    case TRUE: status = TRUE; break;
                    case FALSE: break;
                    case CONFLICT: return CONFLICT; break;
                }
            }
            break;
        case G_NAND:
        case G_NOR:
            if (value != CV[fun]) { break; }
            for (i = 0;i<gptr->No_Fanin();++i) {
                switch (BackwardImply(gptr->Fanin(i), NCV[fun])) {
                    case TRUE: status = TRUE; break;
                    case FALSE: break;
                    case CONFLICT: return CONFLICT; break;
                }
            }
            break;
        default: break;
    }
    return status;
}

//mark and push propagate tree to stack PropagateTree
void CIRCUIT::MarkPropagateTree(GATEPTR gptr)
{
    PropagateTree.push_back(gptr);
    gptr->SetFlag(MARKED);
    for (unsigned i = 0;i<gptr->No_Fanout();++i) {
        if (gptr->Fanout(i)->GetFlag(MARKED)) { continue; }
        MarkPropagateTree(gptr->Fanout(i));
    }
    return;
}

//fault injection
//true: fault is injected successfully and need to do fault propagation
//false: output value is the same with original one or fault is injected in PO
bool CIRCUIT::FaultEvaluate(FAULT* fptr)
{
    GATEPTR igptr(fptr->GetInputGate());
    //store input value
    VALUE ivalue(igptr->GetValue());
    //can not do fault injection
    if (ivalue == X || ivalue == fptr->GetValue()) { return false; }
    else if (ivalue == S1) { igptr->SetValue(D); }
    else if (ivalue == S0) { igptr->SetValue(B); }
    else { return false; } //fault has been injected

    if (!fptr->Is_Branch()) { return true; }
    //for branch fault, the fault has to be propagated to output gate
    GATEPTR ogptr(fptr->GetOutputGate());
    if (ogptr->GetFunction() == G_PO) { return false; }
    VALUE value(Evaluate(ogptr));
    //backup original value to input gate
    igptr->SetValue(ivalue);
    //fault has propagated to output gate
    if (value != ogptr->GetValue()) {
        ogptr->SetValue(value);
        return true;
    }
    return false;
}

//return possible PI decision
GATEPTR CIRCUIT::TestPossible(FAULT* fptr)
{
    GATEPTR decision_gptr;
    GATEPTR ogptr(fptr->GetOutputGate());
    VALUE decision_value;
    if (!ogptr->GetFlag(OUTPUT)) {
        if (ogptr->GetValue() != X) {
            //if no fault injected, return 0
            if (ogptr->GetValue() != B && ogptr->GetValue() != D) { return 0; }
            //search D-frontier
            decision_gptr = FindPropagateGate();
            if (!decision_gptr) { return 0;}
            switch (decision_gptr->GetFunction()) {
                case G_AND:
                case G_NOR: decision_value = S1; break;
                case G_NAND:
                case G_OR: decision_value = S0; break;
                default: return 0;
            }
        }
        else { //output gate == X
            //test if any unknown path can propagate the fault
            if (!TraceUnknownPath(ogptr)) { return 0; }
            if (!fptr->Is_Branch()) { //stem
                decision_value = NotTable[fptr->GetValue()];
                decision_gptr = ogptr;
            }
            else { //branch
                //output gate value is masked by side inputs
                if (fptr->GetInputGate()->GetValue() != X) {
                    switch (ogptr->GetFunction()) {
                        case G_AND:
                        case G_NOR: decision_value = S1; break;
                        case G_NAND:
                        case G_OR: decision_value = S0; break;
                        default: return 0;
                    }
                    decision_gptr = fptr->GetOutputGate();
                }
                //both input and output values are X
                else {
                    decision_value = NotTable[fptr->GetValue()];
                    decision_gptr = fptr->GetInputGate();
                }
            } //end branch
        } //end output gate == X
    } //end if output gate is PO
    else { //reach PO
        if (fptr->GetInputGate()->GetValue() == X) {
            decision_value = NotTable[fptr->GetValue()];
            decision_gptr = fptr->GetInputGate();
        }
        else { return 0; }
    }
    return FindPIAssignment(decision_gptr, decision_value);
}

//find PI decision to set gptr = value
//success: return PI
//fail: return 0
GATEPTR CIRCUIT::FindPIAssignment(GATEPTR gptr, VALUE value)
{
    //search PI desicion
    if (gptr->GetFunction() == G_PI) {
        gptr->SetValue(value);
        return gptr;
    }
    GATEPTR j_gptr(0); //J-frontier
    VALUE j_value(X), cv_out;
    switch (gptr->GetFunction()) {
        case G_AND:
        case G_NAND:
        case G_OR:
        case G_NOR:
            cv_out = CV[gptr->GetFunction()];
            cv_out = gptr->Is_Inversion()? NotTable[cv_out]: cv_out;
            //select one fanin as cv
            if (value == cv_out) {
                j_gptr = FindEasiestControl(gptr);
                j_value = CV[gptr->GetFunction()];
            }
            //select one fanin as ncv 
            else {
                j_gptr = FindHardestControl(gptr);
                j_value = NCV[gptr->GetFunction()];
            }
            break;
        case G_BUF:
        case G_NOT:
            j_value = gptr->Is_Inversion()? NotTable[value]: value;
            j_gptr = gptr->Fanin(0);
            break;
        default:
            break;
    }
    if (j_gptr) { return FindPIAssignment(j_gptr, j_value); }
    return 0;
}

//check if the fault has propagated to PO
bool CIRCUIT::CheckTest()
{
    VALUE value;
    for (unsigned i = 0;i<POlist.size();++i) {
        value = POlist[i]->GetValue();
        if (value == B || value == D) { return true; }
    }
    return false;
}

//search gate from propagate tree to propagate the fault
GATEPTR CIRCUIT::FindPropagateGate()
{
    register unsigned i;
    list<GATEPTR>::iterator gite;
    GATEPTR gptr, fanin;
    for (gite = PropagateTree.begin();gite!=PropagateTree.end();++gite) {
        gptr = *gite;
        if (gptr->GetValue() != X) { continue; }
        for (i = 0;i<gptr->No_Fanin(); ++i) {
            fanin = gptr->Fanin(i);
            if (fanin->GetValue() != D && fanin->GetValue() != B) { continue; }
            if (TraceUnknownPath(gptr)) { return gptr; }
            break;
        }
    }
    return 0;
}

//trace if any unknown path from gptr to PO
bool CIRCUIT::TraceUnknownPath(GATEPTR gptr)
{
    if (gptr->GetFlag(OUTPUT)) { return true; }
    GATEPTR fanout;
    for (unsigned i = 0;i<gptr->No_Fanout();++i) {
        fanout = gptr->Fanout(i);
        if (fanout->GetValue()!=X) { continue; }
        if (TraceUnknownPath(fanout)) { return true; }
    }
    return false;
}

//serach lowest level unknown fanin
GATEPTR CIRCUIT::FindEasiestControl(GATEPTR gptr)
{
    GATEPTR fanin;
    for (unsigned i = 0;i< gptr->No_Fanin();++i) {
        fanin = gptr->Fanin(i);
        if (fanin->GetValue() == X) { return fanin; }
    }
    return 0;
}

//serach highest level unknown fanin
GATEPTR CIRCUIT::FindHardestControl(GATEPTR gptr)
{
    GATEPTR fanin;
    for (unsigned i = gptr->No_Fanin();i>0;--i) {
        fanin = gptr->Fanin(i-1);
        if (fanin->GetValue() == X) { return fanin; }
    }
    return 0;
}

//functor, used to compare the levels of two gates
struct sort_by_level
{
    bool operator()(const GATEPTR gptr1, const GATEPTR gptr2) const
    {
        return (gptr1->GetLevel() < gptr2->GetLevel());
    }
};

//sort fanin by level (from low to high)
void CIRCUIT::SortFaninByLevel()
{
    for (unsigned i = 0;i<Netlist.size();++i) {
        vector<GATE*> &Input_list = Netlist[i]->GetInput_list();
        sort(Input_list.begin(), Input_list.end(), sort_by_level());
    }
    return;
}

//mark and trace all stem fault propagated by fault "val"
void CIRCUIT::TraceDetectedStemFault(GATEPTR gptr, VALUE val)
{
    //no PO stem fault
    if (gptr->GetFunction() == G_PO) { return; }
    //get output fault type
    if (gptr->Is_Inversion()) { val = NotTable[val]; }
    //translate value to flag: S0->ALL_ASSIGNED, S1->MARKED
    FLAGS flag = FLAGS(val);
    if (val != S1 && val != S0) { cerr << "error" << endl; }
    //stem fault has been detected
    if (gptr->GetFlag(flag)) { return; }
    gptr->SetFlag(flag);
    //stop when the gate meets branch
    if (gptr->No_Fanout() > 1) { return; }
    TraceDetectedStemFault(gptr->Fanout(0), val);
    return;
}
