#include <iostream>
#include <ctime>
#include "circuit.h"
#include "GetLongOpt.h"
#include "ReadPattern.h"
#include "genetic.h"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
#include <sstream>
#include <vector>
#include <random>

#include <unordered_map>



using namespace std;

// All defined in readcircuit.l
extern char* yytext;
extern FILE *yyin;
extern CIRCUIT Circuit;
extern int yyparse (void);
extern bool ParseError;

extern void Interactive();

GetLongOpt option;



int SetupOption(int argc, char ** argv)
{
    option.usage("[options] input_circuit_file");
    option.enroll("help", GetLongOpt::NoValue,
            "print this help summary", 0);
    option.enroll("logicsim", GetLongOpt::NoValue,
            "run logic simulation", 0);
    option.enroll("plogicsim", GetLongOpt::NoValue,
            "run parallel logic simulation", 0);
    option.enroll("fsim", GetLongOpt::NoValue,
            "run stuck-at fault simulation", 0);
    option.enroll("stfsim", GetLongOpt::NoValue,
            "run single pattern single transition-fault simulation", 0);
    option.enroll("transition", GetLongOpt::NoValue,
            "run transition-fault ATPG", 0);
    option.enroll("input", GetLongOpt::MandatoryValue,
            "set the input pattern file", 0);
    option.enroll("output", GetLongOpt::MandatoryValue,
            "set the output pattern file", 0);
    option.enroll("bt", GetLongOpt::OptionalValue,
            "set the backtrack limit", 0);
    int optind = option.parse(argc, argv);
    if ( optind < 1 ) { exit(0); }
    if ( option.retrieve("help") ) {
        option.usage();
        exit(0);
    }
    return optind;
}

int main(int argc, char ** argv)
{
    try_1();
    int optind = SetupOption(argc, argv);
    clock_t time_init, time_end;
    time_init = clock();
    //Setup File
    if (optind < argc) {
        if ((yyin = fopen(argv[optind], "r")) == NULL) {
            cout << "Can't open circuit file: " << argv[optind] << endl;
            exit( -1);
        }
        else {
            string circuit_name = argv[optind];
            string::size_type idx = circuit_name.rfind('/');
            if (idx != string::npos) { circuit_name = circuit_name.substr(idx+1); }
            idx = circuit_name.find(".bench");
            if (idx != string::npos) { circuit_name = circuit_name.substr(0,idx); }
            Circuit.SetName(circuit_name);
        }
    }
    else {
        cout << "Input circuit file missing" << endl;
        option.usage();
        return -1;
    }
    cout << "Start parsing input file\n";
    yyparse();
    if (ParseError) {
        cerr << "Please correct error and try Again.\n";
        return -1;
    }
    fclose(yyin);
    Circuit.FanoutList();
    Circuit.SetupIO_ID();
    Circuit.Levelize();
    Circuit.Check_Levelization();
    Circuit.InitializeQueue();

    if (option.retrieve("logicsim")) {
        //logic simulator
        Circuit.InitPattern(option.retrieve("input"));
        Circuit.LogicSimVectors();
    }
    else if (option.retrieve("plogicsim")) {
        //parallel logic simulator
        Circuit.InitPattern(option.retrieve("input"));
        Circuit.ParallelLogicSimVectors();
    }
    else if (option.retrieve("stfsim")) {
        //single pattern single transition-fault simulation
        Circuit.MarkOutputGate();
        Circuit.GenerateAllTFaultList();
        Circuit.InitPattern(option.retrieve("input"));
        Circuit.TFaultSimVectors();
    }
    else if (option.retrieve("transition")) {
        Circuit.MarkOutputGate();
        Circuit.GenerateAllTFaultList();
        Circuit.SortFaninByLevel();
        if (option.retrieve("bt")) {
            Circuit.SetBackTrackLimit(atoi(option.retrieve("bt")));
        }
        Circuit.TFAtpg();
    }
    else {
        Circuit.GenerateAllFaultList();
        Circuit.SortFaninByLevel();
        Circuit.MarkOutputGate();
        if (option.retrieve("fsim")) {
            //stuck-at fault simulator
            Circuit.InitPattern(option.retrieve("input"));
            Circuit.FaultSimVectors();
        }

        else {
            if (option.retrieve("bt")) {
                Circuit.SetBackTrackLimit(atoi(option.retrieve("bt")));
            }
            //stuck-at fualt ATPG
            Circuit.Atpg();
            vector<vector<char> >& PIvector = Circuit.getPIvector();
            vector<vector<char> >& Redun_PIvector = Circuit.getRedunPI();
            cout << "redsize = " << Redun_PIvector.size();
            cout<<"-------------------------------------------------"<<endl;
            Circuit.printPI(PIvector);
            Circuit.printPI(Redun_PIvector);

            double gold_faultCoverage = Circuit.ComputeFaultCoverage(PIvector);

            int xfilling_ratio = Redun_PIvector.size()/PIvector.size();

            //Circuit.printParameters();
            //int count = Circuit.CalSwitchActivity(PIvector);
            //cout<<"count = "<<count<<endl;

            //vector<vector<char> >& reorderedPI = Circuit.reorder(PIvector);
            //Circuit.printPI(reorderedPI);

            
// genetic algorithm - first phase
//******************************************************************************

            //initialization:
            // PIvector_pool is the gene pool with total X genes
            // generate POPSIZE individules
            // the number of gene each individual has is gaussian random distribution N ~ (X/2, X/4)
            // gene is selected randomly from the gene pool and assigned to each individual

            string phase1_evolution = "result_s298_opt_phase1_evolution";
            string phase2_evolution = "result_s298_opt_phase2_evolution";
            srand(time(NULL));
            cout.precision(5);
            int pop_size = POPSIZE;
            int good_individual_num = 0;
            vector<vector<char> >& PIvector_pool = Redun_PIvector;
            int seed = 123456789;
            initialize(seed, PIvector_pool,xfilling_ratio);
//*******************************************************************************
            // compute fault converage for each individual
            for(int m = 0; m < POPSIZE; m++){
                cout << "compute fc "<<m<<endl;
                double fc = Circuit.ComputeFaultCoverage(population[m].gene);
                population[m].fault_coverage = fc;
                // in first phase, the fitness of an individual is its fault_coverage
                population[m].fitness = fc;
            }
            cout<<endl;

            for (int m = 0; m < POPSIZE; m++)
               {
                Oldpopulation[m] = population[m];
               }

            for(int m = 0; m < POPSIZE; m++){
                cout<< "faultCoverage = " <<population[m].fault_coverage<<endl;
            }
            cout << "gold_faultCoverage = " << gold_faultCoverage << endl;
 

 //
                keep_the_best();

//evolution starts ------------------------------------------------------------------
    ofstream OutputStrm;
    OutputStrm.open(phase1_evolution, ofstream::out);
    if(!OutputStrm){
        cout << "Unable to open output file: " << endl;
        cout << "Unsaved output!\n";
        exit(-1);}


    for (int generation = 0; generation < MAXGENS; generation++)
    {
        good_individual_num = 0;
        int individual_to_be_killed = 0;//= selector_p1(seed,pop_size);

        vector<int> newborn_individual = crossover_p1(seed,individual_to_be_killed);
        vector<int> mutated_individual = mutate_p1(seed, Redun_PIvector);
        //report(generation, PIvector);

        //cout << "POPSIZE=" << POPSIZE<< endl;
            for(int m = 0; m < pop_size; m++){
                double fc_old = Oldpopulation[m].fault_coverage; //Circuit.ComputeFaultCoverage(Oldpopulation[m].gene);
                double fc = Circuit.ComputeFaultCoverage(population[m].gene);


                if(fc_old>fc) {
                    population[m] = Oldpopulation[m];
                    population[m].fault_coverage = fc_old;
                    population[m].fitness = fc_old;
                }
                else{
                    population[m].fault_coverage = fc;
                    // in first phase, the fitness of an individual is its fault_coverage
                    population[m].fitness = fc;   
                    Oldpopulation[m] = population[m];               
                }


                if(population[m].fault_coverage >= gold_faultCoverage){
                    good_individual_num ++;
                }
            }
        //evaluate(PIvector);
        //elitist();

        double largest_fc = population[0].fitness;
        int size_p1 = population[0].gene.size();
        double sum_fc = population[0].fitness;
        for(int m = 1; m < pop_size; m++){
            sum_fc = sum_fc + population[m].fitness;
                if(largest_fc < population[m].fitness) {
                    largest_fc = population[m].fitness;
                    size_p1 = population[m].gene.size();
                }
        }
 
        cout<<"run generation: "<<generation<< " largest_fc = "<< largest_fc << ',' << "average_fc = " << (double)sum_fc/pop_size <<','<<"gene_size = "<< size_p1<<endl;

        OutputStrm <<"run generation: "<<generation<< " largest_fc = "<< largest_fc << ',' << "average_fc = " << (double)sum_fc/pop_size <<','<<"gene_size = "<< size_p1<<endl;
        if(good_individual_num == pop_size) break;
    }    

 // add base genes to individual

    for(int m = 0; m < pop_size; m++){
        for(int n = 0; n < PIvector.size(); n++){
            population[m].gene.push_back(PIvector[n]);
            Oldpopulation[m].gene.push_back(PIvector[n]);
        }
    }  

    for(int m = 0; m < POPSIZE; m++){
        double fc = Circuit.ComputeFaultCoverage(population[m].gene);
        population[m].fault_coverage = fc;
        }     

    cout << setw(14) << population[POPSIZE].fitness;
    for(int m = 0; m < pop_size; m++){
        cout << m << " gene size = "<< population[m].gene.size()<<' ';
        cout << m << " fault coverage = "<< population[m].fault_coverage<<endl;

        //cout << m << " gene size = "<< Oldpopulation[m].gene.size()<<' ';
        //cout << m << " fault coverage = "<< Oldpopulation[m].fault_coverage<<endl;  
        // for(int n = 0; n < population[m].gene.size(); n++){
        //     for(int k = 0; k < population[m].gene[n].size(); k++){
        //         //cout<<population[m].gene[n][k];
        //     }
        //     cout <<' ';
        // }
        cout<<endl;

        OutputStrm << m << " gene size = "<< population[m].gene.size()<<' ';
        OutputStrm << m << " fault coverage = "<< population[m].fault_coverage<<endl;
    }
//evolution ends ------------------------------------------------------------------  



 OutputStrm.close();

//phase 2
    for(int m = 0; m < pop_size; m++){
        population[m].fitness = Circuit.CalSwitchActivity(population[m].gene);
        cout << population[m].fitness << endl;
    }
    //population[pop_size] = population[0];
    //keep_the_best_p2();

    for (int m = 0; m < POPSIZE; m++){
        Oldpopulation[m] = population[m];}

    ofstream OutputStrm2;

    OutputStrm2.open(phase2_evolution, ofstream::out);
        if(!OutputStrm){
        cout << "Unable to open output file: " << endl;
        cout << "Unsaved output!\n";
        exit(-1);}

    int final_tests = 0;
    for (int generation = 0; generation < MAXGENS_p2; generation++)
    {
        good_individual_num = 0;
        int individual_to_be_killed = 0;//= selector_p1(seed,pop_size);

        //vector<int> newborn_individual = crossover_p1(seed,individual_to_be_killed);
        vector<int> mutated_individual = mutate_p2(seed, PIvector_pool);
        //report(generation, PIvector);
        //report(generation, PIvector);

            for(int m = 0; m < pop_size; m++){

                //cout << m<<endl;
                //cout << "*******";
                //cout << Oldpopulation[m].gene.size() << ',' << population[m].gene.size()<<endl;
                double fc_old = Oldpopulation[m].fault_coverage; //Circuit.ComputeFaultCoverage(Oldpopulation[m].gene);
                double fc = Circuit.ComputeFaultCoverage(population[m].gene);

                int sw_old = Oldpopulation[m].fitness; //Circuit.CalSwitchActivity(Oldpopulation[m].gene);
                int sw = Circuit.CalSwitchActivity(population[m].gene);

                if(fc < gold_faultCoverage) {
                    population[m] = Oldpopulation[m];
                    population[m].fault_coverage = fc_old;
                    population[m].fitness = sw_old;
                }
                else if(sw_old < sw){
                    population[m] = Oldpopulation[m];
                    population[m].fault_coverage = fc_old;
                    population[m].fitness = sw_old;                 
                }
                else{
                    population[m].fault_coverage = fc;
                    population[m].fitness = sw;  
                    Oldpopulation[m] = population[m];                 
                }

            }

            //selector_p2(seed, pop_size);

            int smallest_sw = population[0].fitness;
            int size = population[0].gene.size();
            double f = population[0].fault_coverage;
            double sum_sw = population[0].fitness;
            int index = 0;
            for(int m = 1; m < pop_size; m++){
                sum_sw = sum_sw + population[m].fitness;
            if(smallest_sw > population[m].fitness) {
                smallest_sw = population[m].fitness;
                size = population[m].gene.size();
                f = population[m].fault_coverage;
                index = m;
            }
        }

        cout<<"run generation: "<<generation<< " index = "<< index <<" smallest_sw = "<< smallest_sw << ',' <<" average sw =" <<(double)sum_sw/pop_size << ','<< "gene_size = "<< size<<" fault_coverage =" << f <<endl;
        OutputStrm2 << cout<<"run generation: "<<generation<< " index = "<< index <<" smallest_sw = "<< smallest_sw << ',' <<" average sw =" <<(double)sum_sw/pop_size << ','<< "gene_size = "<< size<<" fault_coverage =" << f <<endl;
        
        final_tests = index;
        //evaluate(PIvector);
        //elitist_p2();

       // if(good_individual_num == pop_size) break;
    } 
        cout << "final test vector =";
        OutputStrm2 << "final test vector =";
    for(int m = 0; m < population[final_tests].gene.size(); m++){
        for(int n = 0; n < population[final_tests].gene[m].size(); n++){
            cout <<population[final_tests].gene[m][n];
            OutputStrm2<<population[final_tests].gene[m][n];
        }
        cout <<' ';
        OutputStrm2 <<' ';

    }
    cout << endl;
    OutputStrm2<<endl;

// reordering
    vector<vector<char> > reorderedPI = Circuit.reorder(population[final_tests].gene);
    int final_sw = Circuit.CalSwitchActivity(reorderedPI);
    double final_fc = Circuit.ComputeFaultCoverage(reorderedPI);

    int ini_sw = Circuit.CalSwitchActivity(PIvector);
    int ini_tests_size = PIvector.size();

    cout << "reordered sw = " << final_sw << endl;
    cout << "reordered fc = " << final_fc << endl;

    cout << "sw before algorithm = " << ini_sw << endl;
    cout << "initial test number = " << ini_tests_size << endl;

    OutputStrm2 << "reordered sw = " << final_sw << endl;
    OutputStrm2 << "reordered fc = " << final_fc << endl;    

    OutputStrm2 << "sw before algorithm = " << ini_sw << endl;
    OutputStrm2 << "initial test number = " << ini_tests_size << endl;

    for(int m = 0; m < reorderedPI.size(); m++){
        for(int n = 0; n < reorderedPI[m].size(); n++){
            cout <<reorderedPI[m][n] ;
            OutputStrm2<<reorderedPI[m][n] <<' ';
        }
        cout <<' '<<endl;
        OutputStrm2 <<' '<<endl;
    }
    cout << endl;
    int Total_sw_after = Circuit.ComputeTotalsw(reorderedPI);
    int Total_sw_before = Circuit.ComputeTotalsw(PIvector);
    cout <<"Total_sw_after = "<<Total_sw_after << " ,Total_sw_before = " << Total_sw_before << endl;
    OutputStrm2 <<"Total_sw_after = "<<Total_sw_after << " ,Total_sw_before = " << Total_sw_before << endl;
    OutputStrm2.close();

// output file
            // ofstream OutputStrm;

            // OutputStrm.open("credundantFile", ofstream::out);
            // if(!OutputStrm){
            //   cout << "Unable to open output file: " << endl;
            //   cout << "Unsaved output!\n";
            //   exit(-1);}

            //   for(int i = 0; i < PIvector.size(); i++){
            //     //cout<<PIvector[i].size()<<endl;
            //     for(int j = 0; j < PIvector[i].size(); j++){

            //         OutputStrm<<PIvector[i][j]<<',';
            //     }

            //     OutputStrm<<endl;
            //   }

            //   OutputStrm.close();
///////////////////////////////////////////////////////////
        }
    }
    time_end = clock();
    cout << "Total CPU Time = " << double(time_end - time_init)/CLOCKS_PER_SEC << endl;
    cout << endl;


    return 0;
}
