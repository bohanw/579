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
#include <map>

using namespace std;
// 
//  Change any of these parameters to match your needs 
//
# define POPSIZE 20
# define MAXGENS 1000
# define MAXGENS_p2 1000
# define NVARS 30
# define DISTR 10
# define PXOVER 0.8
# define PMUTATION 0.15
# define PMUTATION_P2_DELET 0.5
//
//  Each GENOTYPE is a member of the population, with
//  gene: a string of variables,
//  fitness: the fitness
//  upper: the variable upper bounds,
//  lower: the variable lower bounds,
//  rfitness: the relative fitness,
//  cfitness: the cumulative fitness.
//
struct genotype
{
  vector<vector<char> > gene;
  double fitness;
  double d_fitness;
  double fault_coverage;
  double upper[NVARS];
  double lower[NVARS];
  double rfitness;
  double cfitness;
  int size;
};
struct genotype population[POPSIZE + 1];
struct genotype Oldpopulation[POPSIZE + 1];

struct genotype newpopulation[POPSIZE + 1];

vector<int> crossover_p1(int &seed, int individual_to_be_killed);
void elitist();
void elitist_p2();
void evaluate(vector<vector<char> > &PIvector);
int i4_uniform_ab(int a, int b, int &seed);
void initialize(int &seed, vector<vector<char> > &PIvector,int xfilling_ratio);
void keep_the_best();
void keep_the_best_p2();
vector<int> mutate_p1(int &seed, vector<vector<char> > &PIvector);
vector<int> mutate_p2(int &seed, vector<vector<char> > &PIvector);
double r8_uniform_ab(double a, double b, int &seed);
void report(int generation, vector<vector<char> > &PIvector);
int selector_p1(int &seed, int pop_size);
void selector_p2(int &seed, int pop_size);
void timestamp();
void Xover_p1(int one, int two, int &seed, int individual_to_be_killed);

void try_1();

void try_1(){
  cout << "try_2"<<endl;}


vector<int> crossover_p1(int &seed, int individual_to_be_killed)
{
  const double a = 0.0;
  const double b = 1.0;
  int mem;
  int one;
  int first = 0;
  double x;

  vector<int> newborn_individual;

      for (mem = 0; mem < POPSIZE; ++mem){
        x = r8_uniform_ab(a, b, seed);

        if (x < PXOVER){
          ++first;
          //cout << "newborn_individual 1"<<endl;
          if (first % 2 == 0){

            Xover_p1(one, mem, seed, individual_to_be_killed);

          //cout << "newborn_individual 2"<<endl;
            newborn_individual.push_back(one);
          }
          else{
            one = mem;
          }

        }
      }
  // for(int i = 0; i < newborn_individual.size(); i++){
  //   cout<< "new=" <<newborn_individual[i] << endl;
  // }
  return newborn_individual;
}



//****************************************************************************80

void elitist()

//****************************************************************************80
// 
//  Purpose:
//
//    ELITIST stores the best member of the previous generation.
//
//  Discussion:
//
//    The best member of the previous generation is stored as 
//    the last in the array. If the best member of the current 
//    generation is worse then the best member of the previous 
//    generation, the latter one would replace the worst member 
//    of the current population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//
//    Local, double WORST, the worst fitness value.
//
{
  int i;
  double best;
  int best_mem;
  double worst;
  int worst_mem;

  best = population[0].fitness;
  worst = population[0].fitness;

  for (i = 0; i < POPSIZE - 1; ++i)
  {
    if (population[i + 1].fitness < population[i].fitness)
    {

      if (best <= population[i].fitness)
      {
        best = population[i].fitness;
        best_mem = i;
      }

      if (population[i + 1].fitness <= worst)
      {
        worst = population[i + 1].fitness;
        worst_mem = i + 1;
      }

    }
    else
    {

      if (population[i].fitness <= worst)
      {
        worst = population[i].fitness;
        worst_mem = i;
      }

      if (best <= population[i + 1].fitness)
      {
        best = population[i + 1].fitness;
        best_mem = i + 1;
      }

    }

  }
  // 
  //  If the best individual from the new population is better than 
  //  the best individual from the previous population, then 
  //  copy the best from the new population; else replace the 
  //  worst individual from the current population with the 
  //  best one from the previous generation                     
  //
  if (population[POPSIZE].fitness <= best)
  {
    population[POPSIZE].gene.resize(population[best_mem].gene.size());
    for (i = 0; i < population[best_mem].gene.size(); i++)
    {
      population[POPSIZE].gene = population[best_mem].gene;
    }
    population[POPSIZE].fitness = population[best_mem].fitness;
  }
  else
  {
    population[worst_mem].gene.resize(population[POPSIZE].gene.size());
    for (i = 0; i < population[POPSIZE].gene.size(); i++)
    {
      population[worst_mem].gene = population[POPSIZE].gene;
    }
    population[worst_mem].fitness = population[POPSIZE].fitness;
  }

  return;
}
//****************************************************************************80

void elitist_p2()

//****************************************************************************80
// 
//  Purpose:
//
//    ELITIST stores the best member of the previous generation.
//
//  Discussion:
//
//    The best member of the previous generation is stored as 
//    the last in the array. If the best member of the current 
//    generation is worse then the best member of the previous 
//    generation, the latter one would replace the worst member 
//    of the current population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//
//    Local, double WORST, the worst fitness value.
//
{
  int i;
  double best;
  int best_mem;
  double worst;
  int worst_mem;

  best = population[0].fitness;
  worst = population[0].fitness;

  for (i = 0; i < POPSIZE - 1; ++i)
  {
    if (population[i + 1].fitness > population[i].fitness)
    {

      if (best <= population[i].fitness)
      {
        best = population[i].fitness;
        best_mem = i;
      }

      if (population[i + 1].fitness <= worst)
      {
        worst = population[i + 1].fitness;
        worst_mem = i + 1;
      }

    }
    else
    {

      if (population[i].fitness <= worst)
      {
        worst = population[i].fitness;
        worst_mem = i;
      }

      if (best <= population[i + 1].fitness)
      {
        best = population[i + 1].fitness;
        best_mem = i + 1;
      }

    }

  }
  // 
  //  If the best individual from the new population is better than 
  //  the best individual from the previous population, then 
  //  copy the best from the new population; else replace the 
  //  worst individual from the current population with the 
  //  best one from the previous generation                     
  //
  if (population[POPSIZE].fitness > best)
  {
    population[POPSIZE].gene.resize(population[best_mem].gene.size());
    for (i = 0; i < population[best_mem].gene.size(); i++)
    {
      population[POPSIZE].gene = population[best_mem].gene;
    }
    population[POPSIZE].fitness = population[best_mem].fitness;
  }
  else
  {
    population[worst_mem].gene.resize(population[POPSIZE].gene.size());
    for (i = 0; i < population[POPSIZE].gene.size(); i++)
    {
      population[worst_mem].gene = population[POPSIZE].gene;
    }
    population[worst_mem].fitness = population[POPSIZE].fitness;
  }

  return;
}

void evaluate(vector<vector<char> > &PIvector)

//****************************************************************************80
// 
//  Purpose:
//
//    EVALUATE implements the user-defined valuation function
//
//  Discussion:
//
//    Each time this is changed, the code has to be recompiled.
//    The current function is:  x[1]^2-x[1]*x[2]+x[3]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
{
  int member;
  int geneX;
  int gene_loc;
  int PI_num = PIvector[0].size();
  int i;



  for (member = 0; member < POPSIZE; member++)
  {
    int switching = 0;
    for (gene_loc = 0; gene_loc < PI_num; gene_loc++)
    {
      for (geneX = 0; geneX < population[member].gene.size()-1; geneX++)
      {
        switching = switching + !( population[member].gene[geneX][gene_loc] == population[member].gene[geneX+1][gene_loc]);
      }
    }
    population[member].fitness = switching;
    population[member].size = population[member].gene.size();

  }

  return;
}
    
    //  for (i = 0; i < population[member].gene.size(); i++)
  //  {
  //    x[i + 1] = population[member].gene[i];
  //  }
  //  switching = 0;
  //  for (i = 0; i < population[member].gene.size()-1; i++)
  //  {
  //    switching = switching + !(x[i + 1]== x[i]);
  //  }


//****************************************************************************80

int i4_uniform_ab(int a, int b, int &seed)

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if (seed == 0)
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit(1);
  }
  //
  //  Guarantee A <= B.
  //
  if (b < a)
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * (seed - k * 127773) - k * 2836;

  if (seed < 0)
  {
    seed = seed + i4_huge;
  }

  r = (float)(seed) * 4.656612875E-10;
  //
  //  Scale R to lie between A-0.5 and B+0.5.
  //
  r = (1.0 - r) * ((float)a - 0.5)
    + r   * ((float)b + 0.5);
  //
  //  Use rounding to convert R to an integer between A and B.
  //
  value = round(r);
  //
  //  Guarantee A <= VALUE <= B.
  //
  if (value < a)
  {
    value = a;
  }
  if (b < value)
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

void initialize(int &seed, vector<vector<char> > &PIvector, int xfilling_ratio)
{
  int i;
  ifstream input;
  int j;
  double lbound;
  double ubound;

  int genpool_size = PIvector.size();

  // input.open(filename.c_str());

  // if (!input)
  // {
  //   cerr << "\n";
  //   cerr << "INITIALIZE - Fatal error!\n";
  //   cerr << "  Cannot open the input file!\n";
  //   exit(1);
  // }


/*  string str;
  while (getline(input, str))
  {
    istringstream ss(str);
    int num;
    while (ss >> num)
    {
      cout << num << " ";
    }
    cout << std::endl;
  }
*/
  std::default_random_engine gen;
  //std::mt19937 gen(rd());

  // values near the mean are the most likely
  // standard deviation affects the dispersion of generated values from the mean
  std::normal_distribution<double> d(genpool_size/(2*xfilling_ratio), genpool_size / (4*xfilling_ratio));

  vector<int> store;
  for (i = 0; i < 1000; i++)
  {
    store.push_back(d(gen));
  }

  //std::map<int, int> hist;
  //for (int n = 0; n<10000; ++n) {
  //  ++hist[std::round(d(gen))];
  //}

  // 
  //  Initialize variables within the bounds 
  //

//  for (i = 0; i < NVARS; i++)
//  {
//    input >> lbound >> ubound;

    for (j = 0; j < POPSIZE; j++)
    {
      input >> lbound >> ubound;
      population[j].fitness = 0;
      population[j].d_fitness = 0;
      population[j].rfitness = 0;
      population[j].cfitness = 0;
//      population[j].size = rand() % ((NVARS+ DISTR) - (NVARS - DISTR) + 1) + (NVARS - DISTR);
      double x = d(gen);
      while(x < 1.0) x = d(gen);
      population[j].size = (int)x;

      population[j].size = genpool_size/xfilling_ratio;
      for (i = 0; i < population[j].size; i++)
      {
        int genpool_num = rand() % genpool_size;
        population[j].gene.push_back(PIvector[genpool_num]);
      //  population[j].lower[i] = lbound;
      //  population[j].upper[i] = ubound;
      //  int randomBit = rand() % 2;
      //  population[j].gene.push_back(rand() % 2);
      }
    }


    //  population[j].gene[i] = r8_uniform_ab(lbound, ubound, seed);

    
  input.close();

  return;
}
//****************************************************************************80

void keep_the_best()

{
  int cur_best;
  int mem;
  int i;
  int previous_size;
  cur_best = 0;

  for (mem = 0; mem < POPSIZE; mem++)
  {
    if (population[POPSIZE].fitness < population[mem].fitness)
    {
      cur_best = mem;
      population[POPSIZE].fitness = population[mem].fitness;
    }
  }
  // 
  //  Once the best member in the population is found, copy the genes.
  //
  population[POPSIZE].gene.resize(population[cur_best].gene.size());

  for (i = 0; i < population[cur_best].gene.size(); i++)
  {
    //previous_size = population[cur_best].gene.size();


    population[POPSIZE].gene[i]=population[cur_best].gene[i];
  }

  population[POPSIZE].size = population[POPSIZE].gene.size();

  return;
}

void keep_the_best_p2()

{
  int cur_best;
  int mem;
  int i;
  int previous_size;
  cur_best = 0;

  for (mem = 0; mem < POPSIZE; mem++)
  {
    if (population[POPSIZE].fitness > population[mem].fitness)
    {
      cur_best = mem;
      population[POPSIZE].fitness = population[mem].fitness;
    }
  }
  // 
  //  Once the best member in the population is found, copy the genes.
  //
  population[POPSIZE].gene.resize(population[cur_best].gene.size());

  for (i = 0; i < population[cur_best].gene.size(); i++)
  {
    //previous_size = population[cur_best].gene.size();


    population[POPSIZE].gene[i]=population[cur_best].gene[i];
  }

  population[POPSIZE].size = population[POPSIZE].gene.size();

  return;
}
//****************************************************************************80

vector<int> mutate_p1(int &seed, vector<vector<char> > &PIvector)
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  double lbound;
  double ubound;
  double x;
  double y;
  double z;
  int genpool_size = PIvector.size();

  // return the index of the individual that has been mutated
  vector<int> mutated_individual;
  int has_been_recorded = 0;

  for (i = 0; i < POPSIZE; i++)                                           
  {
    has_been_recorded = 0;

    int num_insert = 3;
    int current_pop_gen_size = population[i].gene.size();

    while(num_insert > 0){
    x = r8_uniform_ab(a, b, seed);                                    // randomly insert a gene into that location
    current_pop_gen_size = population[i].gene.size();
        if ((x < PMUTATION) && (genpool_size != 0) && (current_pop_gen_size != 0))
        {
          int genpool_loc = rand() % genpool_size;
          int gene_loc = rand() % current_pop_gen_size;
          population[i].gene.insert(population[i].gene.begin() + gene_loc, PIvector[genpool_loc]);
          if(has_been_recorded == 0){
            mutated_individual.push_back(i);
            has_been_recorded = 1;
          }
        }

      num_insert --;
      }

    current_pop_gen_size = population[i].gene.size();
    y = r8_uniform_ab(a, b, seed);                                 // randomly delete one of the gene from that indiv
        if ((y < PMUTATION) && (current_pop_gen_size > 2))
        {
          int gene_loc2 = rand() % current_pop_gen_size;
          population[i].gene.erase(population[i].gene.begin() + gene_loc2);

          if(has_been_recorded == 0){
            mutated_individual.push_back(i);
            has_been_recorded = 1;
          }
        }

    z = r8_uniform_ab(a, b, seed);                                 // randomly swap two genes in two individuals

    current_pop_gen_size = population[i].gene.size();
    int swap_pop = rand() % (POPSIZE - 1);
    int swap_size = population[swap_pop].gene.size();

    if ((z < PMUTATION) && (swap_size != 0) && (current_pop_gen_size))
    {
      
      int gene_loc3 = rand() % current_pop_gen_size;
      int gene_loc4 = rand() % population[swap_pop].gene.size();
      vector<char> temp;
      temp = population[i].gene[gene_loc3];
      population[i].gene[gene_loc3] = population[swap_pop].gene[gene_loc4];
      population[swap_pop].gene[gene_loc4] = temp;

      if(has_been_recorded == 0){
        mutated_individual.push_back(i);
        has_been_recorded = 1;
      }
    }

    population[i].size = population[i].gene.size();
  }


  return mutated_individual;
}


vector<int> mutate_p2(int &seed, vector<vector<char> > &PIvector)
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  double lbound;
  double ubound;
  double x;
  double y;
  double z;
  int genpool_size = PIvector.size();

  // return the index of the individual that has been mutated
  vector<int> mutated_individual;
  int has_been_recorded = 0;

  for (i = 0; i < POPSIZE; i++)
  {
    Oldpopulation[i] = population[i];
  }

  for (i = 0; i < POPSIZE; i++)                                           
  {
    has_been_recorded = 0;
    x = r8_uniform_ab(a, b, seed);                                    // randomly insert a gene into that location
    int current_pop_gen_size = population[i].gene.size();


        // if ((x < PMUTATION) && (genpool_size != 0) && (current_pop_gen_size != 0))
        // {
        //   int genpool_loc = rand() % genpool_size;
        //   int gene_loc = rand() % current_pop_gen_size;
        //   population[i].gene.insert(population[i].gene.begin() + gene_loc, PIvector[genpool_loc]);
        //   if(has_been_recorded == 0){
        //     mutated_individual.push_back(i);
        //     has_been_recorded = 1;
        //   }
        // }

                                                         // randomly delete one of the gene from that indiv
      // for(j=0; j< population[i].gene.size(); j++)
      // {
      //      y = r8_uniform_ab(a, b, seed);    
      //    if ((y < PMUTATION) && (current_pop_gen_size != 0))
      //    {
      //      population[i].gene.erase(population[i].gene.begin() + j);
      //   //   if(has_been_recorded == 0){
      //   //     mutated_individual.push_back(i);
      //   //     has_been_recorded = 1;
      //   //   }
      //    }
        
      //   }

    int num_delete = 3;
    while(num_delete > 0){
      current_pop_gen_size = population[i].gene.size();
      y = r8_uniform_ab(a, b, seed);                                 // randomly delete one of the gene from that indiv
        if ((y < PMUTATION_P2_DELET) && (current_pop_gen_size != 0) && (current_pop_gen_size > 2))
        {
          int gene_loc2 = rand() % current_pop_gen_size;
          population[i].gene.erase(population[i].gene.begin() + gene_loc2);
        }

    z = r8_uniform_ab(a, b, seed);                                 // randomly swap two genes of a individual
    num_delete --;
  }


    current_pop_gen_size = population[i].gene.size();
    if ((z < PMUTATION) && (current_pop_gen_size != 0) && (current_pop_gen_size != 1))
    {
      
      int gene_loc3 = rand() % current_pop_gen_size;
      int gene_loc4 = rand() % current_pop_gen_size;
      while(gene_loc3 == gene_loc4) gene_loc4 = rand() % current_pop_gen_size;
     // cout << i << " current_pop_gen_size1 = " << current_pop_gen_size <<endl;
      vector<char> temp;
      temp = population[i].gene[gene_loc3];
      population[i].gene[gene_loc3] = population[i].gene[gene_loc4];
      population[i].gene[gene_loc4] = temp;
    }

    //cout << i<<" current_pop_gen_size2 = " << current_pop_gen_size <<endl;
    population[i].size = population[i].gene.size();
  }


  return mutated_individual;
}
//****************************************************************************80

double r8_uniform_ab(double a, double b, int &seed)

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB returns a scaled pseudorandom R8.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_AB, a number strictly between A and B.
//
{
  int i4_huge = 2147483647;
  int k;
  double value;

  if (seed == 0)
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit(1);
  }

  k = seed / 127773;

  seed = 16807 * (seed - k * 127773) - k * 2836;

  if (seed < 0)
  {
    seed = seed + i4_huge;
  }

  value = (double)(seed) * 4.656612875E-10;

  value = a + (b - a) * value;

  return value;
}
//****************************************************************************80

void report(int generation, vector<vector<char> > &PIvector)
{
  double avg;
  double best_val;
  int i;
  int j;
  double square_sum;
  double stddev;
  double sum;
  double sum_square;
  double best_indiv_size;
  int PI_num = PIvector[0].size();


  if (generation == 0)
  {
    cout << "\n";
    cout << "  Generation       Best      Individual       Average       Standard              Best indivisual\n";
    cout << "  number           value        size          fitness       deviation                  is \n";
  }

  sum = 0.0;
  sum_square = 0.0;

  for (i = 0; i < POPSIZE; i++)
  {
    sum = sum + population[i].fitness;
    sum_square = sum_square + population[i].fitness * population[i].fitness;
  }

  avg = sum / (double)POPSIZE;
  square_sum = avg * avg * POPSIZE;
  stddev = sqrt((sum_square - square_sum) / (POPSIZE - 1));
  best_val = population[POPSIZE].fitness;
  population[POPSIZE].size = population[POPSIZE].gene.size();
  best_indiv_size = population[POPSIZE].size;

  cout << "  " << setw(4) << generation
    << "  " << setw(14) << best_val
    << "  " << setw(10) << best_indiv_size
    << "  " << setw(15) << avg
    << "  " << setw(15) << stddev;

  cout << "  " << setw(10);

  // for (i = 0; i < population[POPSIZE].gene.size(); i++)
  // {
  //   cout << "  " << setw(2);
  //   for (j = 0; j < PI_num; j++)
  //   {
  //     cout <<population[POPSIZE].gene[i][j];
  //   }
  // }
  cout << "\n";
  return;
}
//****************************************************************************80

int selector_p1(int &seed, int pop_size)
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  int mem;
  double p;
  double sum;

  int individual_to_be_killed = -1;

  double min_fault_coverage = 1.1;

  for(int i = 0; i < pop_size; i++){
    if(population[i].fault_coverage < min_fault_coverage){
      individual_to_be_killed = i;
      min_fault_coverage = population[i].fault_coverage;
    }
  }

  return individual_to_be_killed;

  //
  //  Find the total fitness of the population.
  //
  // sum = 0.0;
  // for (mem = 0; mem < POPSIZE; mem++){
  //   sum = sum + population[mem].fitness;
  // }
  // //
  // //  Calculate the relative fitness of each member.
  // //
  // for (mem = 0; mem < POPSIZE; mem++)
  // {
  //   population[mem].rfitness = population[mem].fitness / sum;
  // }
  // // 
  // //  Calculate the cumulative fitness.
  // //
  // population[0].cfitness = population[0].rfitness;
  // for (mem = 1; mem < POPSIZE; mem++)
  // {
  //   population[mem].cfitness = population[mem - 1].cfitness +
  //     population[mem].rfitness;
  // }
  // // 
  // //  Select survivors using cumulative fitness. 
  // //
  // for (i = 0; i < POPSIZE; i++)
  // {
  //   p = r8_uniform_ab(a, b, seed);
  //   if (p < population[0].cfitness)
  //   {
  //     newpopulation[i] = population[0];
  //   }
  //   else
  //   {
  //     for (j = 0; j < POPSIZE; j++)
  //     {
  //       if (population[j].cfitness <= p && p < population[j + 1].cfitness)
  //       {
  //         newpopulation[i] = population[j + 1];
  //       }
  //     }
  //   }
  // }
  // // 
  // //  Overwrite the old population with the new one.
  // //
  // for (i = 0; i < POPSIZE; i++)
  // {
  //   population[i] = newpopulation[i];
  // }


}

void selector_p2(int &seed, int pop_size)
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  int mem;
  double p;
  double sum;

  
   //Find the total fitness of the population.
  
  sum = 0.0;
  for (mem = 0; mem < POPSIZE; mem++){
    if(population[mem].fitness != 0) population[mem].d_fitness = (double)1/population[mem].fitness;
    else population[mem].d_fitness = 0.0;

    sum = sum + population[mem].d_fitness;
  }
  //
  //  Calculate the relative fitness of each member.
  //

  if(sum == 0) {
    cout << "sum == 0 ???????"<<endl;
    return;
  }
  for (mem = 0; mem < POPSIZE; mem++)
  {
    population[mem].rfitness = population[mem].d_fitness / sum;
  }
  // 
  //  Calculate the cumulative fitness.
  //
  population[0].cfitness = population[0].rfitness;
  for (mem = 1; mem < POPSIZE; mem++)
  {
    population[mem].cfitness = population[mem - 1].cfitness +
      population[mem].rfitness;
  }
  // 
  //  Select survivors using cumulative fitness. 
  //
  for (i = 0; i < POPSIZE; i++)
  {
    p = r8_uniform_ab(a, b, seed);
    if (p < population[0].cfitness)
    {
      newpopulation[i] = population[0];
    }
    else
    {
      for (j = 0; j < POPSIZE; j++)
      {
        if (population[j].cfitness <= p && p < population[j + 1].cfitness)
        {
          newpopulation[i] = population[j + 1];
        }
      }
    }
  }
  // 
  //  Overwrite the old population with the new one.
  //
  for (i = 0; i < POPSIZE; i++)
  {
    population[i] = newpopulation[i];
  }
}
//****************************************************************************80

void timestamp()

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    John Burkardt
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

void Xover_p1(int one, int two, int &seed, int individual_to_be_killed)

//****************************************************************************80
// 
//  Purpose:
//
//    XOVER performs crossover of the two selected parents. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int point, the crossover point.
//
//  Parameters:
//
//    Input, int ONE, TWO, the indices of the two parents.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  int i;
  int point;
  double t;
  int size_of_one;
  int size_of_two;
  int new_size;
  // 
  //  Select the crossover point.
  //
  point = i4_uniform_ab(0, NVARS - 1, seed);
  //
  //  Swap genes in positions 0 through POINT-1.
  //

/*
  for (i = 0; i < point; i++)
  {
    t = population[one].gene[i];
    population[one].gene[i] = population[two].gene[i];
    population[two].gene[i] = t;
  }
*/

  for (i = 0; i < POPSIZE; i++)
  {
    Oldpopulation[i] = population[i];
  }

  size_of_one = population[one].size;
  size_of_two = population[two].size;
  new_size = ceil((size_of_one + size_of_two) / 2);

  population[one].gene.resize(new_size);
  population[one].size = new_size;

  for (i = 0; i < population[one].size; i++)
  {
    if (i < size_of_one/2){
      
    }
    else{
      population[one].gene[i] = population[two].gene[i - floor(size_of_one / 2)];
    }
  }



  return;
}