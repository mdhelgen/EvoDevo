#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <iostream>

#include "Experiment.h"

#include "Trace.h"
//#include "ExternTrace.h"

using namespace std;

Trace t;

//global variable for hill parameter
//much easier than passing this through 4 class constructors to get to the DNA object
// consider moving all command line parameter values here
int hillParam = 1;

int main(int argc, char** argv){
  
  //flags set by command line options
  int verbose_flag = 0;
  int graphviz_flag = 0;
  int gnuplot_flag = 0;
  int outputall_flag = 0;
  int csvCell_flag = 0;
  int csvData_flag = 0;
  int usage_flag = 0;
  int gillespie_flag = 0;
  int rungeKutta_flag = 0;

  // used by command line parser
  int c;


  //trace related trace messages
  t.addTraceType("trce",0);
  
  //program arguments
  t.addTraceType("args",0);

  //error messages
  t.addTraceType("error",1);
 
  //object creation / construction
  t.addTraceType("init",0);

  //memory location of created objects
  t.addTraceType("mloc",0);

  //object deletion / destruction
  t.addTraceType("free",0);

  //calculated effect of interactions
  t.addTraceType("efct",0);

  //generational messages
  t.addTraceType("gens",0);

  // runge kutta (data / general messages)
  t.addTraceType("rk-4",1);
  
  // runge kutta (rk vals)
  t.addTraceType("rk-val",1);
 
  
  // runge kutta (calculation of next points) 
  t.addTraceType("rk-new",0);

  // hill / goodwin term calculation (transcription::getEffect)
  t.addTraceType("hill",0);

  // molecule scoring
  t.addTraceType("score",0);
  
  // mutation
  t.addTraceType("mutate",1);

  // polymorphic comparisons (not used?)
  t.addTraceType("typeid",0);

  t.addTraceType("stoch",1);

  int numCells = 2;
  int numGenerations = 10;

  int scoringInterval = 1;

  int maxBasic = 1;
  int maxPTM = 1;
  int maxComp = 1;
  int maxPromoter = 1;

  float minKineticRate = 0;
  float maxKineticRate = 1;

  float initialConcentration = 0;

  float rkTimeLimit = 20;
  float rkTimeStep = .05;


  // this loop parses the command line options. it was mostly adapted from online examples
  while(1){

  //define the command line options and their usage
    static struct option long_options[] =
     {
      {"help", no_argument, &usage_flag, 1},
      {"usage", no_argument, &usage_flag, 1},
      {"graphviz", no_argument, &graphviz_flag, 1},
      {"gnuplot",   no_argument, &gnuplot_flag, 1},
      {"outputall", no_argument, &outputall_flag, 1},
      {"csvCell", no_argument, &csvCell_flag, 1}, 
      {"csvData", no_argument, &csvData_flag, 1},
	  {"deterministic", no_argument, &rungeKutta_flag, 1},
	  {"stochastic", no_argument, &gillespie_flag, 1},

      {"cells",  required_argument, 0, 'c'},
      {"gens",  required_argument, 0, 'g'},
      {"minrate", required_argument, 0, 'a'},
      {"maxrate", required_argument, 0, 'b'},
      {"maxbasic", required_argument, 0, 'd'},
      {"maxptm", required_argument, 0, 'e'},
      {"maxcomp", required_argument, 0, 'f'},
      {"maxprom", required_argument, 0, 'h'},
      {"initconc", required_argument, 0, 'i'},
      {"rklim", required_argument, 0, 'j'},
      {"rkstep", required_argument, 0, 'k'},
      {"interval", required_argument, 0, 'l'},
      {"hill", required_argument, 0, 'm'},

      {0,0,0,0}
     };

 int option_index = 0;

 // which argument is currently seen?
 c = getopt_long (argc, argv, "a:b:c:d:e:f:g:h:i:j:k:l:", long_options, &option_index);

 if (c == -1)
 	break;

// what effect does that argument have
 switch(c)
 {
	case 0:
		if(long_options[option_index].flag != 0)
			break;
		printf("option %s", long_options[option_index].name);
		if (optarg)
			printf(" with arg %s", optarg);
		printf("\n");
		break;
	case 'a':
		minKineticRate = atof(optarg);
		break;
	case 'b':
		maxKineticRate = atof(optarg);
		break;
	case 'c':
		numCells = atoi(optarg);
		break;
	case 'd':
		maxBasic = atoi(optarg);
		break;
	case 'e':
		maxPTM = atoi(optarg);
		break;
	case 'f':
		maxComp = atoi(optarg);
		break;
	case 'g':
		numGenerations = atoi(optarg);
		break;
	case 'h':
		maxPromoter = atoi(optarg);
		break;
	case 'i':
		initialConcentration = atof(optarg);
		break;
	case 'j':
		rkTimeLimit = atof(optarg);
		break;
	case 'k':
		rkTimeStep = atof(optarg);
		break;
	case 'l':
		scoringInterval = atoi(optarg);
		break;
	case 'm':
		hillParam = atoi(optarg);
		break;
	case '?':
		break;
	default:
		abort();
}	
}

if (usage_flag){
      printf("EvoDevo Help\n");
      printf("---------------------------\n");
      printf("Flags:\n");
      printf("  --help            This help message\n");
      printf("  --usage           This help message\n");
      printf("  --graphviz        Output graphviz png files displaying the cell configuration\n");
      printf("  --gnuplot         Output gnuplot png files displaying molecule concentrations over time\n");
      printf("  --outputall       Output data about each cell every generation\n");
      printf("  --csvCell         Output csv data containing cell configuration\n");
      printf("  --csvData         Output csv data containing molecule concentrations\n");
	  printf("  --deterministic   Use deterministic Runge-Kutta solver for solving curves\n");
	  printf("  --stochastic      Use stochastic gillespie algorithm for solving curves\n");
      printf("\n");
      printf("Parameters:\n");
      printf("  --cells <int>        Number of Cells to simulate\n");
      printf("  --gens  <int>        Number of Generations to run for\n");
      printf("  --minrate <float>    Minimum value for random kinetic rates\n");
      printf("  --maxrate <float>    Maximum value for random kinetic rates\n");
      printf("  --maxbasic <int>     Maximum number of Basic Proteins\n");
      printf("  --maxptm <int>       Maximum number of Post Translationally Modified proteins\n");
      printf("  --maxcomp <int>      Maximum number of Protein-Protein Complexes\n");
      printf("  --maxprom <int>      Maximum number of Protein-Promoter Interactions\n");
      printf("  --initconc <float>   Initial concentration of molecules\n");
      printf("  --rklim <float>      Upper limit on time (x-axis) for differential equation solving\n");
      printf("  --rkstep <float>     Step size between points for differential equation solving\n");
      printf("  --interval <int>     Number of generations between equation solving and scoring\n");
      printf("  --hill <int>         Value of Hill Coefficient for DNA Transcription\n");

return 0;
}

// create our experiment with the options from the command line
Experiment e = Experiment(numCells, numGenerations, maxBasic, maxPTM, maxComp, maxPromoter, minKineticRate, maxKineticRate, rkTimeLimit, rkTimeStep, initialConcentration, rungeKutta_flag, gillespie_flag);

//set options related to output
e.setOutputOptions(graphviz_flag, gnuplot_flag, outputall_flag, csvCell_flag, csvData_flag, scoringInterval);

//start the experiment
e.start();


return 0;

}
