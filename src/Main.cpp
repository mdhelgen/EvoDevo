
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <iostream>

#include "Experiment.h"

#include "Trace.h"
//#include "ExternTrace.h"

using namespace std;

Trace t;

int main(int argc, char** argv){

  int verbose_flag = 0;
  int graphviz_flag = 0;
  int gnuplot_flag = 0;
  int outputall_flag = 0;

  int c;


  //trace related trace messages
  t.addTraceType("trce",0);
  
  //program arguments
  t.addTraceType("args",0);
 
  //object creation / construction
  t.addTraceType("init",0);

  //memory location of created objects
  t.addTraceType("mloc",0);

  //object deletion / destruction
  t.addTraceType("free",0);

  //calculated effect of interactions
  t.addTraceType("efct",0);

  //generational messages
  t.addTraceType("gens",1);

  //runge kutta
  t.addTraceType("rk-4",0);

  t.addTraceType("rk-val",0);
  t.addTraceType("rk-new",0);

  t.addTraceType("hill",0);

  t.addTraceType("score",1);
  //mutation
  t.addTraceType("mutate",0);

  //polymorphic comparisons
  t.addTraceType("typeid",0);


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


  while(1){
    static struct option long_options[] =
     {
      {"graphviz", no_argument, &graphviz_flag, 1},
      {"gnuplot",   no_argument, &gnuplot_flag, 1},
      {"outputall", no_argument, &outputall_flag, 1},

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

      {0,0,0,0}
     };

 int option_index = 0;


 c = getopt_long (argc, argv, "a:b:c:d:e:f:g:h:i:j:k:l:", long_options, &option_index);

 if (c == -1)
 	break;
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
	case '?':
		break;
	default:
		abort();
}	
}

if(verbose_flag)
	printf("verbose flag set\n");
	



Experiment e = Experiment(numCells, numGenerations, maxBasic, maxPTM, maxComp, maxPromoter, minKineticRate, maxKineticRate, rkTimeLimit, rkTimeStep, initialConcentration);
e.setOutputOptions(graphviz_flag, gnuplot_flag, outputall_flag, scoringInterval);
e.start();


//t.trace("free","Deleting Experiment object at location %d\n", e);

return 0;

}
