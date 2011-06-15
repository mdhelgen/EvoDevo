
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <iostream>

#include "Experiment.h"

using namespace std;

int main(int argc, char** argv){

static int verbose_flag;
int c;

int numCells = 2;
int numGenerations = 10;


while(1){
static struct option long_options[] =
 {
  {"verbose", no_argument, &verbose_flag, 1},
  {"brief",   no_argument, &verbose_flag, 0},

  {"cells",  required_argument, 0, 'c'},
  {"gens",  required_argument, 0, 'g'},
  {0,0,0,0}
 };

 int option_index = 0;


 c = getopt_long (argc, argv, "abc:g:", long_options, &option_index);

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

	case 'c':
		numCells = atoi(optarg);
		break;
	case 'g':
		numGenerations = atoi(optarg);
		break;
	case '?':
		break;
	default:
		abort();
}	
}

if(verbose_flag)
	printf("verbose flag set\n");
	



Experiment* e = new Experiment(numCells, numGenerations);


return 0;
}
