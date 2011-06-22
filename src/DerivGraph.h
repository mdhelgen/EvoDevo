#ifndef DERIVGRAPH_H_
#define DERIVGRAPH_H_

#include "lemon/list_graph.h"
#include "lemon/concepts/maps.h"
#include <cstdio>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <typeinfo>


#include "Molecule.h"
#include "Interaction.h"

#include "CustomInteractions.h"
#include "CustomMolecules.h"

using namespace std;
using namespace lemon;

class DerivGraph{

public:
	DerivGraph();
	~DerivGraph();
	
	void test();
	void rungeKuttaEvaluate(float);
	ListDigraph* getListDigraph();
	ListDigraph::NodeMap<Molecule*>* getNodeMap();
	ListDigraph::ArcMap<Interaction*>* getArcMap();
	void outputDotImage(int, int);
	
	//mutation methods
	void newBasic();


private:
	
	//graph structure
	ListDigraph* derivs;
	ListDigraph::NodeMap<Molecule*>* molecules;
	ListDigraph::ArcMap<Interaction*>* interactions;

	vector<Molecule*>* MoleculeList;
	vector<Protein*>* ProteinList;
	vector<mRNA*>* mRNAList;
	vector<DNA*>* DNAList;
	vector<Complex*>* ComplexList;


	//null node
	ListDigraph::Node nullnode;


	//helper method
	float getEffect(ListDigraph::Node, ListDigraph::Arc, int, float);
	//utility method
	ListDigraph::Node add(Molecule*);
	ListDigraph::Arc add(Interaction*, ListDigraph::Node, ListDigraph::Node);
	int count;
};


#endif
