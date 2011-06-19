#ifndef DERIVGRAPH_H_
#define DERIVGRAPH_H_

#include "lemon/list_graph.h"
#include "lemon/concepts/maps.h"
#include <cstdio>
#include <vector>
//#include "Molecule.h"
//#include "Interaction.h"
//#include "CustomMolecules.h"
//#include "CustomInteractions.h"


using namespace std;
using namespace lemon;

class Molecule;
class Interaction;
//class Test;

class DerivGraph{

public:
	DerivGraph();
	~DerivGraph();
	
	void test();
	void rungeKuttaEvaluate(float);
	ListDigraph* getListDigraph();
	ListDigraph::NodeMap<Molecule*>* getNodeMap();
	ListDigraph::ArcMap<Interaction*>* getArcMap();

private:
	ListDigraph* derivs;
	ListDigraph::NodeMap<Molecule*>* molecules;
	ListDigraph::ArcMap<Interaction*>* interactions;
	
	float getEffect(ListDigraph::Node, ListDigraph::Arc, int, float);
	ListDigraph::Node add(Molecule*);
	ListDigraph::Arc add(Interaction*, ListDigraph::Node, ListDigraph::Node);
};





#endif
