#ifndef DERIVGRAPH_H_
#define DERIVGRAPH_H_

#include "lemon/list_graph.h"
#include "lemon/concepts/maps.h"
#include <cstdio>
#include "Interaction.h"
#include "CustomInteractions.h"
#include "CustomMolecules.h"
#include <vector>

using namespace std;
using namespace lemon;

class DerivGraph{

public:
	DerivGraph();
	~DerivGraph();
	
	void test();
	void rungeKuttaEvaluate(float);
private:
	ListDigraph* derivs;
	ListDigraph::NodeMap<Molecule*>* molecules;
	ListDigraph::ArcMap<Interaction*>* interactions;
	
	float getEffect(ListDigraph::Node, ListDigraph::Arc, int, float);
	ListDigraph::Node add(Molecule*);
	ListDigraph::Arc add(Interaction*, ListDigraph::Node, ListDigraph::Node);
};





#endif
