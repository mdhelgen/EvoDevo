#ifndef DERIVGRAPH_H_
#define DERIVGRAPH_H_

#include "lemon/list_graph.h"
#include "lemon/concepts/maps.h"
#include <cstdio>
#include "Interaction.h"
#include "CustomInteractions.h"
#include "CustomMolecules.h"

using namespace lemon;

class DerivGraph{

public:
	DerivGraph();
	~DerivGraph();
	
	void test();
private:
	ListDigraph* derivs;
	ListDigraph::NodeMap<Molecule*>* molecules;
	ListDigraph::ArcMap<Interaction*>* interactions;
	
	//ListDigraph::NodeMap<vector<float>*>* rungeKuttaSolution;
        
	float getEffect(ListDigraph::Node, ListDigraph::Arc);
	ListDigraph::Node add(Molecule*);
	ListDigraph::Arc add(Interaction*, ListDigraph::Node, ListDigraph::Node);
};





#endif