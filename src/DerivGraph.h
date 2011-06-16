#ifndef DERIVGRAPH_H_
#define DERIVGRAPH_H_

#include "lemon/list_graph.h"
#include "lemon/concepts/maps.h"
#include <cstdio>
#include "Interaction.h"

using namespace lemon;

class DerivGraph{

public:
	DerivGraph();
	~DerivGraph();
	
	void test();
	float getEffect(ListDigraph::Node, ListDigraph::Arc);
private:
	ListDigraph* derivs;
	ListDigraph::NodeMap<int>* molecules;
	ListDigraph::ArcMap<Interaction*>* interactions;
};





#endif
