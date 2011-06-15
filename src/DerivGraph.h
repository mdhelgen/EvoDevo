#ifndef DERIVGRAPH_H_
#define DERIVGRAPH_H_

#include "lemon/list_graph.h"
#include "lemon/concepts/maps.h"

using namespace lemon;

class DerivGraph{

public:
	DerivGraph();
	~DerivGraph();

private:
	ListDigraph* derivs;
	ListDigraph::NodeMap<int>* molecules;
	ListDigraph::ArcMap<int>* interactions;


};





#endif
