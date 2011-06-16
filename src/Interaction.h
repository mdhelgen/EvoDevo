/**
 * Interaction.h
 *
 * Data structure for an Interaction in the Cell.
 *
 */
#include "lemon/list_graph.h"
using namespace lemon; 

#ifndef INTERACTION_H_
#define INTERACTION_H_

class Interaction{

public:
	Interaction();
	~Interaction();

	virtual float getEffect(ListDigraph* d, ListDigraph::NodeMap<int>* b, ListDigraph::ArcMap<Interaction*>* c, ListDigraph::Node a);
	void setRate(float);
	int arcID;
private:
	float rate;
};



#endif
