/**
 * Interaction.h
 *
 * Data structure for an Interaction in the Cell.
 *
 */
#ifndef INTERACTION_H_
#define INTERACTION_H_

#include "Molecule.h"

#include "lemon/list_graph.h"
using namespace lemon; 



class Interaction{

public:
	Interaction();
	~Interaction();

	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);
	void setRate(float);
	
	float rkApprox(ListDigraph::NodeMap<Molecule*>*, ListDigraph::Node, float, int);
	int arcID;
		
protected:
	float rate;
};



#endif
