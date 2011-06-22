/**
 * Interaction.h
 *
 * Data structure for an Interaction in the Cell.
 *
 */
#ifndef INTERACTION_H_
#define INTERACTION_H_

#include "Molecule.h"
#include "CustomMolecules.h"

#include "lemon/list_graph.h"
using namespace lemon; 


class Interaction{

public:
	Interaction();
	~Interaction();

//	virtual float getEffect(DerivGraph*, ListDigraph::Node, int, float);
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);

	const char* getName();

	float setRate(float);
	float getRate();
	const char* name;
	int arcID;
		
protected:
	float rate;
};



#endif
