/**
 * InteractionType.h
 *
 * InteractionType holds default information about a type of interaction.
 */

#ifndef INTERACTION_TYPE_H_
#define INTERACTION_TYPE_H_

#include "lemon/list_graph.h"


class InteractionType{

public:
	InteractionType();
	~InteractionType();

	virtual float getEffect(ListDigraph::Arc a);

private:
	
	
	//Name associated with the type of interaction
	const char* name;

	//Short name associated with the type of interaction
	const char* shortName;

	//Default factor applied to incoming interaction effects
	float defaultProductFactor;

	//Default factor applied to outgoing interactino effects
	float defaultSubstrateFactor;

	//Default minimum rate for this interaction type
	float defaultMinRate;

	//Default maximum rate for this interaction type
	float defaultMaxRate;

	//Display this type of interaction in the output files
	bool showInOutput;

	//Calculate value updates during integration for this interaction type
	bool calcDuringRK;



}





#endif
