/**
 * Implementation file for Custom Interactions.
 *
 * Interactions may overload the virtual method getEffect() to create a custom effect between Molecules
 *
 */

#include "Interaction.h"
#include <cstdio>
#include "lemon/list_graph.h"

#include "CustomInteractions.h"

#include "ExternTrace.h"


Test::Test(){


	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type Test\n");
	t.trace("mloc","Interaction at location %u\n", (unsigned int)this);
	
	rate = 3.0;
	
	t.trace("init","New Interaction created\n");
}

Test::~Test(){


}


float Test::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a){
	t.trace("mloc","Interaction %u trace location: %u\n",(unsigned int)this, (unsigned int)&t);
	t.trace("efct","Original Node value: %f\n", (*m)[a]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue());

	if(g->source(g->arcFromId(arcID)) == a)
		return 0;
	else
		return (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue() * rate;

}

