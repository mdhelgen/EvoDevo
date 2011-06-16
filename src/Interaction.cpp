#include "Interaction.h"
#include <cstdio>
#include "lemon/list_graph.h"

#include "ExternTrace.h"

Interaction::Interaction(){

	t.trace("init","Creating new Interaction\n");
	t.trace("mloc","Interaction at location %u\n", (unsigned int)this);
	
	
	rate = 5.0;
	
	t.trace("init","New Interaction created\n");
}

Interaction::~Interaction(){

}


float Interaction::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a){
	
	
	t.trace("efct","Original Node value: %f\n", (*m)[a]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue());

	if(g->source(g->arcFromId(arcID)) == a)
		return -1 * (*m)[a]->getValue() * rate;
	else
		return (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue() * rate;

}

void Interaction::setRate(float f){

	rate = f;

}
