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


float Interaction::getEffect(ListDigraph* g, ListDigraph::NodeMap<int>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a){
	t.trace("mloc","Interaction %u trace location: %u\n",(unsigned int)this, (unsigned int)&t);
	t.trace("efct","Original Node value: %d\n", (*m)[a]);
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %d\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]);
	

	return rate;

}

void Interaction::setRate(float f){

	rate = f;

}
