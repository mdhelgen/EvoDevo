/**
 * Interaction base class implementation
 *
 */

#include "Interaction.h"
#include <cstdio>
#include "lemon/list_graph.h"

#include "ExternTrace.h"
/**
 * Interaction::Interaction()
 * 
 * Interaction default constructor.
 */
Interaction::Interaction(){

	t.trace("init","Creating new Interaction\n");
	t.trace("mloc","Interaction at location %d\n", this);
	
	rate = 5.0;
	
	t.trace("init","New Interaction created\n");
}

/**
 * Interaction::~Interaction()
 *
 * Interaction default destructor.
 */
Interaction::~Interaction(){

	t.trace("free","Deleting Interaction at location %d\n", this);
}


float Interaction::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a, int rkIter, float rkStep){
	
	
	t.trace("efct","Original Node value: %f\n", (*m)[a]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue());
	if(g->source(g->arcFromId(arcID)) == a)
		return -1 * rkApprox(m, a, rkStep, rkIter)  * rate;
	else
		return rkApprox(m, g->oppositeNode(a, g->arcFromId(arcID)), rkStep, rkIter) * rate;
}


float Interaction::rkApprox(ListDigraph::NodeMap<Molecule*>* m, ListDigraph::Node a , float stepSize, int k){

	float concentration = (*m)[a]->getValue();
	float rkVal;	
	float kstep;
	switch(k){
	case 0:
		kstep = 1;
	 	rkVal = 0;	
		break;
	case 1:
	case 2:
		kstep = stepSize/2;
		rkVal = (*m)[a]->getrkVal(k-1);
		break;
	case 3:
		kstep = stepSize;
		rkVal = (*m)[a]->getrkVal(k-1);
		break;
	}

	return concentration + (rkVal * kstep); 

}

void Interaction::setRate(float f){

	rate = f;

}
