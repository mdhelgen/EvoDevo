/**
 * Implementation file for Custom Interactions.
 *
 * Interactions may overload the virtual method getEffect() to create a custom effect between Molecules
 *
 */

#include "CustomInteractions.h"

#include "ExternTrace.h"

Transcription::Transcription(){

}

Transcription::~Transcription(){

}

//float Transcription::getEffect(DerivGraph* d, ListDigraph::Node a, int rkIter, float rkStep){


float Transcription::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a, int rkIter, float rkStep){	


	t.trace("efct","Original Node value: %f\n", (*m)[a]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue());

	Molecule* thisMol = (*m)[a];
	Molecule* oppositeMol = (*m)[g->oppositeNode(a, g->arcFromId(arcID))];

	if(g->source(g->arcFromId(arcID)) == a)
		return 0;
	else
		return oppositeMol->rkApprox(rkIter, rkStep) * rate;

}

Degradation::Degradation(){}
Degradation::~Degradation(){}

//float Degradation::getEffect(DerivGraph* d, ListDigraph::Node a, int rkIter, float rkStep){

float Degradation::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a, int rkIter, float rkStep){	
//	ListDigraph* g = d->getListDigraph();
//	ListDigraph::NodeMap<Molecule*>* m = d->getNodeMap();
	t.trace("efct","Original Node value: %f\n", (*m)[a]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue());

	Molecule* thisMol = (*m)[a];
	Molecule* oppositeMol = (*m)[g->oppositeNode(a, g->arcFromId(arcID))];

	if(g->source(g->arcFromId(arcID)) == a)
		return -1 * oppositeMol->rkApprox(rkIter, rkStep) * rate;
	else
		return 0;

}
Test::Test(){

	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type Test\n");
	t.trace("mloc","Interaction at location %u\n", (unsigned int)this);
	
	rate = 3.0;
	
	t.trace("init","New Interaction created\n");
}

Test::~Test(){


}


//float Test::getEffect(DerivGraph* d, ListDigraph::Node a, int rkIter, float rkStep){
float Test::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a, int rkIter, float rkStep){	
	
//	ListDigraph* g = d->getListDigraph();
//	ListDigraph::NodeMap<Molecule*>* m = d->getNodeMap();
//	ListDigraph::ArcMap<Interaction*>* i = d->getArcMap();
	
	t.trace("mloc","Interaction %u trace location: %u\n",(unsigned int)this, (unsigned int)&t);
	t.trace("efct","Original Node value: %f\n", (*m)[a]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue());


	Molecule* oppositeMol = (*m)[g->oppositeNode(a, g->arcFromId(arcID))];

	if(g->source(g->arcFromId(arcID)) == a)
		return 0;
	else
		return oppositeMol->rkApprox(rkIter, rkStep) * rate;


}

