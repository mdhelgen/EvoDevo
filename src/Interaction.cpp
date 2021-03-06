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
        name = "default";	
	rate = .05;
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

/**
 * float Interaction::getEffect(ListDigraph* , NodeMap<Molecule*>* , ArcMap<Interaction*>* , Node , int, float)
 *
 * Get the effect this interaction has on a particular node.
 *
 * This method defines the behavior of an interaction which connects two molecules. The effect on Node a can be dependent on
 * any other molecule, which can be accessed using the ListDigraph, NodeMap, and ArcMap parameters.
 *
 * Runge-Kutta iteratively approximates the change in concentration during a given timestep. The first iteration is based soley on the current concentration,
 * and each further iteration takes the result of the previous iteration into account. The Runge-Kutta data are stored in each molecule, and it is necessary to
 * call Molecule::rkApprox(stepsize, iteration) rather than Molecule::getValue() to get the current Iteration's approximated concentration. 
 *
 * @param g The graph object containing Node-Node relationships. 
 * @param m The NodeMap object containing Node-Molecule mappings.
 * @param i The ArcMap object containing Arc-Interaction mappings.
 * @param n The Node to calculate the effect for
 * @param rkIter The current iteration of Runge-Kutta [0,3]
 * @param rkStep The stepsize of Runge-Kutta
 *
 */
float Interaction::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node n, int rkIter, float rkStep){	
	
	t.trace("efct","Original Node value: %f\n", (*m)[n]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == n) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(n, g->arcFromId(arcID))]->getValue());
	Molecule* thisMol = (*m)[n];
	Molecule* oppositeMol = (*m)[g->oppositeNode(n, g->arcFromId(arcID))];
	

	if(isSourceNode(g, n) == 1)
		
		return -1 * thisMol->rkApprox(rkIter, rkStep) * rate;
	
	else if(isTargetNode(g, n) == 1)
		
		return oppositeMol->rkApprox(rkIter, rkStep) * rate;
	
	else{
		t.trace("error","%s reached error case, not source or target (%p)\n", name, this); 	
		return 0;
	}

}
/*
 * int Interaction::isSourceNode(ListDigraph*, ListDigraph::Node)
 *
 * Helper method to determine whether a particular node is the source or target node in a
 * particular arc.
 *
 * @param g
 * @param n
 *
 * @return 1 if the node is the source, and 0 otherwise.
 */
int Interaction::isSourceNode(ListDigraph* g, ListDigraph::Node n){

	if(g->source(g->arcFromId(arcID)) == n)
		return 1;
	else
		return 0;

}


/*
 * int Interaction::isTargetNode(ListDigraph*, ListDigraph::Node, ListDigraph::Arc)
 *
 * Helper method to determine whether a particular node is the target or source node in a
 * particular arc.
 *
 * @param g
 * @param n
 * @param a
 *
 * @return 1 if the node is the target, -1 if it is the source, and 0 otherwise.
 */
int Interaction::isTargetNode(ListDigraph* g, ListDigraph::Node n){

	if(g->target(g->arcFromId(arcID)) == n)
		return 1;
	else
		return 0;

}


/**
 * void Interaction::setRate(float)
 *
 * Change the kinetic rate of the Interaction
 *
 * @param f the new rate for the interaction
 * 
 * @return the old rate for the interaction
 */
float Interaction::setRate(float f){
	
	float oldRate = f;
	rate = f;
	return oldRate;

}

float Interaction::getRate(){
	return rate;
}

const char* Interaction::getName(){

	return name;

}
