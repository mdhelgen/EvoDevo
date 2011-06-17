/**
 * DerivGraph.cpp
 *
 * Interface to the LEMON graph library.
 *
 * The DerivGraph class should abstract the low level LEMON graph work.
 *   This includes:
 *     Adding molecules
 *     Mutation of the cell
 *     Runge-Kutta solving
 */

#include <iostream>
#include "DerivGraph.h"
#include "Interaction.h"
using namespace std;

#include "ExternTrace.h"

/**
 * DerivGraph::DerivGraph()
 *
 * DerivGraph constructor.
 *
 * The DerivGraph holds LEMON objects such as ListDigraph, NodeMap, and ArcMap.
 * It also holds the data produced by Runge-Kutta and facilitates plotting using Gnuplot.
 *
 * The derivatives describing the concentration of molecules in the cell can be represented as a directed graph.
 *
 * Each Node represents a type of molecule in the cell, and each Arc represents an interaction which has an
 * effect on the Nodes which it connects.
 *
 * Allocates:
 *     1 ListDigraph() object
 *     2 ListDigraph::NodeMap objects
 *     1 ListDigraph::ArcMap object
 */
DerivGraph::DerivGraph(){
    
    t.trace("init","Creating new DerivGraph\n");
    
    t.trace("mloc","DerivGraph location at %u\n",(unsigned int) this);

    //create the directed graph
    derivs = new ListDigraph();
    t.trace("mloc","DerivGraph %u ListDigraph location at %u\n", (unsigned int) this, (unsigned int) derivs);

    //map molecules onto the nodes
    molecules = new ListDigraph::NodeMap<Molecule*>(*derivs);
    t.trace("mloc","DerivGraph %u NodeMap location at %u\n", (unsigned int) this, (unsigned int) molecules);
    
    //map interactions onto the arcs
    interactions = new ListDigraph::ArcMap<Interaction*>(*derivs);
    t.trace("mloc","DerivGraph %u ArcMap location at %u\n", (unsigned int) this, (unsigned int) interactions);
   
    t.trace("init","New DerivGraph created\n");
}

/**
 * DerivGraph::~DerivGraph()
 *
 * DerivGraph Destructor.
 *
 * Frees:
 * 	1 NodeMap object
 * 	   n contained Molecule objects
 * 	1 NodeMap object
 * 	   n contained Vector objects
 * 	1 ArcMap object
 * 	   m contained Interaction objects
 * 	1 ListDigraph object
 * 	
 *
 */
DerivGraph::~DerivGraph(){
  
   //delete all Molecule objects mapped by Nodes
   t.trace("free","Deleting members of NodeMap at location %u\n",(unsigned int) molecules);
   for(ListDigraph::NodeIt it(*derivs); it !=INVALID; ++it){
	t.trace("free","Deleting NodeMap member at location %u\n",(unsigned int) (*molecules)[it]);
   	delete (*molecules)[it];
   }
  
   //delete the Molecule NodeMap
   t.trace("free","Deleting NodeMap object at location %u\n",molecules);
   delete molecules;

   //delete all Interaction objects mapped by Arcs
   t.trace("free","Deleting members of ArcMap at location %u\n", interactions);
   for(ListDigraph::ArcIt it(*derivs); it !=INVALID; ++it){
   	t.trace("free","Deleting ArcMap member at location %d\n", (*interactions)[it]);
	delete (*interactions)[it];
   }

   //delete the Interaction ArcMap
   t.trace("free","Deleting ArcMap object at location %d\n",interactions);
   delete interactions;


   //delete the ListDigraph
   t.trace("free","Deleting ListDigraph object at location %d\n",derivs);
   delete derivs;
}

void DerivGraph::test(){

	ListDigraph::Node A = add(new Molecule());
	ListDigraph::Node B = add(new Molecule());
	ListDigraph::Node C = add(new Molecule());

	(*molecules)[A]->setValue(2);
	(*molecules)[B]->setValue(4);
	(*molecules)[C]->setValue(1);

	ListDigraph::Arc AB = add(new Interaction, A, B);
	ListDigraph::Arc BC = add(new Interaction, B, C);
	ListDigraph::Arc AC = add(new Interaction, A, C);
	ListDigraph::Arc CB = add(new Interaction, C, B);
	
	(*interactions)[AB]->setRate(.01);
	(*interactions)[BC]->setRate(.07);
	(*interactions)[AC]->setRate(.03);
	(*interactions)[CB]->setRate(.09);


	float rkStep = 1.0;
  
	rungeKuttaEvaluate(rkStep); 
}

void DerivGraph::rungeKuttaEvaluate(float rkStep){

	for(int i = 0; i< 10; i++){
	for(int k = 0; k<4; k++){
	for(ListDigraph::ArcIt it(*derivs); it != INVALID; ++it){
		t.trace("rk-4","%f --> %f  (%f) s\n", (*molecules)[derivs->source(it)]->getValue(), 
						      (*molecules)[derivs->target(it)]->getValue(),
						      getEffect(derivs->source(it), it, 0, rkStep));

		t.trace("rk-4","%f --> %f  (%f) t\n", (*molecules)[derivs->source(it)]->getValue(),
						      (*molecules)[derivs->target(it)]->getValue(), 
						      getEffect(derivs->target(it), it, 0, rkStep)); 
		
			(*molecules)[derivs->source(it)]->updateRkVal(k, getEffect(derivs->source(it), it, k, rkStep));
			(*molecules)[derivs->target(it)]->updateRkVal(k, getEffect(derivs->target(it), it, k, rkStep));

	}
}
	for(ListDigraph::NodeIt it(*derivs); it != INVALID; ++it){
		(*molecules)[it]->nextPoint(rkStep);
	}
}
	for(ListDigraph::NodeIt it(*derivs); it != INVALID; ++it){
		(*molecules)[it]->outputRK();	
	}
}

/**
 * float DerivGraph::getEffect(ListDigraph::Node, ListDigraph::Arc)
 *
 * Get the effect that a particular interaction will have on another node.
 *
 * The effect of an Interaction is defined by the Interaction::getEffect method.
 * The result of this method will be a positive or negative change in concentration.
 *
 * @param m The Node containing the Molecule being affected
 * @param i The Arc containing the Interaction taking place
 * @param rkIteration The current iteration of the Runge-Kutta algorithm
 *
 * @return  the effect a particular Arc (Interaction) will have on a Node (Molecule)
 *
 */
float DerivGraph::getEffect(ListDigraph::Node m, ListDigraph::Arc i, int rkIteration, float rkStep){

	return (*interactions)[i]->getEffect(derivs, molecules, interactions, m, rkIteration, rkStep);

}

/**
 * ListDigraph::Node DerivGraph::add(Molecule*)
 *
 * Add a molecule to the ListDigraph and set up the NodeMaps for the new Molecule.
 *
 * @param newMolecule The new molecule being added.
 *
 * @return the Node object in the ListDigraph containing the new Molecule.
 */
ListDigraph::Node DerivGraph::add(Molecule * newMolecule){
	
	//add a new node to the ListDigraph
	ListDigraph::Node newNode = derivs->addNode();
	
	//map the Molecule to the molecule NodeMap
	(*molecules)[newNode] = newMolecule;

	//return the newly created Node
	return newNode;
}

/**
 * ListDigraph::Arc DerivGraph::add(Interaction*, ListDigraph::Node, ListDigraph::Node)
 *
 * Add a new Interaction to the ListDigraph, between the two supplied Nodes
 * 
 * @param newInteraction The new Interaction being added
 * @param from The source node for the new Interaction
 * @param to The target Node for the new Interaction
 *
 * @return the Arc object in the ListDigraph containing the new Interaction.
 */
ListDigraph::Arc DerivGraph::add(Interaction * newInteraction, ListDigraph::Node from, ListDigraph::Node to){

	ListDigraph::Arc newArc = derivs->addArc(from, to);
	(*interactions)[newArc] = newInteraction;
	(*interactions)[newArc]->arcID = derivs->id(newArc);
	return newArc;
}
