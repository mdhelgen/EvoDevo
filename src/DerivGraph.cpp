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
   
   t.trace("free","Deleting NodeMap object at location %u\n",molecules);
   delete molecules;

   //delete all Interaction objects mapped by Arcs
   t.trace("free","Deleting members of ArcMap at location %u\n", interactions);
   for(ListDigraph::ArcIt it(*derivs); it !=INVALID; ++it){
   	t.trace("free","Deleting ArcMap member at location %d\n", (*interactions)[it]);
	delete (*interactions)[it];
   }

   t.trace("free","Deleting ArcMap object at location %d\n",interactions);
   delete interactions;


   t.trace("free","Deleting ListDigraph object at location %d\n",derivs);
   delete derivs;
}

void DerivGraph::test(){

   ListDigraph::Node n1 = add(new Molecule());

   ListDigraph::Node n2 = add(new DNA());

   ListDigraph::Arc a1 = add(new Interaction(), n1, n2);
   
   ListDigraph::Arc a2 = add(new Test(), n1, n2);
   
   ListDigraph::Arc a3 = add(new Interaction(), n2, n1);
   
   ListDigraph::Arc a4 = add(new Test(), n2, n1);


   t.trace("efct","DerivGraph node 1 arc 1 getEffect: %f\n", getEffect(n1, a1));
   t.trace("efct","DerivGraph node 2 arc 1 getEffect: %f\n", getEffect(n2, a1));
   t.trace("efct","DerivGraph node 1 arc 2 getEffect: %f\n", getEffect(n1, a2));
   t.trace("efct","DerivGraph node 2 arc 2 getEffect: %f\n", getEffect(n2, a2));
   
   t.trace("efct","DerivGraph node 1 arc 3 getEffect: %f\n", getEffect(n1, a3));
   t.trace("efct","DerivGraph node 2 arc 3 getEffect: %f\n", getEffect(n2, a3));
   t.trace("efct","DerivGraph node 1 arc 4 getEffect: %f\n", getEffect(n1, a4));
   t.trace("efct","DerivGraph node 2 arc 4 getEffect: %f\n", getEffect(n2, a4));
   


}
float DerivGraph::getEffect(ListDigraph::Node m, ListDigraph::Arc i){


	return (*interactions)[i]->getEffect(derivs, molecules, interactions, m);

}

ListDigraph::Node DerivGraph::add(Molecule * newMolecule){
	
	ListDigraph::Node newNode = derivs->addNode();
	(*molecules)[newNode] = newMolecule;
	return newNode;
}

ListDigraph::Arc DerivGraph::add(Interaction * newInteraction, ListDigraph::Node from, ListDigraph::Node to){

	ListDigraph::Arc newArc = derivs->addArc(from, to);
	(*interactions)[newArc] = newInteraction;
	(*interactions)[newArc]->arcID = derivs->id(newArc);
	return newArc;
}