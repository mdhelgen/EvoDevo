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

   ListDigraph::Node n1 = derivs->addNode();
   (*molecules)[n1] = new Molecule();

   ListDigraph::Node n2 = derivs->addNode();
   (*molecules)[n2] = new DNA();

   ListDigraph::Arc a1 = derivs->addArc(n1, n2);
   (*interactions)[a1] = new Interaction();
   (*interactions)[a1]->arcID = derivs->id(a1);
   
   ListDigraph::Arc a2 = derivs->addArc(n1, n2);
   (*interactions)[a2] = new Test();
   (*interactions)[a2]->arcID = derivs->id(a2);


   ListDigraph::Arc a3 = derivs->addArc(n2, n1);
   (*interactions)[a3] = new Interaction();
   (*interactions)[a3]->arcID = derivs->id(a3);
   
   ListDigraph::Arc a4 = derivs->addArc(n2, n1);
   (*interactions)[a4] = new Test();
   (*interactions)[a4]->arcID = derivs->id(a4);


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
