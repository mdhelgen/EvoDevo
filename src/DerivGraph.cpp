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
    molecules = new ListDigraph::NodeMap<int>(*derivs);
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
   
   t.trace("free","Deleting NodeMap object at location %d\n",molecules);
   delete molecules;


   t.trace("free","Deleting ArcMap object at location %d\n",interactions);
   delete interactions;


   t.trace("free","Deleting ListDigraph object at location %d\n",derivs);
   delete derivs;
}

void DerivGraph::test(){

   ListDigraph::Node n1 = derivs->addNode();
   (*molecules)[n1] = 2;

   ListDigraph::Node n2 = derivs->addNode();
   (*molecules)[n2] = 3;

   ListDigraph::Arc a1 = derivs->addArc(n1, n2);
   (*interactions)[a1] = new Interaction();
   (*interactions)[a1]->arcID = derivs->id(a1);

   t.trace("efct","DerivGraph getEffect: %f\n", getEffect(n1, a1));
   (*interactions)[a1]->setRate(8);
   t.trace("efct","DerivGraph getEffect2: %f\n", getEffect(n2, a1));
}

float DerivGraph::getEffect(ListDigraph::Node m, ListDigraph::Arc i){


	return (*interactions)[i]->getEffect(derivs, molecules, interactions, m);

}
