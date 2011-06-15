#include <iostream>
#include "DerivGraph.h"

using namespace std;

#include "ExternTrace.h"

DerivGraph::DerivGraph(){


    derivs = new ListDigraph();
    molecules = new ListDigraph::NodeMap<int>(*derivs);
    interactions = new ListDigraph::ArcMap<int>(*derivs);
    t->trace("init","New DerivGraph Created\n");
    //molecules = new ListDigraph::NodeMap<int>(derivs);	
    //interactions = new ListDigraph::ArcMap<int>(derivs);
//derivs = new ListDigraph();
	//molecules = new ListDigraph::NodeMap<int>(*derivs);
}


DerivGraph::~DerivGraph(){

   t->trace("free","Deleting NodeMap object at location %d\n",molecules);
   delete molecules;
   t->trace("free","Deleting ArcMap object at location %d\n",interactions);
   delete interactions;
   t->trace("free","Deleting ListDigraph object at location %d\n",derivs);
   delete derivs;
    //delete molecules;
}
