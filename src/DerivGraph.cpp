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
 *     1 ListDigraph::NodeMap objects
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

    count=0;
}

/**
 * DerivGraph::~DerivGraph()
 *
 * DerivGraph Destructor.
 *
 * Frees:
 * 	1 NodeMap object
 * 	   n contained Molecule objects
 * 	1 ArcMap object
 * 	   m contained Interaction objects
 * 	1 ListDigraph object
 */
DerivGraph::~DerivGraph(){

   //delete all Molecule objects mapped by Nodes
   t.trace("free","Deleting members of NodeMap at location %u\n",(unsigned int) molecules);
   for(ListDigraph::NodeIt it(*derivs); it !=INVALID; ++it){
	t.trace("free","Deleting NodeMap member at location %u\n",(unsigned int) (*molecules)[it]);
	t.trace("free","longnname: %s\n", (*molecules)[it]->getLongName());
	t.trace("free","shortname: %s\n", (*molecules)[it]->getShortName());
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

/*
	ListDigraph::Node A = add(new Molecule());
	ListDigraph::Node B = add(new Molecule());
	ListDigraph::Node C = add(new Molecule());

	(*molecules)[A]->setValue(2);
	(*molecules)[B]->setValue(4);
	(*molecules)[C]->setValue(1);

	ListDigraph::Arc AB = add(new Interaction(), A, B);
	ListDigraph::Arc AC = add(new Interaction(), A, C);
	ListDigraph::Arc BC = add(new Interaction(), B, C);
	ListDigraph::Arc CB = add(new Interaction(), C, B);

	(*interactions)[AB]->setRate(.01);
	(*interactions)[AC]->setRate(.03);
	(*interactions)[BC]->setRate(.07);
	(*interactions)[CB]->setRate(.09);
*/
	ListDigraph::Node A = add(new DNA());
	ListDigraph::Node B = add(new mRNA());
	ListDigraph::Node C = add(new Protein());
	ListDigraph::Node D = add(new Null());

	ListDigraph::Node E = add(new DNA());
	ListDigraph::Node F = add(new mRNA());
	ListDigraph::Node G = add(new Protein());

	(*molecules)[A]->setValue(1);
	(*molecules)[B]->setValue(10);
	(*molecules)[C]->setValue(15);

	(*molecules)[E]->setValue(1);
	(*molecules)[F]->setValue(10);
	(*molecules)[G]->setValue(15);



	ListDigraph::Arc AB = add(new Transcription(), A, B);
	ListDigraph::Arc BC = add(new Translation(), B, C);
	ListDigraph::Arc BD = add(new Degradation(), B, D);
	ListDigraph::Arc CD = add(new Degradation(), C, D);
	
	ListDigraph::Arc EF = add(new Transcription(), E, F);
	ListDigraph::Arc FG = add(new Translation(), F, G);
	ListDigraph::Arc FD = add(new Degradation(), F, D);
	ListDigraph::Arc GD = add(new Degradation(), G, D);
	
	(*interactions)[AB]->setRate(.05);
	(*interactions)[BC]->setRate(.05);
	(*interactions)[BD]->setRate(.01);
	(*interactions)[CD]->setRate(.02);

	(*interactions)[EF]->setRate(.05);
	(*interactions)[FG]->setRate(.05);
	(*interactions)[FD]->setRate(.01);
	(*interactions)[GD]->setRate(.02);


	ListDigraph::Node CG = add(new Complex(derivs->id(C), derivs->id(G)));
	(*molecules)[CG]->setValue(5);

	ListDigraph::Arc CGf1 = add(new ForwardComplexation(derivs->id(C), derivs->id(G)), C, CG);
	ListDigraph::Arc CGf2 = add(new ForwardComplexation(derivs->id(C), derivs->id(G)), G, CG);
	ListDigraph::Arc CGr1 = add(new ReverseComplexation(derivs->id(C), derivs->id(G)), CG, C);
	ListDigraph::Arc CGr2 = add(new ReverseComplexation(derivs->id(C), derivs->id(G)), CG, G);
	ListDigraph::Arc CGdeg = add(new Degradation(), CG, D);

	(*interactions)[CGf1]->setRate(.1);
	(*interactions)[CGf2]->setRate(.1);

	(*interactions)[CGr1]->setRate(.08);
	(*interactions)[CGr2]->setRate(.08);

	(*interactions)[CGdeg]->setRate(.04);


	float rkStep = 1.0;
  
	rungeKuttaEvaluate(rkStep); 
}

/**
 * Uses the Runge-Kutta fourth order method to approximate the solutions to the system of differential equations
 *
 * The result of this algorithm is the vector rungeKuttaSolution within each Molecule object containing the approximation of
 * the concentration at each timestep.
 *
 * @param rkStep the timestep (precision) between calculated points
 */
void DerivGraph::rungeKuttaEvaluate(float rkStep){

	for(ListDigraph::NodeIt it(*derivs); it != INVALID; ++it){
		(*molecules)[it]->reset();
	}
	
	//time loop
	for(int i = 0; i< 10; i++){
			
		//this loop runs each of the four runge-kutta iterations
		for(int k = 0; k<4; k++){
			
			//loop through arc in the graph
			for(ListDigraph::ArcIt it(*derivs); it != INVALID; ++it){
				
				//get the effect the interaction has on the source node, and update the rkval
				(*molecules)[derivs->source(it)]->updateRkVal(k, getEffect(derivs->source(it), it, k, rkStep));
				//get the effect the interaction has on the target node, and update the rkval
				(*molecules)[derivs->target(it)]->updateRkVal(k, getEffect(derivs->target(it), it, k, rkStep));

			}
		}

		//after the four rkVals are calcualted for all molecules, the next point can be calculated
		for(ListDigraph::NodeIt it(*derivs); it != INVALID; ++it){
			(*molecules)[it]->nextPoint(rkStep);
		}
	}
	
	//test output, display the values calculated by runge kutta
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

	//get the Interaction for this arc
	Interaction* a = (*interactions)[i];

	//calculate the effect this Interaction will have on the Molecule in Node m
	
	//return a->getEffect(this, m, rkIteration, rkStep);
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

	//add a new Node to the graph
	ListDigraph::Node newNode = derivs->addNode();
	
	//map the new Node to the Molecule 
	(*molecules)[newNode] = newMolecule;

	(*molecules)[newNode]->nodeID = derivs->id(newNode);
	
	newMolecule->setID(count++);

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

	//add a new arc to the graph
	ListDigraph::Arc newArc = derivs->addArc(from, to);
	
	//map the new Arc to the Interaction 
	(*interactions)[newArc] = newInteraction;
	
	//store the arcID of the mapped Arc in the Interaction
	//this gives the Interaction information about its position in the graph
	(*interactions)[newArc]->arcID = derivs->id(newArc);

	//return the newly created Arc
	return newArc;
}

ListDigraph* DerivGraph::getListDigraph(){
	return derivs;
}

ListDigraph::NodeMap<Molecule*>* DerivGraph::getNodeMap(){
	return molecules;
}

ListDigraph::ArcMap<Interaction*>* DerivGraph::getArcMap(){
	return interactions;
}

void DerivGraph::outputDotImage(int cellNum, int gen){
if(0)
return;
	char buf[200];
	sprintf(buf, "circo -Tpng -oCell%dGen%d.png",cellNum, gen);
	FILE* dot = popen(buf,"w");
	
	fprintf(dot,"digraph mol_interactions {\n");
	fflush(dot);

//	fprintf(dot,"rankdir = LR;\n");
//	fflush(dot);

	fprintf(dot,"size=\"8,5\"\n");
	fflush(dot);

	fprintf(dot,"node [shape = circle];\n");
	fflush(dot);

	for(ListDigraph::ArcIt it(*derivs); it != INVALID; ++it){
		fprintf(dot, "%s -> %s [ label = \"%s\" ];\n",(*molecules)[derivs->source(it)]->getShortName(), (*molecules)[derivs->target(it)]->getShortName(), (*interactions)[it]->getName());
		fflush(dot);
}	
	fprintf(dot,"}\n");
	fflush(dot);

	pclose(dot);




}
