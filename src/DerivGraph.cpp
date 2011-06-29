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
	max_rate = 1.0;
	min_rate = .05;
//	maxBasic = 5;
//	maxComp = 2;
//	maxProm = 2;

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
  
  /*
   * Set up the list vectors.
   *
   * These vectors hold references to molecules and of the same type, which are added manually in the DerivGraph mutation
   * methods.
   *
   * The vector lists are convenient when trying to iterate a selected type or get a random member, and are more
   * efficient than iterating the entire graph.
   *
   * These vectors must be deleted in the destructor. The member objects are all pointers to objects deleted elsewhere.
   *
   * NOTE: If a molecule is deleted from the derivgraph, references to that molecule should be removed from the vectors
   */

    //all molecules are added to this list
    MoleculeList = new vector<Molecule*>();
    t.trace("mloc","DerivGraph %u MoleculeList vector at %u\n", (unsigned int) this, (unsigned int) MoleculeList);
    
    //basic proteins are added to this list
    ProteinList = new vector<Protein*>();
    t.trace("mloc","DerivGraph %u ProteinList vector at %u\n", (unsigned int) this, (unsigned int) ProteinList);
    
    //mRNAs are added to this list
    mRNAList = new vector<mRNA*>();
    t.trace("mloc","DerivGraph %u mRNAList vector at %u\n", (unsigned int) this, (unsigned int) mRNAList);
    
    DNAList = new vector<DNA*>();
    t.trace("mloc","DerivGraph %u DNAList vector at %u\n", (unsigned int) this, (unsigned int) DNAList);
    
    ComplexList = new vector<Complex*>();
    t.trace("mloc","DerivGraph %u ComplexList vector at %u\n", (unsigned int) this, (unsigned int) ComplexList);

    PTMList = new vector<PTMProtein*>();
    t.trace("mloc","DerivGraph %u PTMList vector at %u\n", (unsigned int) this, (unsigned int) PTMList);

    InteractionList = new vector<Interaction*>();
    t.trace("mloc","DerivGraph %u InteractionList vector at %u\n", (unsigned int) this, (unsigned int) InteractionList);
    
    TranscriptionList = new vector<Transcription*>();
    t.trace("mloc","DerivGraph %u TranscriptionList vector at %u\n", (unsigned int) this, (unsigned int) TranscriptionList);
    
    TranslationList = new vector<Translation*>();
    t.trace("mloc","DerivGraph %u TranslationList vector at %u\n", (unsigned int) this, (unsigned int) TranslationList);
    
    DegradationList = new vector<Degradation*>();
    t.trace("mloc","DerivGraph %u DegradationList vector at %u\n", (unsigned int) this, (unsigned int) DegradationList);
    
    ForwardComplexationList = new vector<ForwardComplexation*>();
    t.trace("mloc","DerivGraph %u ForwardComplexationList vector at %u\n", (unsigned int) this, (unsigned int) ForwardComplexationList);
    
    ReverseComplexationList = new vector<ReverseComplexation*>();
    t.trace("mloc","DerivGraph %u ReverseComplexationList vector at %u\n", (unsigned int) this, (unsigned int) ReverseComplexationList);

    ForwardPTMList = new vector<ForwardPTM*>();
    t.trace("mloc","DerivGraph %u ForwardPTMList vector at %u\n", (unsigned int) this, (unsigned int) ForwardPTMList);

    ReversePTMList = new vector<ReversePTM*>();
    t.trace("mloc","DerivGraph %u ReversePTMList vector at %u\n", (unsigned int) this, (unsigned int) ReversePTMList);

    PromoterBindList = new vector<PromoterBind*>();
    t.trace("mloc","DerivGraph %u PromoterBindList vector at %u\n", (unsigned int) this, (unsigned int) PromoterBindList);


    t.trace("init","New DerivGraph created\n");

    count=0;

    nullnode = add(new NullNode());
    (*molecules)[nullnode]->setID(count++);




/** Test code /
    newBasic();
    newBasic();
    newBasic();
    newPTM();
    newPTM();
    newPTM();
    newPTM();
    newPTM();
    newPTM();
    newPTM();
    newPTM();
    newPTM();
    newPTM();
    newPTM();
    newPTM();
    newPTM();
    newPTM();
    newPTM();
//test();


	rungeKuttaEvaluate(.5);	
*/
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

   //delete the various molecule lists
   delete ProteinList;
   delete mRNAList;
   delete DNAList;
   delete ComplexList;
   delete MoleculeList;
   delete PTMList;

   //delete the various interaction lists
   delete InteractionList;
   delete TranscriptionList;
   delete TranslationList;
   delete DegradationList;
   delete ForwardComplexationList;
   delete ReverseComplexationList;
   delete ForwardPTMList;
   delete ReversePTMList;
   delete PromoterBindList;


   //delete all Molecule objects mapped by Nodes
   t.trace("free","Deleting members of NodeMap at location %u\n",(unsigned int) molecules);
   for(ListDigraph::NodeIt it(*derivs); it !=INVALID; ++it){
	
	t.trace("free","Deleting NodeMap member at location %u\n",(unsigned int) (*molecules)[it]);
	
	//output some information about the molecule being deleted
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
	
	
	ListDigraph::Node A = add(new Molecule());
	ListDigraph::Node B = add(new Molecule());
	ListDigraph::Node C = add(new Molecule());
	(*molecules)[A]->setID(count++);
	(*molecules)[B]->setID(count++);
	(*molecules)[C]->setID(count++);

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

	rungeKuttaEvaluate(1.0);

return;
/*
	count = 1;
	ListDigraph::Node A = add(new DNA());
	(*molecules)[A]->setID(count++);
	ListDigraph::Node B = add(new mRNA());
	(*molecules)[B]->setID((*molecules)[A]->getID());
	ListDigraph::Node C = add(new Protein());
	(*molecules)[C]->setID((*molecules)[B]->getID());
	ListDigraph::Node D = add(new NullNode());
	(*molecules)[D]->setID(0);

	ListDigraph::Node E = add(new DNA());
	(*molecules)[E]->setID(count++);
	ListDigraph::Node F = add(new mRNA());
	(*molecules)[F]->setID((*molecules)[E]->getID());
	ListDigraph::Node G = add(new Protein());
	(*molecules)[G]->setID((*molecules)[F]->getID());

	(*molecules)[A]->setValue(1);
	(*molecules)[B]->setValue(10);
	(*molecules)[C]->setValue(5);

	(*molecules)[E]->setValue(1);
	(*molecules)[F]->setValue(10);
	(*molecules)[G]->setValue(5);



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
	
	(*molecules)[CG]->setID(count++);

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

//	float rkStep = 1.0;
  
//	rungeKuttaEvaluate(rkStep); 
*/
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

	//reset the runge-kutta internal variables for all molecules
	for(ListDigraph::NodeIt it(*derivs); it != INVALID; ++it)
		(*molecules)[it]->reset();

	

	//time loop
	for(float i = 0; i< 10; i+=rkStep){
			

		//each iteration of this loop refines the approximation based on the previous calculations	
		for(int k = 0; k<4; k++){
			
			//for every interaction in the graph
			for(ListDigraph::ArcIt it(*derivs); it != INVALID; ++it){
				
				//calculate the effect this interaction has on the source molecule
				(*molecules)[derivs->source(it)]->updateRkVal(k, getEffect(derivs->source(it), it, k, rkStep));
				//calculate the effect this interaction has on the target molecule
				(*molecules)[derivs->target(it)]->updateRkVal(k, getEffect(derivs->target(it), it, k, rkStep));

			}
		}

		//after the four rkVals are calcualted for all molecules, the next point can be computed
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
 *
 * The effect of an Interaction is defined by the Interaction::getEffect method.
 * The result of this method will be a positive or negative change in concentration.
 *
 * @param m The Node containing the Molecule being affected
 * @param i The Arc containing the Interaction taking place
 * @param rkIteration The current iteration of the Runge-Kutta algorithm
 * @param rkStep The timestep advanced by Runge-Kutta
 *
 * @return  the effect a particular Arc (Interaction) will have on a Node (Molecule)
 *
 */
float DerivGraph::getEffect(ListDigraph::Node m, ListDigraph::Arc i, int rkIteration, float rkStep){

	//get the effect the chosen interaction will have on the chosen molecule
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
	
	//store the nodeid in the molecule
	//this allow sthe molecule to find its position in the graph structure
	(*molecules)[newNode]->nodeID = derivs->id(newNode);
	
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
	
	//store the arcid in the interaction
	//this allows the interaction to find its position in the graph structure
	(*interactions)[newArc]->arcID = derivs->id(newArc);

	//return the newly created Arc
	return newArc;
}

/**
 * DerivGraph::newBasic()
 *
 * Create a new DNA, mRNA, and protein in the cell.
 *
 * DNA ---> mRNA ----> Protein
 *            |           |
 *            v           v
 *           Deg         Deg
 *
 *
 *
 */
void DerivGraph::newBasic(){

	t.trace("mutate","DerivGraph %u, new Basic Protein\n",(unsigned int)this);
/*
	if (DNAList->size() >= maxBasic)
	{
		t.trace("mutate","Basic Protein count is at limit\n");
		return;
	}
*/
	//create a new DNA, MRNA, and Protein
	ListDigraph::Node d = add(new DNA());
	ListDigraph::Node m = add(new mRNA());
	ListDigraph::Node p = add(new Protein());

	//create the interactions between the newly created basic system
	ListDigraph::Arc txn = add(new Transcription(), d, m);
	ListDigraph::Arc tsln = add(new Translation(), m, p);
	ListDigraph::Arc mdeg = add(new Degradation(), m, nullnode);
	ListDigraph::Arc pdeg = add(new Degradation(), p, nullnode);
	
	
	//give the DNA the next unique ID, and assign the same ID to the rest of the system
	// ex. d5 -> m5 -> p5   instead of d5 -> m6 -> p7
	(*molecules)[d]->setID(count++);
	(*molecules)[m]->setID((*molecules)[d]->getID());
	(*molecules)[p]->setID((*molecules)[d]->getID());


	//add molecule references to the appropriate lists
	DNAList->push_back( (DNA*) (*molecules)[d]);
	mRNAList->push_back( (mRNA*) (*molecules)[m]);
	ProteinList->push_back( (Protein*) (*molecules)[p]);
	MoleculeList->push_back( (*molecules)[d]);
	MoleculeList->push_back( (*molecules)[m]);
	MoleculeList->push_back( (*molecules)[p]);


	//add interaction references to the appropriate lists
	TranscriptionList->push_back( (Transcription*) (*interactions)[txn]);
	TranslationList->push_back( (Translation*) (*interactions)[tsln]);
	DegradationList->push_back( (Degradation*) (*interactions)[mdeg]);
	DegradationList->push_back( (Degradation*) (*interactions)[pdeg]);


	//the listsizes help to easily verify objects were created and added
	t.trace("mutate","DNAList.size() = %d\n", DNAList->size());
	t.trace("mutate","mRNAList.size() = %d\n", mRNAList->size());
	t.trace("mutate","ProteinList.size() = %d\n", ProteinList->size());

	t.trace("mutate","TranscriptionList.size() = %d\n", TranscriptionList->size());
	t.trace("mutate","TranslationList.size() = %d\n", TranslationList->size());
	t.trace("mutate","DegradationList.size() = %d\n", DegradationList->size());

}

/**
 * Randomly select a forward interaction and modify its rate.
 *
 *
 */
void DerivGraph::forwardRateChange(){

	//get the total number of forward interactions
	int totalSize = 0;
	totalSize += TranslationList->size();
	totalSize += ForwardComplexationList->size();
	totalSize += ForwardPTMList->size();
	
	t.trace("mutate","size = %d (%d + %d)\n", totalSize-1, TranslationList->size(), ForwardComplexationList->size());

	Interaction* selectedInteraction;
	
	//select a random integer between 0 and the total number of forward interactions
	unsigned int selectedIndex = r.randInt(totalSize - 1);
	t.trace("mutate","selectedIndex = %d\n", selectedIndex);

	//
	if(selectedIndex < TranslationList->size())
	{
		t.trace("mutate","TranslationList[%d]\n", selectedIndex);
		selectedInteraction = (*TranslationList)[selectedIndex];
	}
	else if(selectedIndex >= TranslationList->size())
	{
		selectedIndex -= TranslationList->size();
		t.trace("mutate","ForwardComplexation[%d]\n",selectedIndex);
		selectedInteraction = (*ForwardComplexationList)[selectedIndex];
	}

	ListDigraph::Arc selectedArc= derivs->arcFromId(selectedInteraction->arcID);

	Molecule* source = (*molecules)[derivs->source(selectedArc)];
	Molecule* target = (*molecules)[derivs->target(selectedArc)];

float newRate = min_rate + r.rand(max_rate - min_rate);
	
	t.trace("mutate","%s -> %s new rate: %f (old rate: %f)\n",source->getShortName(), target->getShortName(), newRate, selectedInteraction->getRate());

	selectedInteraction->setRate(newRate);

}

void DerivGraph::reverseRateChange(){

	int totalSize = 0;
	totalSize += ReverseComplexationList->size();
	totalSize += ReversePTMList->size();
	
	if(totalSize < 1)
	{
		t.trace("mutate","Reverse rate change failure: no reverse rates\n");
		return;
	}

	t.trace("mutate","size = %d (%d + %d)\n", totalSize-1, ReverseComplexationList->size(), ReversePTMList->size());

	Interaction* selectedInteraction;

	unsigned int selectedIndex = r.randInt(totalSize - 1);
	t.trace("mutate","selectedIndex = %d\n", selectedIndex);
	
	if(selectedIndex < ReverseComplexationList->size())
	{
		t.trace("mutate","ReverseComplexationList[%d]\n", selectedIndex);
		selectedInteraction = (*ReverseComplexationList)[selectedIndex];
	}
	else if(selectedIndex >= ReverseComplexationList->size())
	{
		selectedIndex -= ReverseComplexationList->size();
		t.trace("mutate","ReversePTM[%d]\n",selectedIndex);
		selectedInteraction = (*ReversePTMList)[selectedIndex];
	}

	ListDigraph::Arc selectedArc= derivs->arcFromId(selectedInteraction->arcID);

	Molecule* source = (*molecules)[derivs->source(selectedArc)];
	Molecule* target = (*molecules)[derivs->target(selectedArc)];
	float newRate = min_rate + r.rand(max_rate - min_rate);
	t.trace("mutate","%s -> %s new rate: %f (old rate: %f)\n",source->getShortName(), target->getShortName(), newRate, selectedInteraction->getRate());

	selectedInteraction->setRate(newRate);
}

void DerivGraph::degradationRateChange(){
	
	int selectedIndex = r.randInt(DegradationList->size()-1);
	float newRate = min_rate + r.rand(max_rate - min_rate);

	Interaction* selectedInteraction = (*DegradationList)[selectedIndex];

	ListDigraph::Arc selectedArc= derivs->arcFromId(selectedInteraction->arcID);

	Molecule* source = (*molecules)[derivs->source(selectedArc)];
	Molecule* target = (*molecules)[derivs->target(selectedArc)];


	t.trace("mutate","%s -> %s new rate: %f (old rate: %f)\n",source->getShortName(), target->getShortName(), newRate, selectedInteraction->getRate());
	selectedInteraction->setRate(newRate);

}

DNA* DerivGraph::histoneMod(){

	//choose a random index from the DNAList vector	
	int selectedIndex = r.randInt(DNAList->size()-1);
	//choose a new value for the histone mod [0,2]
	float newHistoneModValue = r.rand(2.0);

	//set the selected DNA's histone mod value to the new number
	(*DNAList)[selectedIndex]->setHistoneModValue(newHistoneModValue);
	t.trace("mutate","Histone Mod: DNAList[%d] -> %s. New Value = %f\n",selectedIndex, (*DNAList)[selectedIndex]->getShortName(), newHistoneModValue);
	return (*DNAList)[selectedIndex];
}

void DerivGraph::newPTM(){


	int PTMSelected = -1;

	int totalSize = 0;
	totalSize += ProteinList->size();
	totalSize += PTMList->size();
	
	t.trace("mutate","size = %d (%d + %d)\n", totalSize-1, ProteinList->size(), PTMList->size());

	Molecule* selectedMolecule;

	unsigned int selectedIndex = r.randInt(totalSize - 1);
	t.trace("mutate","selectedIndex = %d\n", selectedIndex);
	
	if(selectedIndex < ProteinList->size())
	{
		t.trace("mutate","ProteinList[%d]\n", selectedIndex);
		selectedMolecule = (Molecule*) (*ProteinList)[selectedIndex];
		PTMSelected = 0;
	}
	else if(selectedIndex >= ProteinList->size())
	{
		selectedIndex -= ProteinList->size();
		t.trace("mutate","PTMList[%d]\n",selectedIndex);
		selectedMolecule = (Molecule*) (*PTMList)[selectedIndex];
		PTMSelected = 1;
	}


	ListDigraph::Node selectedNode = derivs->nodeFromId(selectedMolecule->nodeID);
	
	ListDigraph::Node newPTM;
	
	newPTM = add(new PTMProtein());
	
	//copy the ptm counts to the new PTMProtein
	for(int i = 0; i< 4; i++)
		((PTMProtein*)(*molecules)[newPTM])->setPTMCount(i, selectedMolecule->getPTMCount(i));

	//add the new PTM to the relevant lists
	PTMList->push_back( (PTMProtein*) (*molecules)[newPTM]);
	MoleculeList->push_back( (*molecules)[newPTM]);
	
	(*molecules)[newPTM]->setID(count++);
	
	((PTMProtein*)(*molecules)[newPTM])->addRandPTM(r.randInt(3));
	
	ListDigraph::Arc PTM_f = add(new ForwardPTM(), selectedNode, newPTM);
	ListDigraph::Arc PTM_r = add(new ReversePTM(), newPTM, selectedNode);
	ListDigraph::Arc PTM_d = add(new Degradation(), newPTM, nullnode);

	DegradationList->push_back( (Degradation*) (*interactions)[PTM_d]);
	ForwardPTMList->push_back( (ForwardPTM*) (*interactions)[PTM_f]);
	ReversePTMList->push_back( (ReversePTM*) (*interactions)[PTM_r]);

	t.trace("mutate","OldPTM: %s\n",(PTMProtein*) selectedMolecule->getLongName());
	t.trace("mutate","NewPTM: %s\n",(PTMProtein*) (*molecules)[newPTM]->getLongName());

}

void DerivGraph::newComplex(){
/*
	if(ComplexList->size() >= maxComp)
	{
		t.trace("mutate","Total Complex protein count is at limit\n");
		return;

	}
*/

	unsigned int i1 = -1;
	unsigned int i2 = -1;
	Protein* p1;
	Protein* p2;
	int totalSize = 0;
	totalSize += ProteinList->size();
	totalSize += ComplexList->size();

	if(totalSize < 2)
	{
		t.trace("mutate","New Complex failed. Not enough proteins\n");
		return;
	}
	while(i1 == i2)
	{
		i1 = r.randInt(totalSize-1);
		i2 = r.randInt(totalSize-1);

	}

	t.trace("mutate","Complex: %d - %d (%d + %d)\n",i1,i2,ProteinList->size(), ComplexList->size());


	if(i1 < ProteinList->size())
	{
		t.trace("mutate","ProteinList[%d]\n", i1);
		p1 = (*ProteinList)[i1];
	}
	else if(i1 >= ProteinList->size())
	{
		i1 -= ProteinList->size();
		t.trace("mutate","ComplexList[%d]\n",i1);
		p1 = (*ComplexList)[i1];

	}
	
	if(i2 < ProteinList->size())
	{
		t.trace("mutate","ProteinList[%d]\n", i2);
		p2 = (*ProteinList)[i2];
	}
	else if(i2 >= ProteinList->size())
	{
		i2 -= ProteinList->size();
		t.trace("mutate","ComplexList[%d]\n",i2);
		p2 = (*ComplexList)[i2];

	}
	
	int id1 = p1->nodeID;
	int id2 = p2->nodeID;
	ListDigraph::Node n1 = derivs->nodeFromId(id1);
	ListDigraph::Node n2 = derivs->nodeFromId(id2);


	int a = 0;
	int b = 0;
	for(unsigned int c = 0; c < ComplexList->size(); c++)
	{
		a = (*ComplexList)[c]->getComponentId(1);
		b = (*ComplexList)[c]->getComponentId(2);
		if( (a == id1 && b == id2) || (a == id2 && b == id1))
		{
			t.trace("mutate","New Complex failed. Already Exists\n");
			return;
		}
	}


	ListDigraph::Node comp = add(new Complex(id1, id2));
	(*molecules)[comp]->setID(count++);

	ListDigraph::Arc f1 = add(new ForwardComplexation(id1, id2), n1, comp); 
	ListDigraph::Arc f2 = add(new ForwardComplexation(id1, id2), n2, comp); 
	ListDigraph::Arc r1 = add(new ReverseComplexation(id1, id2), comp, n1); 
	ListDigraph::Arc r2 = add(new ReverseComplexation(id1, id2), comp, n2); 
	ListDigraph::Arc deg = add(new Degradation(), comp, nullnode);	

	ComplexList->push_back( (Complex*) (*molecules)[comp]);
	MoleculeList->push_back((*molecules)[comp]);

	ForwardComplexationList->push_back( (ForwardComplexation*) (*interactions)[f1]);
	ForwardComplexationList->push_back( (ForwardComplexation*) (*interactions)[f2]);
	ReverseComplexationList->push_back( (ReverseComplexation*) (*interactions)[r1]);
	ReverseComplexationList->push_back( (ReverseComplexation*) (*interactions)[r2]);
	DegradationList->push_back( (Degradation*) (*interactions)[f1]);

	t.trace("mutate","New complex created\n");


}

void DerivGraph::newPromoter(){
/*
	if(PromoterBindList->size() >= maxProm){
		t.trace("mutate","Promoter count is at limit\n");
		return;
	}
*/
	int selectionIndex = r.randInt(DNAList->size() -1);
	if( (*DNAList)[selectionIndex]->promoterId >= 0)
	{
		t.trace("mutate","New Promoter Failed: already taken\n");
		return;
	}
	DNA* d = (*DNAList)[selectionIndex];

	
	int selectionIndex2 = r.randInt(ProteinList->size() -1);
	Protein* p = (*ProteinList)[selectionIndex2];
	
	ListDigraph::Node nd = derivs->nodeFromId(d->nodeID);
	ListDigraph::Node np = derivs->nodeFromId(p->nodeID);
	
	float fwd = 0;
	float rev = 1;
	while(rev > fwd)
	{
		fwd = min_rate + r.rand(max_rate - min_rate);
		rev = min_rate + r.rand(max_rate - min_rate);
	}
	t.trace("mutate","gene: %s protein: %s kf: %f kr: %f\n",d->getShortName(),p->getShortName(), fwd, rev);

	ListDigraph::Arc a = add(new PromoterBind(fwd, rev), np, nd);
	d->promoterId = derivs->id(a);

	PromoterBindList->push_back( (PromoterBind*) (*interactions)[a]);
	

}

/**
 * void DerivGraph::outputDotImage(int, int)
 *
 * Output a png image of the current graph structure using GraphViz.
 *
 * A process running GraphViz is forked and a pipe opened to its standard in.
 * The general layout of the output file can be changed below.
 * The Node and Arc names are defined within the Molecule and Interaction classes.
 *
 * @param cellNum the cell number to put in the filename
 * @param gen the generation number to put in the filename
 *
 */
void DerivGraph::outputDotImage(int cellNum, int gen){
	
	char buf[200];
	sprintf(buf, "dot -Gsize=\"20,20\" -Tpng -o../output/Cell%dGen%d.png",cellNum, gen);

	//popen forks and execs and returns a pipe to the new process stdin
	FILE* dot = popen(buf,"w");
	
	fprintf(dot,"digraph mol_interactions {\n");
	fflush(dot);

	fprintf(dot,"size=\"8,5\"\n");
	fflush(dot);

	fprintf(dot,"node [shape = circle];\n");
	fflush(dot);

	fprintf(dot,"edge [len =2 ] ;\n");
	fflush(dot);

	//iterate all of the Arcs and add them to the visualization. Nodes are implicitly defined by the source and target of the interactions.
	for(ListDigraph::ArcIt it(*derivs); it != INVALID; ++it){
		fprintf(dot, "%s -> %s [ label = \"%s (%f)\"];\n",(*molecules)[derivs->source(it)]->getShortName(), (*molecules)[derivs->target(it)]->getShortName(), (*interactions)[it]->getName(), (*interactions)[it]->getRate());
		fflush(dot);
}	

	fprintf(dot, "overlap=scale\n");
	fflush(dot);

	fprintf(dot,"}\n");
	fflush(dot);

	//close cleanly
	pclose(dot);

}
void DerivGraph::outputDataPlot(int cellNum, int gen, float step){
	
	FILE* gnuplot = popen("gnuplot","w");

	//FILE* gnuplot = fopen("test.txt","w");
	for(unsigned int i =  0; i < MoleculeList->size(); i++){
	
		fprintf(gnuplot, "set term png\n");
		fflush(gnuplot);
	
		fprintf(gnuplot, "set xlabel \"time\"\n");
		fflush(gnuplot);

		fprintf(gnuplot, "set ylabel \"concentration\"\n");
		fflush(gnuplot);

		fprintf(gnuplot, "set format x \"%c03.2f\"\n",'%');
		fflush(gnuplot);

		fprintf(gnuplot, "set output \"../output/%sc%dg%d.plot.png\"\n",(*MoleculeList)[i]->getShortName(),cellNum, gen);
		fflush(gnuplot);
	
		fprintf(gnuplot, "plot \\\n");
		fflush(gnuplot);
	
		fprintf(gnuplot, "\"-\" using 2:($1==%d ? $3 : 1/0) t \"%s\" pt 1 \n",i,(*MoleculeList)[i]->getLongName());
		fflush(gnuplot);
	
		float t;
		t = 0;
		for(unsigned int j = 0; j < (*MoleculeList)[i]->getRungeKuttaSolution()->size(); j++)
		{	
			float k = (*MoleculeList)[i]->getRungeKuttaSolution()->at(j);
			fprintf(gnuplot, "%d %f %f\n",i, t, k);
			fflush(gnuplot);
			t+=step;
		}	
	
		fprintf(gnuplot, "exit\n");
		fflush(gnuplot);
	}
	pclose(gnuplot);
}
