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
    
    t.trace("mloc","DerivGraph location at %p\n",  this);

    //create the directed graph
    derivs = new ListDigraph();
    t.trace("mloc","DerivGraph %p ListDigraph location at %p\n",   this,   derivs);

    //map molecules onto the nodes
    molecules = new ListDigraph::NodeMap<Molecule*>(*derivs);
    t.trace("mloc","DerivGraph %p NodeMap location at %p\n",   this,   molecules);
    
    //map interactions onto the arcs
    interactions = new ListDigraph::ArcMap<Interaction*>(*derivs);
    t.trace("mloc","DerivGraph %p ArcMap location at %p\n",   this,   interactions);
  
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
    t.trace("mloc","DerivGraph %p MoleculeList vector at %p\n",   this,   MoleculeList);
    
    //basic proteins are added to this list
    ProteinList = new vector<Protein*>();
    t.trace("mloc","DerivGraph %p ProteinList vector at %p\n",   this,   ProteinList);
    
    //mRNAs are added to this list
    mRNAList = new vector<mRNA*>();
    t.trace("mloc","DerivGraph %p mRNAList vector at %p\n",   this,   mRNAList);
    
    DNAList = new vector<DNA*>();
    t.trace("mloc","DerivGraph %p DNAList vector at %p\n",   this,   DNAList);
    
    ComplexList = new vector<Complex*>();
    t.trace("mloc","DerivGraph %p ComplexList vector at %p\n",   this,   ComplexList);

    PTMList = new vector<PTMProtein*>();
    t.trace("mloc","DerivGraph %p PTMList vector at %p\n",   this,   PTMList);

    InteractionList = new vector<Interaction*>();
    t.trace("mloc","DerivGraph %p InteractionList vector at %p\n",   this,   InteractionList);
    
    TranscriptionList = new vector<Transcription*>();
    t.trace("mloc","DerivGraph %p TranscriptionList vector at %p\n",   this,   TranscriptionList);
    
    TranslationList = new vector<Translation*>();
    t.trace("mloc","DerivGraph %p TranslationList vector at %p\n",   this,   TranslationList);
    
    DegradationList = new vector<Degradation*>();
    t.trace("mloc","DerivGraph %p DegradationList vector at %p\n",   this,   DegradationList);
    
    ForwardComplexationList = new vector<ForwardComplexation*>();
    t.trace("mloc","DerivGraph %p ForwardComplexationList vector at %p\n",   this,   ForwardComplexationList);
    
    ReverseComplexationList = new vector<ReverseComplexation*>();
    t.trace("mloc","DerivGraph %p ReverseComplexationList vector at %p\n",   this,   ReverseComplexationList);

    ForwardPTMList = new vector<ForwardPTM*>();
    t.trace("mloc","DerivGraph %p ForwardPTMList vector at %p\n",   this,   ForwardPTMList);

    ReversePTMList = new vector<ReversePTM*>();
    t.trace("mloc","DerivGraph %p ReversePTMList vector at %p\n",   this,   ReversePTMList);

    PromoterBindList = new vector<PromoterBind*>();
    t.trace("mloc","DerivGraph %p PromoterBindList vector at %p\n",   this,   PromoterBindList);


    t.trace("init","New DerivGraph created\n");

    //the count is used to name new molecules
    count=0;

    //create the null node targeted by degradation interactions
    nullnode = add(new NullNode());
    (*molecules)[nullnode]->setID(count++);



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
   t.trace("free","Deleting members of NodeMap at location %p\n",  molecules);
   for(ListDigraph::NodeIt it(*derivs); it !=INVALID; ++it){
	
	t.trace("free","Deleting NodeMap member at location %p\n",  (*molecules)[it]);
	
	//output some information about the molecule being deleted
	t.trace("free","longnname: %s\n", (*molecules)[it]->getLongName());
	t.trace("free","shortname: %s\n", (*molecules)[it]->getShortName());
	
	delete (*molecules)[it];
   }
  
   //delete the Molecule NodeMap
   t.trace("free","Deleting NodeMap object at location %p\n",molecules);
   delete molecules;

   //delete all Interaction objects mapped by Arcs
   t.trace("free","Deleting members of ArcMap at location %p\n", interactions);
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
	
	
	ListDigraph::Node A = add(new Protein());
	ListDigraph::Node B = add(new Protein());
	ListDigraph::Node C = add(new Protein());
	(*molecules)[A]->setID(count++);
	(*molecules)[B]->setID(count++);
	(*molecules)[C]->setID(count++);

	(*molecules)[A]->setValue(2);
	(*molecules)[B]->setValue(4);
	(*molecules)[C]->setValue(1);

	ListDigraph::Arc AB = add(new ForwardPTM(), A, B);
	ListDigraph::Arc AC = add(new ForwardPTM(), A, C);
	ListDigraph::Arc BC = add(new ForwardPTM(), B, C);
	ListDigraph::Arc CB = add(new ForwardPTM(), C, B);

	(*interactions)[AB]->setRate(.01);
	(*interactions)[AC]->setRate(.03);
	(*interactions)[BC]->setRate(.07);
	(*interactions)[CB]->setRate(.09);

	rungeKuttaEvaluate(rkTimeStep, rkTimeLimit);

return;

}

/**
 * Uses the Runge-Kutta fourth order method to approximate the solutions to the system of differential equations
 *
 * The result of this algorithm is the vector rungeKuttaSolution within each Molecule object containing the approximation of
 * the concentration at each timestep.
 *
 * @param rkStep the timestep (precision) between calculated points
 */
void DerivGraph::rungeKuttaEvaluate(float rkStep, float rkLimit){

	//reset the runge-kutta internal variables for all molecules
	for(ListDigraph::NodeIt it(*derivs); it != INVALID; ++it)
		(*molecules)[it]->reset();

	

	//time loop
	for(float i = 0; i< rkLimit; i+=rkStep){

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

	//test output, display the values calculated by runge kutta for each molecule to stdout
	//for(ListDigraph::NodeIt it(*derivs); it != INVALID; ++it){
	//	(*molecules)[it]->outputRK();	
	//}
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


	(*molecules)[newNode]->setValue(defaultInitialConcentration);

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

	(*interactions)[newArc]->setRate(minKineticRate + r.rand(maxKineticRate-minKineticRate));
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

	t.trace("mutate","DerivGraph %p, new Basic Protein\n", this);
	if ((int)DNAList->size() >= maxBasic)
	{
		t.trace("mutate","Basic Protein count is at limit\n");
		return;
	}
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
	MoleculeList->push_back( (*molecules)[d]);

	mRNAList->push_back( (mRNA*) (*molecules)[m]);
	MoleculeList->push_back( (*molecules)[m]);

	ProteinList->push_back( (Protein*) (*molecules)[p]);
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
 * void DerivGraph::forwardRateChange()
 *
 * Randomly select a forward interaction and modify its rate.
 * 
 * Forward interactions are of type Translation, ForwardComplex
 * 
 *
 * Random Selection across multiple arrays logic:
 * Arr1: [0 1 2 3] (size = 4)
 * Arr2: [0 1 2]   (size = 3)
 * Arr3: [0 1]     (size = 2)
 *
 * Implicitly: [0 1 2 3|0 1 2|0 1] (totalSize = 9)
 *(subscripts)  0 1 2 3 4 5 6 7 8
 *
 * 1. Select random number between 0 and totalSize-1
 *
 * 2. Is number >= 0 and < Arr1.size? 
 *     (yes) Select Arr1[number]
 *     (no) continue
 *
 * 3. Is number >= Arr1.size and < Arr1.size + Arr2.size?
 *     (yes) Select Arr2[number - Arr1.size]
 *     (no) continue
 *
 * 4. Is number >= Arr1.size + Arr2.size and < Arr1.size + Arr2.size + Arr3.size?
 *     (yes) Select Arr3[number - Arr1.size - Arr2.size]
 *
 * Can extend to any number of arrays 
 */
void DerivGraph::forwardRateChange(){

	//get the total number of forward interactions
	int totalSize = 0;
	totalSize += TranslationList->size();
	totalSize += ForwardComplexationList->size();
	totalSize += ForwardPTMList->size();
	
	t.trace("mutate","size = %d (%d + %d + %d)\n", totalSize-1, TranslationList->size(), ForwardComplexationList->size(), ForwardPTMList->size());

	Interaction* selectedInteraction;
	
	//if a forwardComplexationInteraction is chosen, its pair reaction must be changed too
	int complexInteractionPairID = 0;

	
	//select a random integer between 0 and the total number of forward interactions
	unsigned int randIndex = r.randInt(totalSize - 1);
	t.trace("mutate","randIndex = %d\n", randIndex);

	//index falls within the TranslationList
	if(randIndex >= 0 && randIndex < TranslationList->size())
	{
		t.trace("mutate","TranslationList[%d]\n", randIndex);
		selectedInteraction = (*TranslationList)[randIndex];
	}
	//index falls within the ForwardComplexationList
	else if(randIndex >= TranslationList->size() && randIndex < TranslationList->size() + ForwardComplexationList->size())
	{
		t.trace("mutate","ForwardComplexation[%d]\n",randIndex - TranslationList->size());
		selectedInteraction = (*ForwardComplexationList)[randIndex - TranslationList->size()];
		complexInteractionPairID = ((ForwardComplexation*) selectedInteraction)->pairArcID;
	}
	//index falls within the ForwardPTMList
	else if(randIndex >= TranslationList->size() + ForwardComplexationList->size() && randIndex < TranslationList->size() + ForwardComplexationList->size() + ForwardPTMList->size())
	{
		t.trace("mutate","ForwardPTM[%d]\n",randIndex - TranslationList->size() - ForwardComplexationList->size());
		selectedInteraction = (*ForwardPTMList)[randIndex - TranslationList->size() - ForwardComplexationList->size()];
	}
		
	//get the arc which holds the interaction
	ListDigraph::Arc selectedArc= derivs->arcFromId(selectedInteraction->arcID);

	Molecule* source = (*molecules)[derivs->source(selectedArc)];
	Molecule* target = (*molecules)[derivs->target(selectedArc)];

	//get a random kinetic rate
	float newRate = minKineticRate + r.rand(maxKineticRate - minKineticRate);
	
	t.trace("mutate","%s -> %s new rate: %f (old rate: %f)\n",source->getShortName(), target->getShortName(), newRate, selectedInteraction->getRate());

	//set the chosen interaction rate to the newly generated rate
	selectedInteraction->setRate(newRate);

	//if the selectedInteraction was a forward Complex, change the pair interaction so the rates remain the same	
	if(complexInteractionPairID){
		(*interactions)[derivs->arcFromId(complexInteractionPairID)]->setRate(newRate);

		t.trace("mutate","%s -> %s pair interaction also changed\n", (*molecules)[derivs->source(derivs->arcFromId(complexInteractionPairID))]->getShortName(), (*molecules)[derivs->target(derivs->arcFromId(complexInteractionPairID))]->getShortName());
	}

}
/**
 * void DerivGraph::reverseRateChange()
 *
 * Randomly select a reverse interaction and modify its rate.
 *
 * Reverse interactions are of type ReverseComplexation, ReversePTM
 * 
 * See DerivGraph::forwardRateChange() for detailed explaination of how the rate is randomly selected 
 */
void DerivGraph::reverseRateChange(){

	int totalSize = 0;
	totalSize += ReverseComplexationList->size();
	totalSize += ReversePTMList->size();
	
	//it is possible that no rates exist at this point
	if(totalSize < 1)
	{
		t.trace("mutate","Reverse rate change failure: no reverse rates\n");
		return;
	}

	t.trace("mutate","size = %d (%d + %d)\n", totalSize-1, ReverseComplexationList->size(), ReversePTMList->size());

	Interaction* selectedInteraction;

	int complexInteractionPairID = 0;

	//get a random number between 1 and the total number of reverse reactions
	unsigned int randIndex = r.randInt(totalSize - 1);
	t.trace("mutate","randIndex = %d\n", randIndex);

	//if the index falls within the ReverseComplexationList	
	if(randIndex >= 0 && randIndex < ReverseComplexationList->size())
	{
		t.trace("mutate","ReverseComplexationList[%d]\n", randIndex);
		selectedInteraction = (*ReverseComplexationList)[randIndex];
		complexInteractionPairID = ((ForwardComplexation*) selectedInteraction)->pairArcID;
	}
	//if the index falls within the ReversePTMList
	else if(randIndex >= ReverseComplexationList->size() && randIndex < ReverseComplexationList->size() + ReversePTMList->size())
	{
		t.trace("mutate","ReversePTM[%d]\n",randIndex - ReverseComplexationList->size());
		selectedInteraction = (*ReversePTMList)[randIndex - ReverseComplexationList->size()];
	}
	
	//get the Arc holding the interaction
	ListDigraph::Arc selectedArc= derivs->arcFromId(selectedInteraction->arcID);
	//get the source and target molecules
	Molecule* source = (*molecules)[derivs->source(selectedArc)];
	Molecule* target = (*molecules)[derivs->target(selectedArc)];
	
	//select a random rate between the minimum and maxium values
	float newRate = minKineticRate + r.rand(maxKineticRate - minKineticRate);
	t.trace("mutate","%s -> %s new rate: %f (old rate: %f)\n",source->getShortName(), target->getShortName(), newRate, selectedInteraction->getRate());
	
	//set the chosen interaction to the new rate
	selectedInteraction->setRate(newRate);
	
	//if the selectedInteraction was a reverse Complex, change the pair interaction so the rates remain the same	
	if(complexInteractionPairID){
		(*interactions)[derivs->arcFromId(complexInteractionPairID)]->setRate(newRate);

		t.trace("mutate","%s -> %s pair interaction also changed\n", (*molecules)[derivs->source(derivs->arcFromId(complexInteractionPairID))]->getShortName(), (*molecules)[derivs->target(derivs->arcFromId(complexInteractionPairID))]->getShortName());
	}

}
/**
 * void DerivGraph::degradationRateChange()
 *
 * Randomly select a degradation interaction and modify its rate.
 *
 * Degradation interactions are of type Degradation
 */
void DerivGraph::degradationRateChange(){
	
	int selectedIndex = r.randInt(DegradationList->size()-1);
	float newRate = minKineticRate + r.rand(maxKineticRate - minKineticRate);

	Interaction* selectedInteraction = (*DegradationList)[selectedIndex];

	ListDigraph::Arc selectedArc= derivs->arcFromId(selectedInteraction->arcID);

	Molecule* source = (*molecules)[derivs->source(selectedArc)];
	Molecule* target = (*molecules)[derivs->target(selectedArc)];


	t.trace("mutate","%s -> %s new rate: %f (old rate: %f)\n",source->getShortName(), target->getShortName(), newRate, selectedInteraction->getRate());
	selectedInteraction->setRate(newRate);

}
/**
 * void DerivGraph::histoneMod()
 *
 * Randomly select a DNA molecule, and set the Histone factor to a random value [0,2].
 *
 * The histone value is a constant multiplied factor applied to the rate of mRNA production
 * by a DNA molecule. It is initialized at 1.0, and is randomly assigned a value between [0,2].
 *
 * A value [0,1) results in repression of mRNA production.
 * A value (1,2] results in activation of mRNA production.
 */
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
/**
 * void DerivGraph::newPTM()
 *
 * Randomly select a Protein molecule to which a PTM should be applied. The new PTM Protein
 * has the same counts of modifications, with a random index incremented by one to reflect
 * the new value. A PTM can be applied to a basic protein or a previously existing PTM.
 *
 * TODO: check if a PTM already exists on a protein and prevent duplicate PTM's
 */
void DerivGraph::newPTM(){

	if((int)PTMList->size() >= maxPTM)
	{
		t.trace("mutate","PTM count is at limit\n");
		return;
	}

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
/**
 * void DerivGraph::newComplex()
 *
 * Randomly select two Protein molecules to be complexed together. 
 *
 * If the two selected proteins already exist in a complex reaction together, the mutation
 * will fail.
 *
 * Molecules which can complex together are Protein, and ComplexProteins
 *
 * TODO: add PTM's to the possible complex
 */
void DerivGraph::newComplex(){
	if((int)ComplexList->size() >= maxComp)
	{
		t.trace("mutate","Total Complex protein count is at limit\n");
		return;

	}

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

	//create the pair forward complexation interactions
	ListDigraph::Arc f1 = add(new ForwardComplexation(id1, id2), n1, comp); 
	ListDigraph::Arc f2 = add(new ForwardComplexation(id1, id2), n2, comp); 
	
	//set each interactions pairArcID, so changes to one can easily be made to the other
	((ForwardComplexation*)(*interactions)[f1])->setPairArcID(derivs->id(f2));
	((ForwardComplexation*)(*interactions)[f2])->setPairArcID(derivs->id(f1));
	
	t.trace("mutate","f1 arc id: %d, f2 arc id: %d\n", derivs->id(f1), derivs->id(f2));
	t.trace("mutate","f1 pair arc: %d, f2 pair arc: %d\n", ((ForwardComplexation*)(*interactions)[f1])->pairArcID, ((ForwardComplexation*)(*interactions)[f2])->pairArcID);

	//create the pair of reverse complexation interactions
	ListDigraph::Arc r1 = add(new ReverseComplexation(id1, id2), comp, n1); 
	ListDigraph::Arc r2 = add(new ReverseComplexation(id1, id2), comp, n2); 
	
	//set each interaction's pairArcID, so changes to one c an easily be made to the other
	((ReverseComplexation*)(*interactions)[r1])->setPairArcID(derivs->id(r2));
	((ReverseComplexation*)(*interactions)[r2])->setPairArcID(derivs->id(r1));
	
	t.trace("mutate","r1 arc id: %d, r2 arc id: %d\n", derivs->id(r1), derivs->id(r2));
	t.trace("mutate","r1 pair arc: %d, r2 pair arc: %d\n", ((ReverseComplexation*)(*interactions)[r1])->pairArcID, ((ReverseComplexation*)(*interactions)[r2])->pairArcID);

	
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

/**
 * void DerivGraph::newPromoter()
 *
 * Select a random protein and DNA to add a Protein-Promoter interaction to.
 *
 * The Protein-Promoter interaction is used in conjunction with the Hill model
 * of cooperativity, and affects the Goodwin term used by DNA to calculate 
 * the rate of production.
 *
 */
void DerivGraph::newPromoter(){

	if((int)PromoterBindList->size() >= maxProm){
		t.trace("mutate","Promoter count is at limit\n");
		return;
	}
	
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
		fwd = minKineticRate + r.rand(maxKineticRate - minKineticRate);
		rev = minKineticRate + r.rand(maxKineticRate - minKineticRate);
	}
	t.trace("mutate","gene: %s protein: %s kf: %f kr: %f\n",d->getShortName(),p->getShortName(), fwd, rev);

	ListDigraph::Arc a = add(new PromoterBind(fwd, rev), np, nd);
	d->promoterId = derivs->id(a);

	PromoterBindList->push_back( (PromoterBind*) (*interactions)[a]);
	

}

Molecule* DerivGraph::getBestMolecule(int CellID){


	//rungeKuttaEvaluate(rkTimeStep, rkTimeLimit);

	Molecule* bestMolecule = 0;
	int maxScore = -1;
	int s = -1;


	for(unsigned int i =  0; i < MoleculeList->size(); i++){
		s = (*MoleculeList)[i]->getScore();
		if(s > maxScore){
			maxScore = s;
			bestMolecule = (*MoleculeList)[i];
		}

	}
	t.trace("score","Cell %d best molecule is %s (%d)\n",CellID, bestMolecule->getShortName(), maxScore);
	return bestMolecule;
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
void DerivGraph::outputDotImage(const char* prefix, int pid, int cellNum, int gen){
	
	char buf[200];
	sprintf(buf, "dot -Gsize=\"20,20\" -Tpng -o%s/%d/cell%d/Cell%dGen%d.png",prefix, pid, cellNum, cellNum, gen);

	//popen forks and execs and returns a pipe to the new process stdin
	FILE* dot = popen(buf,"w");
	
	fprintf(dot,"digraph mol_interactions {\n");
	fflush(dot);

	//fprintf(dot,"size=\"8,5\"\n");
	//fflush(dot);

	fprintf(dot,"node [shape = ellipse];\n");
	fflush(dot);

	fprintf(dot,"edge [len =2 ] ;\n");
	fflush(dot);

	//iterate all of the Arcs and add them to the visualization. Nodes are implicitly defined by the source and target of the interactions.
	for(ListDigraph::ArcIt it(*derivs); it != INVALID; ++it){
		fprintf(dot, "\"%s (%d)\" -> \"%s (%d)\" [ label = \"%s (%f)\",  penwidth= %f];\n",(*molecules)[derivs->source(it)]->getShortName(),(*molecules)[derivs->source(it)]->getScore(), (*molecules)[derivs->target(it)]->getShortName(),(*molecules)[derivs->target(it)]->getScore(), (*interactions)[it]->getName(), (*interactions)[it]->getRate(), 0.5 + 2*(*interactions)[it]->getRate()/maxKineticRate);
		fflush(dot);
}	

	fprintf(dot, "overlap=scale\n");
	fflush(dot);

	fprintf(dot,"}\n");
	fflush(dot);

	//close cleanly
	pclose(dot);

}
/**
 * void DerivGraph::outputDataPlot(int, int, float)
 * 
 * Output a png image of the concentration data of molecules plotted by Gnuplot.
 *
 * For each molecule in the MoleculeList, a process running gnuplot is forked to which
 * data from Runge-Kutta is fed to produce a plot.
 * 
 * @param cellNum the cell number to put in the filename
 * @param gen the generation number to put in the filename
 * @param step the stepSize used between the rungeKuttaSolution data points
 */
void DerivGraph::outputDataPlot(const char* prefix, int pid, int cellNum, int gen, float step){
	
	FILE* gnuplot = popen("gnuplot > /dev/null 2>&1","w");

//	FILE* test = fopen("test.txt","w");
	for(unsigned int i =  0; i < MoleculeList->size(); i++){
		if((*MoleculeList)[i]->getScore() < 3)
			continue;
		fprintf(gnuplot, "set term png size 2048,1536\n");
		fflush(gnuplot);
	
		fprintf(gnuplot, "set xlabel \"time\"\n");
		fflush(gnuplot);

		fprintf(gnuplot, "set ylabel \"concentration\"\n");
		fflush(gnuplot);

		fprintf(gnuplot, "set format x \"%c03.2f\"\n",'%');
		fflush(gnuplot);

		fprintf(gnuplot, "set output \"%s/%d/cell%d/%sc%dg%d.plot.png\"\n",prefix, pid, cellNum,(*MoleculeList)[i]->getShortName(),cellNum, gen);
		fflush(gnuplot);
	
		fprintf(gnuplot, "plot \\\n");
		fflush(gnuplot);
	
		fprintf(gnuplot, "\"-\" using 2:($1==%d ? $3 : 1/0) t \"%s\" pt 1 with linespoints\n",i,(*MoleculeList)[i]->getLongName());
		fflush(gnuplot);
	
		float t;
		t = 0;
		for(unsigned int j = 0; j < (*MoleculeList)[i]->getRungeKuttaSolution()->size(); j++)
		{	
		
		if(j % 10 != 0){
			t+=step;
			continue;
		}
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

void DerivGraph::outputDataCsv(const char* prefix, int pid, int CellID, int currentGen, float rkTimeStep){



}

/**
 * DerivGraph::setLimits(int, int, int, int)
 *
 * Set the occurrence limits for mutation types
 *
 * @param max_basic Maximum basic proteins allowed
 * @param max_ptm Maximum number of post translationally modified proteins allowed
 * @param max_comp Maximum number of complexed proteins allowed
 * @param max_promoter Maximum number of protein-promoter interactions allowed
 */
void DerivGraph::setLimits(int max_basic, int max_ptm, int max_comp, int max_promoter){
	maxBasic = max_basic;
	maxPTM = max_ptm;
	maxComp = max_comp;
	maxProm = max_promoter;
	t.trace("args","Max Basic: %d\n", maxBasic);
	t.trace("args","Max PTM: %d\n", maxPTM);
	t.trace("args","Max Complex: %d\n", maxComp);
	t.trace("args","Max Promoters: %d\n", maxProm);

}
/**
 * DerivGraph::setRungeKuttaEval(float, float)
 *
 * Set the parameters for runge-kutta evalutation
 *
 * @param rk_time_step The timestep between points (t, conc) calculated by Runge-Kutta
 * @param rk_time_limit The upper time limit for Runge-Kutta calculation (time = x axis)
 */
void DerivGraph::setRungeKuttaEval(float rk_time_step, float rk_time_limit){
	rkTimeStep = rk_time_step;
	rkTimeLimit = rk_time_limit;
}

/**
 * DerivGraph::setKineticRateLimits(float, float)
 *
 * Assign lower and upper bounds to the randomly generated kinetic rates
 *
 * @param min_kinetic_rate The lower bound on random kinetic rates
 * @param max_kinetic_rate The upper bound on random kinetic rates
 */
void DerivGraph::setKineticRateLimits(float min_kinetic_rate, float max_kinetic_rate){
	minKineticRate = min_kinetic_rate;
	maxKineticRate = max_kinetic_rate;

}

/**
 * DerivGraph::setDefaultInitialConcentration(float)
 *
 * Set the default initial concentration for molecules
 *
 * @param initial_conc The initial concentration for new molecules
 */
void DerivGraph::setDefaultInitialConc(float initial_conc){

	defaultInitialConcentration = initial_conc;

}
