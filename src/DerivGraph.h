#ifndef DERIVGRAPH_H_
#define DERIVGRAPH_H_

#include "lemon/list_graph.h"
#include "lemon/concepts/maps.h"
#include <cstdio>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <typeinfo>
#include "MersenneTwister.h"

#include "Molecule.h"
#include "Interaction.h"

#include "CustomInteractions.h"
#include "CustomMolecules.h"

using namespace std;
using namespace lemon;


class DerivGraph{

public:
	DerivGraph();
	~DerivGraph();
	
	MTRand r;
	
	void test();
	void rungeKuttaEvaluate(float);
	void outputDotImage(int, int);
	void outputDataPlot(int, int, float);
	Molecule* getBestMolecule(int);

	//deprecated?
	ListDigraph* getListDigraph();
	ListDigraph::NodeMap<Molecule*>* getNodeMap();
	ListDigraph::ArcMap<Interaction*>* getArcMap();
	
	//mutation methods
	void newBasic();
	void forwardRateChange();
	void reverseRateChange();
	void degradationRateChange();
	DNA* histoneMod();
	void newComplex();	
	void newPromoter();
	void newPTM();	

private:

	float max_rate;
	float min_rate;

	int maxComp;
	int maxBasic;
	int maxProm;

	//graph structure
	ListDigraph* derivs;
	ListDigraph::NodeMap<Molecule*>* molecules;
	ListDigraph::ArcMap<Interaction*>* interactions;

	// molecule lists
	vector<Molecule*>* MoleculeList;
	vector<Protein*>* ProteinList;
	vector<mRNA*>* mRNAList;
	vector<DNA*>* DNAList;
	vector<Complex*>* ComplexList;
	vector<PTMProtein*>* PTMList;

	// interaction lists
	vector<Interaction*>* InteractionList;
	vector<Transcription*>* TranscriptionList;
	vector<Translation*>* TranslationList;
	vector<Degradation*>* DegradationList;
	vector<ForwardComplexation*>* ForwardComplexationList;
	vector<ReverseComplexation*>* ReverseComplexationList;
	vector<ForwardPTM*>* ForwardPTMList;
	vector<ReversePTM*>* ReversePTMList;
	vector<PromoterBind*>* PromoterBindList;

	//null node
	ListDigraph::Node nullnode;

	//helper method
	float getEffect(ListDigraph::Node, ListDigraph::Arc, int, float);
	//utility method
	ListDigraph::Node add(Molecule*);
	ListDigraph::Arc add(Interaction*, ListDigraph::Node, ListDigraph::Node);
	int count;
};


#endif
