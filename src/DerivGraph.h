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

struct prop_entry{
	float pe_propensity;
	Interaction* pe_interaction;
};

class DerivGraph{

public:
	DerivGraph();
	~DerivGraph();
	
	MTRand r;
	
	void test();
	void rungeKuttaEvaluate(float, float);
	void gillespieEvaluate(float);

	void outputDotImage(const char*, int, int, int);
	void outputDataPlot(const char*, int, int, int, float);
        void outputDataCsv(const char*, int , int, int, float);
	void outputInteractionCsv(const char*, int, int, int);

	void gillespieOutputDataCsv(const char*, int, int, int, float);

	Molecule* getBestMolecule(int);

	void setLimits(int, int, int, int);
	void setKineticRateLimits(float, float);
	void setRungeKuttaEval(float, float);
	void setDefaultInitialConc(float);


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

	float minKineticRate;
	float maxKineticRate;

	float defaultInitialConcentration;

	int maxComp;
	int maxBasic;
	int maxProm;
	int maxPTM;

	float rkTimeStep;
	float rkTimeLimit;

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

	//stochastic
	vector<ListDigraph::Arc> Propensities;

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
