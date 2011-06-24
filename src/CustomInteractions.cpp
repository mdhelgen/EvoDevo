/**
 * Implementation file for Custom Interactions.
 *
 * Interactions may overload the virtual method getEffect() to create a custom effect between Molecules
 *
 */

#include "CustomInteractions.h"

#include "ExternTrace.h"

Transcription::Transcription(){
	name="txn";
}

Transcription::~Transcription(){

}



float Transcription::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a, int rkIter, float rkStep){	


	t.trace("efct","Original Node value: %f\n", (*m)[a]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue());

	Molecule* thisMol = (*m)[a];
	Molecule* oppositeMol = (*m)[g->oppositeNode(a, g->arcFromId(arcID))];

	

	if(g->source(g->arcFromId(arcID)) == a)
		return 0;
	else
	{
		int prom_id = ((DNA*)oppositeMol)->promoterId;
		if(prom_id == -1)
			return oppositeMol->rkApprox(rkIter, rkStep) * rate;
		PromoterBind* pb = (PromoterBind*)(*i)[g->arcFromId(prom_id)];
		float f = pb->kf;
		float r = pb->kr;
		int h = ((DNA*)oppositeMol)->hill;
		return (1/(1+pow(f/r,h))) * oppositeMol->rkApprox(rkIter, rkStep) * rate;
	}
}


Degradation::Degradation(){
	name="deg";
}
Degradation::~Degradation(){}


float Degradation::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a, int rkIter, float rkStep){	
	
	t.trace("efct","Original Node value: %f\n", (*m)[a]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue());

	Molecule* thisMol = (*m)[a];
	Molecule* oppositeMol = (*m)[g->oppositeNode(a, g->arcFromId(arcID))];

	if(g->source(g->arcFromId(arcID)) == a)
		return -1 * thisMol->rkApprox(rkIter, rkStep) * rate;
	else
		return 0;

}

Translation::Translation(){

	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type Translation\n");
	t.trace("mloc","Interaction at location %u\n", (unsigned int)this);
	
	name="tsln";	

	t.trace("init","New Interaction created\n");
}

Translation::~Translation(){}


ForwardComplexation::ForwardComplexation(int n1, int n2){
	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type Complexation\n");
	t.trace("mloc","Interaction at location %u\n", (unsigned int) this);
	name="f_cmplx";

	firstNodeID = n1;
	secondNodeID = n2;
	
	t.trace("init","New Interaction created\n");


}
ForwardComplexation::~ForwardComplexation(){}

float ForwardComplexation::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a, int rkIter, float rkStep){	
	
	t.trace("efct","Original Node value: %f\n", (*m)[a]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue());

	Molecule* thisMol = (*m)[a];
	Molecule* oppositeMol = (*m)[g->oppositeNode(a, g->arcFromId(arcID))];

	Molecule* compMol1 = (*m)[g->nodeFromId(firstNodeID)];
	Molecule* compMol2 = (*m)[g->nodeFromId(secondNodeID)];


	//effect on source node
	if(g->source(g->arcFromId(arcID)) == a)
		return -1 * rate * compMol1->rkApprox(rkIter, rkStep) * compMol2->rkApprox(rkIter, rkStep);
	//effect on target node
	else
		return .5 * rate * compMol1->rkApprox(rkIter, rkStep) * compMol2->rkApprox(rkIter, rkStep);

}

ReverseComplexation::ReverseComplexation(int n1, int n2){
	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type Complexation\n");
	t.trace("mloc","Interaction at location %u\n", (unsigned int) this);

	name="r_cmplx";
	firstNodeID = n1;
	secondNodeID = n2;
	
	t.trace("init","New Interaction created\n");


}
ReverseComplexation::~ReverseComplexation(){}

float ReverseComplexation::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a, int rkIter, float rkStep){	
	
	t.trace("efct","Original Node value: %f\n", (*m)[a]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue());

	
	Molecule* oppositeMol = (*m)[g->oppositeNode(a, g->arcFromId(arcID))];

	Molecule* compMol1 = (*m)[g->nodeFromId(firstNodeID)];
	Molecule* compMol2 = (*m)[g->nodeFromId(secondNodeID)];


	//effect on source node
	if(g->source(g->arcFromId(arcID)) == a)
		return -1 * .5 * rate * (*m)[a]->rkApprox(rkIter, rkStep);
	//effect on target node
	else
		return  rate * oppositeMol->rkApprox(rkIter, rkStep); 

}

ForwardPTM::ForwardPTM(){

	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type ForwardPTM\n");
	t.trace("mloc","Interaction at location %u\n", (unsigned int)this);
	
	name="f_ptm";	

	t.trace("init","New Interaction created\n");
}

ForwardPTM::~ForwardPTM(){}
ReversePTM::ReversePTM(){

	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type ReversePTM\n");
	t.trace("mloc","Interaction at location %u\n", (unsigned int)this);
	
	name="r_ptm";	

	t.trace("init","New Interaction created\n");
}

ReversePTM::~ReversePTM(){}

TestInt::TestInt(){

	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type TestInt\n");
	t.trace("mloc","Interaction at location %u\n", (unsigned int)this);
	
	name="test";	
	t.trace("init","New Interaction created\n");
}

TestInt::~TestInt(){


}


float TestInt::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a, int rkIter, float rkStep){	
	
	
	t.trace("efct","Original Node value: %f\n", (*m)[a]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue());


	Molecule* oppositeMol = (*m)[g->oppositeNode(a, g->arcFromId(arcID))];

	if(g->source(g->arcFromId(arcID)) == a)
		return 0;
	else
		return oppositeMol->rkApprox(rkIter, rkStep) * rate;


}

PromoterBind::PromoterBind(float fwdRate, float revRate){
	name="pro";
	kf = fwdRate;
	kr = revRate;
	rate = kf - kr;
}
PromoterBind::~PromoterBind(){}


float PromoterBind::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node a, int rkIter, float rkStep){	
	
	t.trace("efct","Original Node value: %f\n", (*m)[a]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == a) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(a, g->arcFromId(arcID))]->getValue());

	Molecule* thisMol = (*m)[a];
	Molecule* oppositeMol = (*m)[g->oppositeNode(a, g->arcFromId(arcID))];

	if(g->source(g->arcFromId(arcID)) == a)
		return -1 * oppositeMol->rkApprox(rkIter, rkStep) * (kf-kr);
	else
		return 0;

}
