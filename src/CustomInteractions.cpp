/**
/a
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



float Transcription::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node n, int rkIter, float rkStep){	


	t.trace("efct","Original Node value: %f\n", (*m)[n]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == n) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(n, g->arcFromId(arcID))]->getValue());

	Molecule* thisMol = (*m)[n];
	Molecule* oppositeMol = (*m)[g->oppositeNode(n, g->arcFromId(arcID))];

	
	if(isSourceNode(g, n) == 1)
		return 0;
	else if(isTargetNode(g, n) == 1)
	{
		int prom_id = ((DNA*)oppositeMol)->promoterId;
		if(prom_id == -1)
			return oppositeMol->rkApprox(rkIter, rkStep) * rate;
		PromoterBind* pb = (PromoterBind*)(*i)[g->arcFromId(prom_id)];
		Molecule* repressor = (*m)[g->source(g->arcFromId(prom_id))];

		float f = pb->kf;
		float r = pb->kr;
		int h = ((DNA*)oppositeMol)->hill;
		t.trace("hill","f:%f r:%f h:%d value:%f\n",f,r,h,(1/(1+pow(f/r,h))));
		return (1/(1+(f/r)*pow(repressor->rkApprox(rkIter,rkStep),h))) * oppositeMol->rkApprox(rkIter, rkStep) * rate;
	}
	else{
		t.trace("error", "%s getEffect reached error case, not source or target (%p)\n", name, this);
		return 0;
	}	
	
}


Degradation::Degradation(){
	name="deg";
}
Degradation::~Degradation(){}


float Degradation::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node n, int rkIter, float rkStep){	
	
	t.trace("efct","Original Node value: %f\n", (*m)[n]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == n) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(n, g->arcFromId(arcID))]->getValue());

	Molecule* thisMol = (*m)[n];
	Molecule* oppositeMol = (*m)[g->oppositeNode(n, g->arcFromId(arcID))];

	if(isSourceNode(g, n) == 1)	
		return -1 * thisMol->rkApprox(rkIter, rkStep) * rate;
	else if(isTargetNode(g, n) == 1)
		return 0;
	else{
		t.trace("error", "%s getEffect reached error case, not source or target (%p)\n", name, this);
		return 0;
	}	
}

Translation::Translation(){

	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type Translation\n");
	t.trace("mloc","Interaction at location %p\n", this);
	
	name="tsln";	

	t.trace("init","New Interaction created\n");
}

Translation::~Translation(){}

float Translation::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node n, int rkIter, float rkStep){	
	
	t.trace("efct","Original Node value: %f\n", (*m)[n]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == n) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(n, g->arcFromId(arcID))]->getValue());
	t.trace("efct","isSourceNode() == %d\n", isSourceNode(g,n));
	Molecule* thisMol = (*m)[n];
	Molecule* oppositeMol = (*m)[g->oppositeNode(n, g->arcFromId(arcID))];

	//if the effect is being calculated for the source node
	if(isSourceNode(g, n) == 1)
		return 0;
	//effect is being calculated for the target node
	else if(isTargetNode(g, n) == 1)
		return oppositeMol->rkApprox(rkIter, rkStep) * rate;

	else{ 
		t.trace("error", "%s getEffect reached error case, not source or target (%p)\n", name, this);
		return 0;
	}
}

ForwardComplexation::ForwardComplexation(int n1, int n2){
	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type Complexation\n");
	t.trace("mloc","Interaction at location %p\n", this);
	name="f_cmplx";

	firstNodeID = n1;
	secondNodeID = n2;
	
	t.trace("init","New Interaction created\n");


}
ForwardComplexation::~ForwardComplexation(){}

float ForwardComplexation::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node n, int rkIter, float rkStep){	
	
	t.trace("efct","Original Node value: %f\n", (*m)[n]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == n) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(n, g->arcFromId(arcID))]->getValue());

	Molecule* thisMol = (*m)[n];
	Molecule* oppositeMol = (*m)[g->oppositeNode(n, g->arcFromId(arcID))];

	Molecule* compMol1 = (*m)[g->nodeFromId(firstNodeID)];
	Molecule* compMol2 = (*m)[g->nodeFromId(secondNodeID)];


	//effect on source node
	if(isSourceNode(g, n) == 1)
		return -1 * rate * compMol1->rkApprox(rkIter, rkStep) * compMol2->rkApprox(rkIter, rkStep);
	//effect on target node
	else if(isTargetNode(g, n) == 1)
		return .5 * rate * compMol1->rkApprox(rkIter, rkStep) * compMol2->rkApprox(rkIter, rkStep);
	else{
		t.trace("error", "%s getEffect reached error case, not source or target (%p)\n", name, this);
		return 0;	
	}

}

void ForwardComplexation::setPairArcID(int i){
	pairArcID = i;
}

ReverseComplexation::ReverseComplexation(int n1, int n2){
	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type Complexation\n");
	t.trace("mloc","Interaction at location %p\n", this);

	name="r_cmplx";
	firstNodeID = n1;
	secondNodeID = n2;
	
	t.trace("init","New Interaction created\n");


}
ReverseComplexation::~ReverseComplexation(){}

float ReverseComplexation::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node n, int rkIter, float rkStep){	
	
	t.trace("efct","Original Node value: %f\n", (*m)[n]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == n) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(n, g->arcFromId(arcID))]->getValue());

	
	Molecule* oppositeMol = (*m)[g->oppositeNode(n, g->arcFromId(arcID))];

	Molecule* compMol1 = (*m)[g->nodeFromId(firstNodeID)];
	Molecule* compMol2 = (*m)[g->nodeFromId(secondNodeID)];


	//effect on source node
	if(isSourceNode(g, n) == 1)
		return -1 * .5 * rate * (*m)[n]->rkApprox(rkIter, rkStep);
	//effect on target node
	else if(isTargetNode(g, n) == 1)
		return  rate * oppositeMol->rkApprox(rkIter, rkStep); 
	else{
		t.trace("error", "%s getEffect reached error case, not source or target (%p)\n", name, this);
		return 0;	
	}

}

void ReverseComplexation::setPairArcID(int i){
	pairArcID = i;
}

ForwardPTM::ForwardPTM(){

	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type ForwardPTM\n");
	t.trace("mloc","Interaction at location %p\n", this);
	
	name="f_ptm";	

	t.trace("init","New Interaction created\n");
}

ForwardPTM::~ForwardPTM(){}
ReversePTM::ReversePTM(){

	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type ReversePTM\n");
	t.trace("mloc","Interaction at location %p\n", this);
	
	name="r_ptm";	

	t.trace("init","New Interaction created\n");
}

ReversePTM::~ReversePTM(){}

TestInt::TestInt(){

	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type TestInt\n");
	t.trace("mloc","Interaction at location %p\n", this);
	
	name="test";	
	t.trace("init","New Interaction created\n");
}

TestInt::~TestInt(){


}


float TestInt::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node n, int rkIter, float rkStep){	
	
	
	t.trace("efct","Original Node value: %f\n", (*m)[n]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == n) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(n, g->arcFromId(arcID))]->getValue());


	Molecule* oppositeMol = (*m)[g->oppositeNode(n, g->arcFromId(arcID))];

	if(isSourceNode(g, n) == 1)
		return 0;
	else if(isTargetNode(g, n) == 1)
		return oppositeMol->rkApprox(rkIter, rkStep) * rate;
	else{
		t.trace("error", "%s getEffect reached error case, not source or target (%p)\n", name, this);
		return 0;	
	}

}

PromoterBind::PromoterBind(float fwdRate, float revRate){
	name="pro";
	kf = fwdRate;
	kr = revRate;
	rate = kf - kr;
}
PromoterBind::~PromoterBind(){}


float PromoterBind::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node n, int rkIter, float rkStep){	
	
	t.trace("efct","Original Node value: %f\n", (*m)[n]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == n) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(n, g->arcFromId(arcID))]->getValue());

	Molecule* thisMol = (*m)[n];
	Molecule* oppositeMol = (*m)[g->oppositeNode(n, g->arcFromId(arcID))];

	if(isSourceNode(g, n) == 1)
		return -1 * oppositeMol->rkApprox(rkIter, rkStep) * (kf-kr);
	else if(isTargetNode(g, n) == 1)
		return 0;
	else{
		t.trace("error", "%s getEffect reached error case, not source or target (%p)\n", name, this);
		return 0;
	}	

}
