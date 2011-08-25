/**
 * Implementation file for Custom Interactions.
 *
 * Interactions may overload the virtual method getEffect() to create a custom effect between Molecules
 *
 */

#include "CustomInteractions.h"

#include "ExternTrace.h"

/**
 * Transcription::Transcription()
 *
 * Derived Interaction class for Transcription Interactions
 *
 * DNA >==(Transcription)==> mRNA
 *
 * DNA is not consumed in this process
 */
Transcription::Transcription(){
	name="txn";
}

Transcription::~Transcription(){}

/**
 * float Transcription::getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float)
 * 
 * Overload of virtual method Interaction::getEffect
 *
 * Transcription has a DNA source and an mRNA target. DNA is not removed as mRNA is created. The amount of mRNA created may be repressed
 *  by a protein which is "bound" to the DNA.
 *
 * @param g Reference to the Digraph object. Can be used to make complex effects based on local configuration around the current node
 * @param m Reference to the NodeMap object. Can be used to get the Molecule objects from the graph contaniers.
 * @param i Reference to the ArcMap object. Can be used to get the Interaction objects from the graph containers.
 * @param n The Node being affected by the interaction. Note that this can be the source or target.
 * @param rkIter The current iteration of Runge-Kutta evaluation.
 * @param rkStep The stepsize of the Runge-Kutta evaluation.
 *
 * @return The effect (change) of this interaction has on the concentration of node n at the current iteration of runge-kutta
 */
float Transcription::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node n, int rkIter, float rkStep){	


	t.trace("efct","Original Node value: %f\n", (*m)[n]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == n) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(n, g->arcFromId(arcID))]->getValue());

	Molecule* thisMol = (*m)[n];
	Molecule* oppositeMol = (*m)[g->oppositeNode(n, g->arcFromId(arcID))];
	
	
	
	if(isSourceNode(g, n) == 1){
		
		int prom_id = ((DNA*)thisMol)->promoterId;
		if(prom_id == -1)
			return 0;

		PromoterBind* pb = (PromoterBind*)(*i)[g->arcFromId(prom_id)];
		Molecule* repressor = (*m)[g->source(g->arcFromId(prom_id))];
		float f = pb->kf;
		return -1 * f * oppositeMol->rkApprox(rkIter, rkStep) * repressor->rkApprox(rkIter, rkStep);
	}
	else if(isTargetNode(g, n) == 1)
	{
		return oppositeMol->rkApprox(rkIter, rkStep) * rate;
	/*	
		PromoterBind* pb = (PromoterBind*)(*i)[g->arcFromId(prom_id)];
		Molecule* repressor = (*m)[g->source(g->arcFromId(prom_id))];

		float f = pb->kf;
		float r = pb->kr;
		int h = ((DNA*)oppositeMol)->hill;
		t.trace("hill","f:%f r:%f h:%d value:%f\n",f,r,h,(1/(1+pow(f/r,h))));
		return (1/(1+(f/r)*pow(repressor->rkApprox(rkIter,rkStep),h))) * oppositeMol->rkApprox(rkIter, rkStep) * rate;
	*/
	}
	else{
		t.trace("error", "%s getEffect reached error case, not source or target (%p)\n", name, this);
		return 0;
	}	
	
}

/**
 * Degradation::Degradation() 
 * 
 * Derived Interaction class for Degradation Interactions
 * 
 * [any molecule] >==(degradation)==> NullNode
 *
 * The null node always has a concentration of 0. 
 */
Degradation::Degradation(){
	name="deg";
}
Degradation::~Degradation(){}

/**
 * float Degradation::getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float)
 * 
 * Overload of virtual method Interaction::getEffect
 *
 * @param g Reference to the Digraph object. Can be used to make complex effects based on local configuration around the current node
 * @param m Reference to the NodeMap object. Can be used to get the Molecule objects from the graph contaniers.
 * @param i Reference to the ArcMap object. Can be used to get the Interaction objects from the graph containers.
 * @param n The Node being affected by the interaction. Note that this can be the source or target.
 * @param rkIter The current iteration of Runge-Kutta evaluation.
 * @param rkStep The stepsize of the Runge-Kutta evaluation.
 *
 * @return The effect (change) of this interaction has on the concentration of node n at the current iteration of runge-kutta
 */
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

/**
 * Translation::Translation()
 *
 * Derived class for Translation Interactions.
 *
 * mRNA >==(Translation)==> Protein
 *
 * mRNA is not consumed in this process
 */
Translation::Translation(){

	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type Translation\n");
	t.trace("mloc","Interaction at location %p\n", this);
	
	name="tsln";	

	t.trace("init","New Interaction created\n");
}

Translation::~Translation(){}

/**
 * float Translation::getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float)
 * 
 * Overload of virtual method Interaction::getEffect
 *
 * @param g Reference to the Digraph object. Can be used to make complex effects based on local configuration around the current node
 * @param m Reference to the NodeMap object. Can be used to get the Molecule objects from the graph contaniers.
 * @param i Reference to the ArcMap object. Can be used to get the Interaction objects from the graph containers.
 * @param n The Node being affected by the interaction. Note that this can be the source or target.
 * @param rkIter The current iteration of Runge-Kutta evaluation.
 * @param rkStep The stepsize of the Runge-Kutta evaluation.
 *
 * @return The effect (change) of this interaction has on the concentration of node n at the current iteration of runge-kutta
 */
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

/**
 * ForwardComplexation::ForwardComplexation(int, int)
 * 
 * Derived Interaction class for ForwardComplexation Interactions.
 *
 * A Protein-Protein Complex is made up of two sub-proteins. A ForwardComplexation interaction exists between each protein and the complex.
 *
 * The NodeID's of the two proteins involved are saved in each ForwardComplexation interaction, because the effect of the interaction is based
 *  on the concentration of both proteins, and each interaction should have the same effect. The NodeID's allow a simple way to retreive these 
 *  proteins during ForwardComplexation::getEffect() 
 *
 * P1 >==(ForwardComplexation)==> 
 *                                 Complex
 * P2 >==(ForwardComplexation)==> 
 *
 */
ForwardComplexation::ForwardComplexation(){
	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type Complexation\n");
	t.trace("mloc","Interaction at location %p\n", this);
	name="f_cmplx";

	
	t.trace("init","New Interaction created\n");


}
ForwardComplexation::~ForwardComplexation(){}

/**
 * float ForwardComplexation::getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float)
 * 
 * Overload of virtual method Interaction::getEffect
 *
 * @param g Reference to the Digraph object. Can be used to make complex effects based on local configuration around the current node
 * @param m Reference to the NodeMap object. Can be used to get the Molecule objects from the graph contaniers.
 * @param i Reference to the ArcMap object. Can be used to get the Interaction objects from the graph containers.
 * @param n The Node being affected by the interaction. Note that this can be the source or target.
 * @param rkIter The current iteration of Runge-Kutta evaluation.
 * @param rkStep The stepsize of the Runge-Kutta evaluation.
 *
 * @return The effect (change) of this interaction has on the concentration of node n at the current iteration of runge-kutta
 */
float ForwardComplexation::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node n, int rkIter, float rkStep){	
	
	t.trace("efct","Original Node value: %f\n", (*m)[n]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == n) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(n, g->arcFromId(arcID))]->getValue());

	Molecule* thisMol = (*m)[n];
	Molecule* oppositeMol = (*m)[g->oppositeNode(n, g->arcFromId(arcID))];
	Molecule* pairMol = (*m)[g->source(g->arcFromId(pairArcID))];


	//effect on source node
	if(isSourceNode(g, n) == 1)
		return -1 * rate * thisMol->rkApprox(rkIter, rkStep) * pairMol->rkApprox(rkIter, rkStep);
	//effect on target node
	else if(isTargetNode(g, n) == 1)
		return .5 * rate * oppositeMol->rkApprox(rkIter, rkStep) * pairMol->rkApprox(rkIter, rkStep);
	else{
		t.trace("error", "%s getEffect reached error case, not source or target (%p)\n", name, this);
		return 0;	
	}

}

/**
 * void ForwardComplexation::setPairArcID(int i)
 *
 * Save the Lemon Graph library's arcID corresponding the other ForwardComplexation in the pair. (see the constructor)
 *
 * @param i the arcID of the arc holding the other interaction 
 *
 */
void ForwardComplexation::setPairArcID(int i){
	pairArcID = i;
}


/**
 * ReverseComplexation::ReverseComplexation(int, int)
 * 
 * Derived Interaction class for ReverseComplexation Interactions.
 *
 * A Protein-Protein Complex is made up of two sub-proteins. A ReverseComplexation interaction exists between the complex and each protein.
 *
 * The NodeID's of the two proteins involved are saved in each ReverseComplexation interaction, because the effect of the interaction is based
 *  on the concentration of both proteins, and each interaction should have the same effect. The NodeID's allow a simple way to retreive these 
 *  proteins during ReverseComplexation::getEffect() 
 *
 *       >==(ReverseComplexation)==> P1
 *  Complex                          
 *       >==(ReverseComplexation)==> P2
 *
 */
ReverseComplexation::ReverseComplexation(){
	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type Complexation\n");
	t.trace("mloc","Interaction at location %p\n", this);

	name="r_cmplx";
	
	t.trace("init","New Interaction created\n");


}
ReverseComplexation::~ReverseComplexation(){}

/**
 * float ReverseComplexation::getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float)
 * 
 * Overload of virtual method Interaction::getEffect
 *
 * @param g Reference to the Digraph object. Can be used to make complex effects based on local configuration around the current node
 * @param m Reference to the NodeMap object. Can be used to get the Molecule objects from the graph contaniers.
 * @param i Reference to the ArcMap object. Can be used to get the Interaction objects from the graph containers.
 * @param n The Node being affected by the interaction. Note that this can be the source or target.
 * @param rkIter The current iteration of Runge-Kutta evaluation.
 * @param rkStep The stepsize of the Runge-Kutta evaluation.
 *
 * @return The effect (change) of this interaction has on the concentration of node n at the current iteration of runge-kutta
 */
float ReverseComplexation::getEffect(ListDigraph* g, ListDigraph::NodeMap<Molecule*>* m, ListDigraph::ArcMap<Interaction*>* i, ListDigraph::Node n, int rkIter, float rkStep){	
	
	t.trace("efct","Original Node value: %f\n", (*m)[n]->getValue());
	t.trace("efct","Interaction Rate: %f\n", rate);
	t.trace("efct","Interaction Dir: %s\n", (g->source(g->arcFromId(arcID)) == n) ? "outgoing" : "incoming");
	t.trace("efct","Opposite Node value: %f\n", (*m)[g->oppositeNode(n, g->arcFromId(arcID))]->getValue());

	
	Molecule* oppositeMol = (*m)[g->oppositeNode(n, g->arcFromId(arcID))];
	Molecule* pairMol = (*m)[g->source(g->arcFromId(pairArcID))];

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

/**
 * void ReverseComplexation::setPairArcID(int i)
 *
 * Save the Lemon Graph library's arcID corresponding the other ReverseComplexation in the pair. (see the constructor)
 *
 * @param i the arcID of the arc holding the other interaction 
 *
 */
void ReverseComplexation::setPairArcID(int i){
	pairArcID = i;
}

/**
 * ForwardPTM::ForwardPTM() 
 * 
 * Derived class for ForwardPTM interactions.
 *
 * A ForwardPTM interaction converts a Protein into a Post-translationally modified protein.
 * 
 * Protein >==(ForwardPTM)==> PTM
 */
ForwardPTM::ForwardPTM(){

	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type ForwardPTM\n");
	t.trace("mloc","Interaction at location %p\n", this);
	
	name="f_ptm";	

	t.trace("init","New Interaction created\n");
}

ForwardPTM::~ForwardPTM(){}

/**
 * ReversePTM::ReversePTM() 
 * 
 * Derived class for ReversePTM interactions.
 *
 * A ReversePTM interaction converts a Post-translationally modified protein back into it's original form.
 * 
 * PTM >==(ReversePTM)==> Protein
 */
ReversePTM::ReversePTM(){

	t.trace("init","Creating new Interaction\n");
	t.trace("cust","Custom Interaction type ReversePTM\n");
	t.trace("mloc","Interaction at location %p\n", this);
	
	name="r_ptm";	

	t.trace("init","New Interaction created\n");
}

ReversePTM::~ReversePTM(){}

/**
 * PromoterBind::PromoterBind
 * 
 * Derived class for PromoterBind interactions.
 *
 * A PromoterBind is an interaction between a Protein and a DNA, where the Protein binds to the DNA. This sequesters the protein (decreases concentration)
 *  and effects the rate at which the DNA produces mRNA in its Transcription interaction. The rate at which the protein is sequestered is the net total of
 *  the forward - reverse rates. These rates affect the DNA as well. 
 *
 *  Protein >==(PromoterBind)==> DNA
 *
 * @param fwdRate The rate of binding to the DNA
 * @param revRate The rate of unbinding to the DNA
 */
PromoterBind::PromoterBind(float fwdRate, float revRate){
	name="pro";
	kf = fwdRate;
	kr = revRate;
	rate = kf - kr;
	
	// -1 = repression, 1 = activation
	// modified by setAsRepression() / setAsActivation()
	promoterType = 0;
}
PromoterBind::~PromoterBind(){}

/**
 * float PromoterBind::getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float)
 * 
 * Overload of virtual method Interaction::getEffect
 *
 * @param g Reference to the Digraph object. Can be used to make complex effects based on local configuration around the current node
 * @param m Reference to the NodeMap object. Can be used to get the Molecule objects from the graph contaniers.
 * @param i Reference to the ArcMap object. Can be used to get the Interaction objects from the graph containers.
 * @param n The Node being affected by the interaction. Note that this can be the source or target.
 * @param rkIter The current iteration of Runge-Kutta evaluation.
 * @param rkStep The stepsize of the Runge-Kutta evaluation.
 *
 * @return The effect (change) of this interaction has on the concentration of node n at the current iteration of runge-kutta
 */
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
		return kr * (1 - thisMol->rkApprox(rkIter, rkStep));
	else{
		t.trace("error", "%s getEffect reached error case, not source or target (%p)\n", name, this);
		return 0;
	}	

}
/**
 * Determines if the promoter is repressing the gene it is binding to.
 * 
 * @return 1 if the promoter is acting as a repressor, and 0 otherwise.
*/
int PromoterBind::isRepression(){
	
	//promoterType: -1 = repression, 1 = activation
	if(promoterType == -1)
		return 1;
	return 0;
}

/**
 * Determines if the promoter is activating the gene it is binding to. 
 *
 * @return 1 if the promoter is acting as an activator, and 0 otherwise.
*/
int PromoterBind::isActivation(){

	//promoterType: -1 = repression, 1 = activation
	if(promoterType == 1)
		return 1;
	return 0;
}

/**
 * Sets the promoter interaction as a repression type. 
 *
 * The concentration of the protein being bound will cause the DNA being bound to to translate at less than it's basal rate (See Transcription::getEffect)
*/
void PromoterBind::setAsRepression(){
	promoterType = -1;
	name = "rep";
}

/**
 * Sets the promoter interaction as an activation type. 
 *
 * The concentration of the protein being bound will cause the DNA being bound to to translate at more than it's basal rate (See Transcription::getEffect)
*/
void PromoterBind::setAsActivation(){
	name = "act";
	promoterType = 1;
}
