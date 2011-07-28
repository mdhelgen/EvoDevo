#ifndef CUSTOM_INTERACTIONS_H_
#define CUSTOM_INTERACTIONS_H_

class DerivGraph;




#include <cmath>

#include "Interaction.h"
#include <cstdio>
#include "lemon/list_graph.h"
using namespace lemon;



class Transcription : public Interaction{
 public:
	Transcription();
	~Transcription();
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);
};

class Degradation : public Interaction{
 public:
	Degradation();
	~Degradation();
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);
};


class Translation : public Interaction{
 public:
	Translation();
	~Translation();
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);
};

class ForwardComplexation : public Interaction{
 public:
	ForwardComplexation();
	~ForwardComplexation();
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);
	void setPairArcID(int);
	int pairArcID;
};
class ReverseComplexation : public Interaction{
 public:
	ReverseComplexation();
	~ReverseComplexation();
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);
	void setPairArcID(int);
	int pairArcID;

};

class ForwardPTM : public Interaction{
 public:
	ForwardPTM();
	~ForwardPTM();
};

class ReversePTM : public Interaction{
 public:
	ReversePTM();
	~ReversePTM();

};

class PromoterBind : public Interaction{
 public:
	PromoterBind(float, float);
	~PromoterBind();
	float kf;
	float kr;
	
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);
};



#endif
