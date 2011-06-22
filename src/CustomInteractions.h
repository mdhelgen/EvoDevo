#ifndef CUSTOM_INTERACTIONS_H_
#define CUSTOM_INTERACTIONS_H_

class DerivGraph;




#include <cmath>

#include "Interaction.h"
#include <cstdio>
#include "lemon/list_graph.h"
using namespace lemon;


class Test : public Interaction{

public:
	Test();
	~Test();
	
//	virtual float getEffect(DerivGraph*, ListDigraph::Node, int, float);
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);


private:
	float kf;
	float kr;
	int hill;

};

class Transcription : public Interaction{
public:
	Transcription();
	~Transcription();
	//virtual float getEffect(DerivGraph*, ListDigraph::Node, int, float);
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);
};

class Degradation : public Interaction{
public:
	Degradation();
	~Degradation();
	//virtual float getEffect(DerivGraph*, ListDigraph::Node, int, float);
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);
};


class Translation : public Interaction{
public:
	Translation();
	~Translation();

};

class ForwardComplexation : public Interaction{
public:
	ForwardComplexation(int, int);
	~ForwardComplexation();
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);
	int firstNodeID;
	int secondNodeID;
};
class ReverseComplexation : public Interaction{
public:
	ReverseComplexation(int, int);
	~ReverseComplexation();
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);
	int firstNodeID;
	int secondNodeID;
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



#endif
