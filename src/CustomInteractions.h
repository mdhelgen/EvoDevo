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


class Txn : public Interaction{


};

class Tsln : public Interaction{


};

class Deg : public Interaction{

};


#endif
