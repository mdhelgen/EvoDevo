#ifndef CUSTOM_INTERACTIONS_H_
#define CUSTOM_INTERACTIONS_H_

#include <cmath>

#include "Interaction.h"
#include <cstdio>
#include "lemon/list_graph.h"
using namespace lemon;

class DerivGraph;

class Test : public Interaction{

public:
	Test();
	~Test();
	
	virtual float getEffect(DerivGraph*, ListDigraph::Node, int, float);
	//virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int, float);


private:
	float kf;
	float kr;
	int hill;

};


class Txn : public Interaction{


};

class Tsln : public Interaction{


};

class Deg : public Interaction{

};


#endif
