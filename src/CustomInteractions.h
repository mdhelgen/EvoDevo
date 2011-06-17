#ifndef CUSTOM_INTERACTIONS_H_
#define CUSTOM_INTERACTIONS_H_

#include <cmath>

#include "Interaction.h"

class Test : public Interaction{

public:
	Test();
	~Test();
	
	virtual float getEffect(ListDigraph*, ListDigraph::NodeMap<Molecule*>*, ListDigraph::ArcMap<Interaction*>*, ListDigraph::Node, int);


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
