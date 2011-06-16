#ifndef CUSTOM_INTERACTIONS_H_
#define CUSTOM_INTERACTIONS_H_

#include <cmath>

#include "Interaction.h"

class Test : public Interaction{

public:
	Test();
	~Test();
	
	virtual float getEffect(ListDigraph* d, ListDigraph::NodeMap<int>* b, ListDigraph::ArcMap<Interaction*>* c, ListDigraph::Node a);


private:
	float kf;
	float kr;
	int hill;




};

#endif
