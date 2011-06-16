#ifndef CUSTOM_MOLECULES_H_
#define CUSTOM_MOLECULES_H_

#include <cmath>

#include "Molecule.h"

class DNA : public Molecule{

public:
	DNA();
	~DNA();
	
	float getValue();


private:
	float kf;
	float kr;
	int hill;




};

#endif
