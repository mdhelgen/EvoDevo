#ifndef CUSTOM_MOLECULES_H_
#define CUSTOM_MOLECULES_H_

#include "Molecule.h"

class DNA : public Molecule{

public:
	DNA();
	~DNA();
	
	float getValue();


private:
	float hill;




};

#endif
