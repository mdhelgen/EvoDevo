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

class Null : public Molecule{
public:
	Null();
	~Null();

	virtual float getValue();

};

class mRNA : public Molecule{
public:
	mRNA();
	~mRNA();

};

#endif
