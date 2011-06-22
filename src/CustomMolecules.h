#ifndef CUSTOM_MOLECULES_H_
#define CUSTOM_MOLECULES_H_

#include <cmath>

#include "Molecule.h"


class DNA : public Molecule{

public:
	DNA();
	~DNA();
	
	float getValue();
	float rkApprox(int, float);
	void setHistoneModValue(float);
private:
	float kf;
	float kr;
	int hill;
	float histoneModValue;
};

class NullNode : public Molecule{
public:
	NullNode();
	~NullNode();

	virtual float getValue();

};

class mRNA : public Molecule{
public:
	mRNA();
	~mRNA();

};

class Protein : public Molecule{
public:
	Protein();
	~Protein();

};

class Complex : public Molecule{

public:
	Complex(int, int);
	~Complex();

private:

	int id1;
	int id2;

};
#endif
