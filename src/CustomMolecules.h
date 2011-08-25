#ifndef CUSTOM_MOLECULES_H_
#define CUSTOM_MOLECULES_H_

#include <cmath>
#include <map>
#include "Molecule.h"

class PTMProtein : public Molecule{

public:

	PTMProtein();
	PTMProtein(PTMProtein* );
	
	~PTMProtein();
	char* getLongName();	
	void addRandPTM(int);
	
	void setPTMCount(int, int);
};
class DNA : public Molecule{

public:
	DNA();
	~DNA();
	
	float getValue();
	float rkApprox(int, float);
	void setHistoneModValue(float);
	virtual	void setValue(float);
	int promoterId;
	int hill;
private:
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

class Complex : public Protein{

public:
	Complex(int, int);
	~Complex();
	int getComponentId(int);
private:

	int id1;
	int id2;

};
#endif
