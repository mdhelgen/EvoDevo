/**
 * Molecule.h
 *
 * Data structure for a Molecule in the Cell.
 *
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_

class Molecule{

public:
	Molecule();
	~Molecule();

	virtual float getValue();

private:
	float initialConcentration;
	float currentConcentration;


};



#endif
