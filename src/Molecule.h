/**
 * Molecule.h
 *
 * Data structure for a Molecule in the Cell.
 *
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <cstdio>
#include <vector>
#include "DerivGraph.h"

using namespace std;

class DerivGraph;

class Molecule{

public:
	Molecule();
	virtual ~Molecule();

	virtual float getValue();
	void updateRkVal(int, float);
	void nextPoint(float);
	void setValue(float);
	void outputRK();
	float getrkVal(int);
	float rkApprox(int, float);
private:
	float initialConcentration;
	float currentConcentration;
	float rkVal[4]; 

	vector<float> rungeKuttaSolution;
};



#endif
