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
#include <typeinfo>
#include <cstring>
#include "MersenneTwister.h"

using namespace std;


class Molecule{

public:
	Molecule();
	virtual ~Molecule();

	virtual float getValue();
	void updateRkVal(int, float);
	void nextPoint(float);
	virtual	void setValue(float);
	void outputRK();
	float getrkVal(int);
	vector<float>* getRungeKuttaSolution();
	virtual float rkApprox(int, float);
	virtual	char* getShortName();
	virtual char* getLongName();
	void setID(int);
	int getID();
	void reset();
	int nodeID;
	int wasPTM;
	int n;

	int stoch_numMols;

	int getScore();
	int PTMArray[4];
	int getPTMCount(int, int);
	MTRand r;
	
	virtual int getPTMCount(int);



protected:	
	float initialConcentration;
	float currentConcentration;
	float rkVal[4];
	int numChanges;
	int prevDir;
	int currentDir;


	char buf[200];	
	const char* longName;
	const char* shortName;
	int moleculeID;
	vector<float> rungeKuttaSolution;
	
	vector<float> maxima;
	vector<float> minima;
};



#endif
