/**
 * CustomMolecules implementation file.
 *
 * Custom Molecules allow modification of the default behavior of molecules.
 *
 * Each molecule must be defined in the CustomMolecules.h header file.
 *
 */
#include "CustomMolecules.h"
#include "ExternTrace.h"
/**
 * DNA::DNA()
 *
 * Default Constructor
 *
 * Derived from Molecule.
 *
 *
 */

extern int hillParam;
DNA::DNA(){

	t.trace("init","Molecule %p type:DNA\n", this);	
	promoterId = -1;
	currentConcentration = 1;
	hill = hillParam;
	histoneModValue = 1;
	currentDir = 0;
	prevDir = 0;
	longName = "DNA";
	shortName = "g";
}


/**
 * DNA::~DNA()
 * 
 * Default Destructor
 */
DNA::~DNA(){


}

/**
 * Overload of virtual method Molecule::getValue 
 *
 * 
 *
 * @return Goodwin term describing the probability that the DNA is available for transcription.
 */
float DNA::getValue(){


	
	return histoneModValue;
	/*
	if(kf == -1 && kr == -1)
		return histoneModValue * 1.0;
	
	return histoneModValue * (1/(1+pow(kf/kr,hill)));
	*/
}
float DNA::rkApprox(int rkstep, float step)
{
	return currentConcentration * DNA::getValue();
}
void DNA::setHistoneModValue(float newVal){
	histoneModValue = newVal;
}

void DNA::setValue(float newValue){

	newValue = 1;	
	initialConcentration = newValue;
	currentConcentration = newValue;
	rungeKuttaSolution.push_back(newValue);

}
void DNA::setG(float g){
	
	initialConcentration = g;
	currentConcentration = g;
	rungeKuttaSolution.push_back(g);
	reset();
}

NullNode::NullNode(){
	
	t.trace("init","Molecule %p type:NulNode\n", this);	
	longName = "NullNode";
	shortName = "n";
	currentDir = 0;
	currentConcentration = 0;
	prevDir = 0;
}
NullNode::~NullNode(){}
float NullNode::getValue(){
	return 0;
}

mRNA::mRNA(){
	
	t.trace("init","Molecule %p type:mRNA\n", this);	
	longName = "mRNA";
	shortName = "m";
	currentDir = 0;
	prevDir = 0;
	currentConcentration = 0;

}
mRNA::~mRNA(){
}

Protein::Protein(){

	t.trace("init","Molecule %p type:Protein\n", this);	
	longName = "Protein";
	shortName = "p";
	currentDir = 0;
	prevDir = 0;
	currentConcentration = 0;
	numChanges = 0;
}
Protein::~Protein(){}

Complex::Complex(int n1, int n2){

	t.trace("init","Molecule %p type:Complex\n", this);	
	
	currentDir = 0;
	prevDir = 0;

	longName="Complex";
	shortName="c";
	id1 = n1;
	id2 = n2;
	currentConcentration = 0;
}
Complex::~Complex(){

}


int Complex::getComponentId(int i){
	if (i==1)
		return id1;
	if (i==2)
		return id2;
	return 0;

}


PTMProtein::PTMProtein(){

	t.trace("init","Molecule %p type:PTM\n",  this);

	longName = "PTM";
	shortName = "ptm";

	currentDir = 0;
	prevDir = 0;

	PTMArray[0] = 0;
	PTMArray[1] = 0;
	PTMArray[2] = 0;
	PTMArray[3] = 0;

	wasPTM = 0;
}

PTMProtein::PTMProtein(PTMProtein* c )
{
	PTMArray[0] = c->PTMArray[0];
	PTMArray[1] = c->PTMArray[1];
	PTMArray[2] = c->PTMArray[2];
	PTMArray[3] = c->PTMArray[3];

}
void PTMProtein::addRandPTM(int i){

	PTMArray[i]++;

}

void PTMProtein::setPTMCount(int index, int count){
	
	PTMArray[index] = count;


}

char* PTMProtein::getLongName(){

        memset(buf, '\0', 80);
        sprintf(buf, "%s %d (%d,%d,%d,%d)", longName, moleculeID, PTMArray[0],PTMArray[1],PTMArray[2],PTMArray[3]);
        return buf;
}

PTMProtein::~PTMProtein(){
}



