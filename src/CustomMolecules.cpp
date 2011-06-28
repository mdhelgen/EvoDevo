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
DNA::DNA(){

	t.trace("init","\n");
	Molecule::Molecule();
	t.trace("init","Molecule %u type:DNA\n", (unsigned int) this);	
	promoterId = -1;
	hill = 2;
	histoneModValue = 1;

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
	return DNA::getValue();
}
void DNA::setHistoneModValue(float newVal){
	histoneModValue = newVal;
}

NullNode::NullNode(){
	Molecule::Molecule();
	
	t.trace("init","Molecule %u type:NulNode\n", (unsigned int) this);	
	longName = "NullNode";
	shortName = "n";
}
NullNode::~NullNode(){}
float NullNode::getValue(){
	return 0;
}

mRNA::mRNA(){
	t.trace("init","\n");
	Molecule::Molecule();
	t.trace("init","Molecule %u type:mRNA\n", (unsigned int) this);	
	longName = "mRNA";
	shortName = "m";
}
mRNA::~mRNA(){
}

Protein::Protein(){
	t.trace("init","\n");
	Molecule::Molecule();
	t.trace("init","Molecule %u type:Protein\n", (unsigned int) this);	
	longName = "Protein";
	shortName = "p";

}
Protein::~Protein(){}

Complex::Complex(int n1, int n2){

	t.trace("init","\n");
	Molecule::Molecule();
	t.trace("init","Molecule %u type:Complex\n", (unsigned int) this);	
	
	longName="Complex";
	shortName="c";
	id1 = n1;
	id2 = n2;
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

	t.trace("init","Molecule %u type:PTM\n", (unsigned int) this);

	longName = "PTM";
	shortName = "ptm";

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



char* PTMProtein::getLongName(){

        memset(buf, '\0', 80);
        sprintf(buf, "%s %d (%d,%d,%d,%d)", longName, moleculeID, PTMArray[0],PTMArray[1],PTMArray[2],PTMArray[3]);
        return buf;
}

PTMProtein::~PTMProtein(){
}



