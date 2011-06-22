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

	Molecule::Molecule();
	t.trace("init","Molecule %u type:DNA\n", (unsigned int) this);	
	kf = 0;
	kr = .35;
	hill = 2;

	longName = "DNA";
	shortName = "d";
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

	return (1/(1+pow(kf/kr,hill)));

}

float DNA::rkApprox(int rkstep, float step)
{
	return DNA::getValue();
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
	Molecule::Molecule();
	t.trace("init","Molecule %u type:mRNA\n", (unsigned int) this);	
	longName = "mRNA";
	shortName = "m";
}
mRNA::~mRNA(){
}

Protein::Protein(){
	Molecule::Molecule();
	t.trace("init","Molecule %u type:Protein\n", (unsigned int) this);	
	longName = "Protein";
	shortName = "p";
}
Protein::~Protein(){}

Complex::Complex(int n1, int n2){
	
	Molecule::Molecule();
	t.trace("init","Molecule %u type:Complex\n", (unsigned int) this);	
	
	longName="Complex";
	shortName="c";
	id1 = n1;
	id2 = n2;
}
Complex::~Complex(){

}

