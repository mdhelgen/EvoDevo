/**
 * CustomMolecules implementation file.
 *
 * Custom Molecules allow modification of the default behavior of molecules.
 *
 * Each molecule must be defined in the CustomMolecules.h header file.
 *
 */
#include "CustomMolecules.h"

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

Null::Null(){
	longName = "NullNode";
	shortName = "n";
}
Null::~Null(){}
float Null::getValue(){
	return 0;
}

mRNA::mRNA(){
	longName = "mRNA";
	shortName = "m";
}
mRNA::~mRNA(){
}

Protein::Protein(){
	longName = "Protein";
	shortName = "p";
}
Protein::~Protein(){}

Complex::Complex(int n1, int n2){
	longName="Complex";
	shortName="c";
	id1 = n1;
	id2 = n2;
}
Complex::~Complex(){

}

