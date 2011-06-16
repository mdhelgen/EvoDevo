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

	
	kf = .4;
	kr = .35;
	hill = 2;

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
