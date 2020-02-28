////////////////////////////////
//INPUT MODULE - HEADER FILE////
////////////////////////////////
//----------------------------------- THIS HEADER FILE CONTAINS INPUT CONSTANTS -----------------------------------------//


#ifndef INPUT_H
#define INPUT_H

//C++ libraries
#include<cmath>
#include<fstream>


namespace input{
//------------SYSTEM PARAMETERS-------------//
extern int nx;//transfer to input file
extern int ny;//transfer to input file
extern int nz;//trasnfer to input file
extern int no_of_atoms;//make the code calculate it
//------------------------------------------//


//------------MATERIAL PARAMETERS-----------//
extern int n_materials; //transfer
extern std::vector<int>natoms_material;//transfer
extern std::vector<double>exchange_list;//transfer | containts the exchange per material as well as the interlayers exchange
extern std::vector<double>spin_moment_list;//transfer | set the spin magnetic moments of each layer
extern std::vector<double>anis_const_list;//transfer
extern std::vector<double> anis_dir_x;//transfer
extern std::vector<double> anis_dir_y;//transfer
extern std::vector<double> anis_dir_z;//trasnfer
extern double norm_factor; //normalising factor transfer
//------------------------------------------//


//------------APPLIED FIELD-----------------//
extern double hx; //transfer
extern double hy;//transfer
extern double hz;//transfer
extern double H;//amplitude of the field, transfer
//------------------------------------------//


//------------SPIN'S INITIAL ORIENTATION----//
extern std::vector<double>s0x;//transfer
extern std::vector<double>s0y;//transfer
extern std::vector<double>s0z;//transfer
//------------------------------------------//


//------------LAGRANGE PARAMETERS-----------//
extern double lambda_x;//transfer
extern double lambda_y;//trasnfer
extern double lambda_z;//transfer
extern double l1;//transfer
//------------------------------------------//


//------------ORIENTATION CONSTRAINT--------//
extern double v0_constraint_x;//transfer
extern double v0_constraint_y;//transfer
extern double v0_constraint_z;//transfer
//------------------------------------------//
}
#endif // INPUT_H
