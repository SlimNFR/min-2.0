//INPUT MODULE - SOURCE FILE//
					//////////////////////////////


//C++ libraries
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>

//Code defined libraries
#include "input.h"
#include "energy.h"
#include "utils.h"
#include "structure.h"



//DEFINING THE MODULES' FUNCTIONS & INITIALISING THE EXTERNAL VARIABLES

namespace input{
//------------SYSTEM PARAMETERS-------------//
int nx = 1;
int ny = 1;
int nz = 1;
int no_of_atoms = nx*ny*nz;
//------------------------------------------//


//------------MATERIAL PARAMETERS-----------//
int n_materials = 1; 
std::vector<int>natoms_material{1};
std::vector<double>exchange_list{1.0,1,2,1,2};
std::vector<double>spin_moment_list{0.0,0.0,0.0};
std::vector<double>anis_const_list{0.038,0.0,0.0};
std::vector<double> anis_dir_x{0,0,0};
std::vector<double> anis_dir_y{0,0,0};
std::vector<double> anis_dir_z{1,0,0};
double norm_factor = 1.0;
//------------------------------------------//


//------------APPLIED FIELD-----------------//
double hx = 0.0; 
double hy = 0.0;
double hz = 0.0;
double H  = 0.0;
//------------------------------------------//


//------------SPIN'S INITIAL ORIENTATION----//
std::vector<double>s0x{0.0,0,0};
std::vector<double>s0y{0.0,0,0};
std::vector<double>s0z{1.0,0,0};
//------------------------------------------//


//------------LAGRANGE PARAMETERS-----------//
double lambda_x = 0.2;
double lambda_y = 0;
double lambda_z = 0.02;
double l1 = 0.001;
//------------------------------------------//


//------------ORIENTATION CONSTRAINT--------//
double v0_constraint_x;
double v0_constraint_y;
double v0_constraint_z = 0.5;
//------------------------------------------//
}
