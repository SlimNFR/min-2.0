//STRUCTURE MODULE - HEADER FILE//
//////////////////////////////////


#ifndef STRUCTURE_H
#define STRUCTURE_H


//C++ libraries
#include <vector>
#include <math.h>
#include <iostream>

//Code defined libraries
#include "structure.h"
//#include "input.h"
//#include "utils.h"


namespace spin {

//namespace variables

extern std::vector<double> sx; //saves the x component of spins
extern std::vector<double> sy; //saves the y component of spins
extern std::vector<double> sz; //saves the z component of spins

extern std::vector<double> lambda; //saves the Lagrange parameters

extern std::vector<double> global; //saves the global magnetisation


//namespace functions

//---function to generate initial spin configuration
void init_config_f(int n_atoms,
                   std::vector<double> &x0,
                   std::vector<double> &y0,
                   std::vector<double> &z0,
                   std::vector<double> &x,
                   std::vector<double> &y,
                   std::vector<double> &z,
		   std::vector<int> &mat_id);

//---function to initialise the lagrangian constraints
void init_constr_f(double lx,
                   double ly,
                   double lz,
                   double l1,
                   std::vector<double> &constraint);
}

namespace atom {

//---namespace variables
extern std::vector<double> pos_x; //saves the x coordinates of the spins
extern std::vector<double> pos_y; //saves the y coordinates of the spins
extern std::vector<double> pos_z; //saves the z coordinates of the spins

extern std::vector<int> material_id;//saves the material id for each spin

//---namespace functions

//---function to generate the spatial system
void generate_f(int n_atoms,
		int nx_atoms,
		int ny_atoms,
		int nz_atoms,
		std::vector<int> &nmat_atoms,
		std::vector<double> &x, 
                std::vector<double> &y,
                std::vector<double> &z,
                std::vector<int> &mat_id);
}

#endif // STRUCTURE_H
