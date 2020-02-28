//STRUCTURE MODULE - SOURCE FILE//
//////////////////////////////////


//C++ libraries
#include <vector>
#include <math.h>
#include <iostream>
#include <cstdlib>
//Code defined libraries

namespace spin {

//namespace variables

std::vector<double> sx; 
std::vector<double> sy; 
std::vector<double> sz; 

std::vector<double> lambda;

std::vector<double> global={0,0,0};

//namespace functions
void init_config_f(int n_atoms, 
		   std::vector<double> &x0,
		   std::vector<double> &y0,
		   std::vector<double> &z0, 
		   std::vector<double> &x,
		   std::vector<double> &y,
		   std::vector<double> &z,
		   std::vector<int> &mat_id)  
{
 x.reserve(n_atoms);
 y.reserve(n_atoms);
 z.reserve(n_atoms);				  
 
 for(int count_atom=0; count_atom<n_atoms; count_atom++)
 {
  int idx = mat_id[count_atom];
  x[count_atom] = x0[idx];
  y[count_atom] = y0[idx];
  z[count_atom] = z0[idx];  	 
 }   
}

void init_constr_f(double lx,
		   double ly,
		   double lz,
		   double l1,
		   std::vector<double> &constraint)
{
// std::cout<<"Initialise the Lagrange parameters(globally)"<<std::endl;
 constraint.reserve(4);
 constraint[0] = lx;
 constraint[1] = ly;  //setup the vectorial Lagrange constraint
 constraint[2] = lz;

 constraint[3] = l1; //setup the scalar Lagrange constraint
    
}

}


namespace atom{


//---namespace variables
 
 std::vector<double> pos_x;
 std::vector<double> pos_y;
 std::vector<double> pos_z;

 std::vector<int> material_id;

//---namespace functions

void generate_f(int n_atoms,
		int nx_atoms,
		int ny_atoms,
		int nz_atoms,
		std::vector<int> &nmat_atoms,
		std::vector<double> &x,
		std::vector<double> &y,
		std::vector<double> &z,
		std::vector<int> &mat_id)
{
  x.reserve(n_atoms);
  y.reserve(n_atoms);
  z.reserve(n_atoms);
  mat_id.reserve(n_atoms);
	
  int count_atom = 0;
  int count_material = 0; 
  int base = 0;

  //---set up the real atom's coordinates
  //---and assign each atom's material id
  for(int idz=0; idz<nz_atoms; idz++)
  {
   for(int idy=0; idy<ny_atoms; idy++)
   {
    for(int idx=0; idx<nx_atoms; idx++)
    {
     if(nmat_atoms[count_material]==0){count_material++;}
     x[count_atom] = idx;
     y[count_atom] = idy;
     z[count_atom] = idz;
     mat_id[count_atom] = count_material;
     //---check if all atoms in a material are counted
     if((count_atom-base) == nmat_atoms[count_material]-1)//if I looped over all atoms in one material then...
     {
      base = base + nmat_atoms[count_material];//..add all those atom to the base variable
      count_material++;//increment material id
     }
     count_atom++;
    }
   }
  }
}
}
