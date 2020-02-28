//ENERGY MODULE - SOURCE FILE//
//////////////////////////////


//C++ libraries
#include <vector>
#include <math.h>
#include <iostream>

//Code defined libraries
#include "input.h"
#include "energy.h"
#include "utils.h"
#include "structure.h"



//DEFINING THE MODULES' FUNCTIONS & INITIALISING THE EXTERNAL VARIABLES

namespace energy{

//namespace variables
double zeeman;
double exchange;
double uniaxial_anisotropy;
double lagrange_term;
double total;
double real;

//namespace functions
/*
void print_f()
{

    std::cout<<"Anisotropy energy: "<<energy::uniaxial_anisotropy<<" J"<<std::endl;
    std::cout<<"Exchange energy: " <<energy::exchange<<" J"<<std::endl;
    std::cout<<"Lagrange energy: " <<energy::lagrangian_term<<" J"<<std::endl;
    std::cout<<"Zeeman energy: " <<energy::zeeman<<" J"<<std::endl;
    std::cout<<"Total energy(+ Lagrange term): "<<energy::total<<" J"<<"\n";
    std::cout<<"Real energy(- Lagrange term): "<<energy::real_energy<<" J"<<"\n\n";

}

*/
void compute_f()
{
 energy::zeeman = energy::zeeman_f(input::no_of_atoms,
		                   input::hx,
		                   input::hy,
               			   input::hz,
              			   input::H,
		                   input::spin_moment_list,
		                   spin::sx,
               			   spin::sy,
				   spin::sz,
			           atom::material_id);
 
 
 
 //compute the exchange energy
 energy::exchange = energy::exchange_f(input::no_of_atoms,
		                       spin::sx,
		                       spin::sy,
		                       spin::sz,
		                       atom::material_id,
		                       utils::i_list,
                  		       utils::start,
                  		       utils::end,
                  		       utils::Jmatrix);
 


 //compute the uniax. anisotropy energy
 energy::uniaxial_anisotropy = energy::uniaxial_anisotropy_f(input::no_of_atoms,
                             				     input::anis_const_list,
                             			             input::anis_dir_x,
                             				     input::anis_dir_y,
				                             input::anis_dir_z,
			                                     spin::sx,
                            				     spin::sy,
				                             spin::sz,
                             			             atom::material_id);
 
 //compute the lagrange term energy
 energy::lagrange_term = energy::lagrange_term_f(input::no_of_atoms,
				                 input::v0_constraint_z,
			                         spin::sx,
			                         spin::sy,
				                 spin::sz,
		                                 spin::lambda,
                         	                 spin::global);
 energy::total = energy::total_f(energy::zeeman,
		 	         energy::exchange,
		  		 energy::uniaxial_anisotropy,
		  		 energy::lagrange_term);

  energy::real  = (energy::total - energy::lagrange_term);
}

double zeeman_f(int n_atoms,
		int field_x,
		int field_y,
		int field_z,
		int field_amp,
		std::vector<double> &spin_mom,	
		std::vector<double> &sx,
		std::vector<double> &sy,
		std::vector<double> &sz,
		std::vector<int> &mat_id
		)
{
// std::cout<<"Computing the Zeeman energy .."<<std::endl;
 double en_zeeman = 0.0;
 for(int count_atom = 0; count_atom<n_atoms; count_atom++)
 {
  int idx = mat_id[count_atom];

  en_zeeman -= (sx[count_atom]*field_x +
		sy[count_atom]*field_y +
	        sz[count_atom]*field_z )*spin_mom[idx]*field_amp;
   
 }

    return en_zeeman;

}


double exchange_f(int n_atoms,
		  std::vector<double> &sx,
		  std::vector<double> &sy,
		  std::vector<double> &sz,
		  std::vector<int> &mat_id,
		  std::vector<int> &int_list,
 		  std::vector<int> &start,
		  std::vector<int> &end,
		  std::vector<std::vector<double> > &Jmat
		  )
{
 double exch_en = 0.0;
 for(int atom=0;atom<n_atoms;atom++)
 {
  for(int i=start[atom]; i<end[atom]; i++)
  {
   int neigh=int_list[i];//get neighbour
   int m_neigh=mat_id[neigh];//get neighour's material id
   int m_atom=mat_id[atom];//get atom's material id
   exch_en -= Jmat[m_atom][m_neigh]*(sx[atom]*sx[neigh] +
		  		     sy[atom]*sy[neigh] +
				     sz[atom]*sz[neigh]);
  }
 }
 //the method considers both i-->j & j-->i interactions, hence the normalisation
 exch_en *= 0.5;
 return exch_en;
}



double uniaxial_anisotropy_f(int n_atoms,
			     std::vector<double> &anis_cst,
			     std::vector<double> &anis_dir_x,
			     std::vector<double> &anis_dir_y,
			     std::vector<double> &anis_dir_z,
			     std::vector<double> &sx,
			     std::vector<double> &sy,
			     std::vector<double> &sz,
			     std::vector<int> &mat_id
			     )
{
//    std::cout<<"Computing the anisotropy energy .."<<std::endl; 
 double anisotropy_energy = 0.0;

 for(int count_atom = 0; count_atom<n_atoms; count_atom++)
 {
  int idx = mat_id[count_atom];
  anisotropy_energy -= anis_cst[idx]*(pow( sx[count_atom]*anis_dir_x[idx] +
                                           sy[count_atom]*anis_dir_y[idx] +
                                           sz[count_atom]*anis_dir_z[idx], 2.0));  
 }
    return anisotropy_energy;

}

double lagrange_term_f(int n_atoms,
			 double v0z_constraint,
			 std::vector<double> &sx,
			 std::vector<double> &sy,
			 std::vector<double> &sz,
			 std::vector<double> &constraint,
			 std::vector<double> &spin_global
			 )
{

    //std::cout<<"Computing the Lagrange-term energy .."<<std::endl;
    double lagrange_energy = 0.0;
    utils::compute_global_magnetisation(n_atoms,sx,sy,sz,spin_global); // call to function which computes the global magn. vector
    lagrange_energy -= n_atoms * (constraint[0] * spin_global[0]*( spin_global[2] - v0z_constraint) +
                                  constraint[1] * spin_global[1]*( spin_global[2] - v0z_constraint) +
                                  constraint[2] * spin_global[2]*( spin_global[2] - v0z_constraint) + //compute the vectorial component
				  constraint[3] * (spin_global[2] - v0z_constraint)); //add the scalar term;

    return lagrange_energy;
}

double total_f(double zeeman,
	       double exchange,
	       double anisotropy,
	       double lagrange)
{

 double en_total = 0.0;

 en_total = zeeman +
	    exchange +
	    anisotropy +
	    lagrange;
 return en_total;

}
 	

}
