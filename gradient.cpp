//GRADIENT MODULE - SOURCE FILE//
////////////////////////////////

//C++ libraries
#include <vector>
#include <math.h>
#include <iostream>

//Code defined libraries
#include "gradient.h"
#include "structure.h"
#include "utils.h"
#include "input.h"

//DEFINING THE MODULES' FUNCTIONS & INITIALISING THE EXTERNAL VARIABLES

namespace gradient {

//namespace variables

std::vector<double> zeeman_x;
std::vector<double> zeeman_y;
std::vector<double> zeeman_z;

std::vector<double> exchange_x;
std::vector<double> exchange_y;
std::vector<double> exchange_z;

std::vector<double> anisotropy_x;
std::vector<double> anisotropy_y;
std::vector<double> anisotropy_z;

std::vector<double> lagrange_x;
std::vector<double> lagrange_y;
std::vector<double> lagrange_z;

std::vector<double> lambda_param; //array to store the gradient of the energy with respect to the lambda parameters (the constraints)

std::vector<double> total_x;
std::vector<double> total_y;
std::vector<double> total_z;

std::vector<double> real_x;
std::vector<double> real_y;
std::vector<double> real_z;

std::vector<double>torque_x; 
std::vector<double>torque_y; //<torque> = J (Joule)
std::vector<double>torque_z;

std::vector<double>torque_mod;

std::vector<double>spin_dot_gradient;


double magnitude_spin_term;
double magnitude_lambda_term;

//namespace functions


void zeeman_f(int n_atoms,
	      int field_x,
              int field_y,
              int field_z,
              int field_amp,
	      std::vector<double> &spin_mom,
	      std::vector<double> &grad_zeeman_x,
	      std::vector<double> &grad_zeeman_y,
	      std::vector<double> &grad_zeeman_z,
	      std::vector<int> &mat_id
	      )
{
  std::fill(grad_zeeman_x.begin(), grad_zeeman_x.end(),0.0);
  std::fill(grad_zeeman_y.begin(), grad_zeeman_y.end(),0.0);
  std::fill(grad_zeeman_z.begin(), grad_zeeman_z.end(),0.0);

 for(int count_atom = 0; count_atom<n_atoms; count_atom++)
 { 
  int index=mat_id[count_atom];
  grad_zeeman_x[count_atom] = - spin_mom[index]*field_amp*field_x;
  grad_zeeman_y[count_atom] = - spin_mom[index]*field_amp*field_y;
  grad_zeeman_z[count_atom] = - spin_mom[index]*field_amp*field_z;
 }
}

void exchange_f(int n_atoms,
		std::vector<double> &sx,
		std::vector<double> &sy,
		std::vector<double> &sz,
		std::vector<double> &grad_exch_x,
		std::vector<double> &grad_exch_y,
		std::vector<double> &grad_exch_z,
		std::vector<int> &mat_id,
		std::vector<int> &int_list,
		std::vector<int> &start,
		std::vector<int> &end,
		std::vector<std::vector<double> > &Jmat
		)
{

 std::fill(grad_exch_x.begin(), grad_exch_x.end(),0.0);
 std::fill(grad_exch_y.begin(), grad_exch_y.end(),0.0);
 std::fill(grad_exch_z.begin(), grad_exch_z.end(),0.0);


// grad_exch_y[0]=0;
// grad_exch_y[1]=0;
// std::cout<<"BEFORE"<<gradient::exchange_y[1]<<"\n";
// std::cout<<"AFTER"<<grad_exch_y[1]<<"\n";
 for(int atom=0;atom<n_atoms;atom++)
 {
  for(int i=start[atom]; i<end[atom]; i++)
  {
   int neigh=int_list[i];//get neighbour
   int m_neigh=mat_id[neigh];//get neighour's material id
   int m_atom=mat_id[atom];//get atom's material id

   //the 0.5 factor in the gradient is consistent with the Hamiltonian used throughout the code;
   //the exchange energy is calculated as follows: E_{ex} = -J\sum{i<j}s_i*s_j

   grad_exch_x[atom] -= 0.5*Jmat[m_atom][m_neigh]*sx[neigh];
   grad_exch_y[atom] -= 0.5*Jmat[m_atom][m_neigh]*sy[neigh];
 //  std::cout<<" GRAD_EXCH_Y"<<grad_exch_y[atom]<<" "<<sy[neigh]<<" "<<Jmat[m_atom][m_neigh]<<"\n";
   grad_exch_z[atom] -= 0.5*Jmat[m_atom][m_neigh]*sz[neigh];

  } 
 }
}

void uniaxial_anisotropy_f(int n_atoms,
			   std::vector<double> &anis_cst,
			   std::vector<double> &anis_dir_x,
			   std::vector<double> &anis_dir_y,
			   std::vector<double> &anis_dir_z,
			   std::vector<double> &grad_anis_x,
			   std::vector<double> &grad_anis_y,
			   std::vector<double> &grad_anis_z,
			   std::vector<double> &sx,
			   std::vector<double> &sy,
			   std::vector<double> &sz,
			   std::vector<int> &mat_id)
{
//  std::cout<<"Computing the anisotropy gradient: "<<"\n";

 std::fill(grad_anis_x.begin(), grad_anis_x.end(),0.0);
 std::fill(grad_anis_y.begin(), grad_anis_y.end(),0.0);
 std::fill(grad_anis_z.begin(), grad_anis_z.end(),0.0);

 for(int count_atom = 0;count_atom<n_atoms;count_atom++)
 {
 int index = mat_id[count_atom];
 grad_anis_x[count_atom] = -2*anis_cst[index]*(sx[count_atom]*anis_dir_x[index]+
	           		 	       sy[count_atom]*anis_dir_y[index]+
				   	       sz[count_atom]*anis_dir_z[index])*anis_dir_x[index];

 grad_anis_y[count_atom] = -2*anis_cst[index]*(sx[count_atom]*anis_dir_x[index]+
                                   	       sy[count_atom]*anis_dir_y[index]+
                                   	       sz[count_atom]*anis_dir_z[index])*anis_dir_y[index];

 grad_anis_z[count_atom] = -2*anis_cst[index]*(sx[count_atom]*anis_dir_x[index]+
                                   	       sy[count_atom]*anis_dir_y[index]+
                                               sz[count_atom]*anis_dir_z[index])*anis_dir_z[index];
 }
}



void constraints_wrt_spins_f(int n_atoms,
			     double v0z_constraint,
			     std::vector<double> &sx,
			     std::vector<double> &sy,
			     std::vector<double> &sz,
			     std::vector<double> &constraint,
			     std::vector<double> &global,
			     std::vector<double> &grad_lagrange_x,
			     std::vector<double> &grad_lagrange_y,
			     std::vector<double> &grad_lagrange_z)
{
 std::fill(grad_lagrange_x.begin(), grad_lagrange_x.end(),0.0);
 std::fill(grad_lagrange_y.begin(), grad_lagrange_y.end(),0.0);
 std::fill(grad_lagrange_z.begin(), grad_lagrange_z.end(),0.0);

 utils::compute_global_magnetisation(n_atoms,sx,sy,sz,global);
 for(int i = 0; i < n_atoms; i++)
 { 
  grad_lagrange_x[i] = +n_atoms*constraint[0]*(global[2]-v0z_constraint);
  grad_lagrange_y[i] = +n_atoms*constraint[1]*(global[2]-v0z_constraint);
  grad_lagrange_z[i] = +n_atoms*constraint[2]*(2*global[2] - v0z_constraint) - n_atoms*constraint[2];
 }
}

void constraints_wrt_lambda_f(int n_atoms,
			      double v0z_constraint,
			      std::vector<double> &sx,
			      std::vector<double> &sy,
			      std::vector<double> &sz,
			      std::vector<double> &global,
			      std::vector<double> &grad_l_param)
{
 std::fill(grad_l_param.begin(), grad_l_param.end(),0);
 utils::compute_global_magnetisation(n_atoms,sx,sy,sz,global);
 grad_l_param[0] = -n_atoms*global[0]*(global[2] - v0z_constraint);
 grad_l_param[1] = -n_atoms*global[1]*(global[2] - v0z_constraint);
 grad_l_param[2] = -n_atoms*global[2]*(global[2] - v0z_constraint);
 grad_l_param[3] = -n_atoms*(global[2] - v0z_constraint);//additional constraint to make the "v0 = - 0.5" case work
}

void torque_f(int n_atoms,
	      std::vector<double> &sx,
	      std::vector<double> &sy,
	      std::vector<double> &sz,
	      std::vector<double> &grad_total_x,
	      std::vector<double> &grad_total_y,
	      std::vector<double> &grad_total_z,
	      std::vector<double> &torque_x,
	      std::vector<double> &torque_y,
	      std::vector<double> &torque_z,
	      std::vector<double> &torque_mod 
	      )
{
 std::fill(torque_x.begin(), torque_x.end(),0.0);
 std::fill(torque_y.begin(), torque_y.end(),0.0);
 std::fill(torque_z.begin(), torque_z.end(),0.0);
 std::fill(torque_mod.begin(), torque_mod.end(),0.0);

 for(int i=0; i<n_atoms; i++)
 {
  torque_x[i] = sy[i]*grad_total_z[i] - sz[i]*grad_total_y[i]; 
  torque_y[i] = sz[i]*grad_total_x[i] - sx[i]*grad_total_z[i]; 
  torque_z[i] = sx[i]*grad_total_y[i] - sy[i]*grad_total_x[i];

  torque_mod[i]  = sqrt(torque_x[i]*torque_x[i] +
                        torque_y[i]*torque_y[i] +
    	 	        torque_z[i]*torque_z[i] );
//  torque_x[i] /=torque_mod[i];
//  torque_y[i] /=torque_mod[i];
//  torque_z[i] /=torque_mod[i];
 }
 //std::cout<<"TORQUE"<<torque_x[0]<<" "<<torque_y[0]<<" "<<torque_z[0]<<"\n";
 //std::cout<<"TORQUE"<<torque_x[1]<<" "<<torque_y[1]<<" "<<torque_z[1]<<"\n";
 //std::cout<<"TORQUE"<<torque_x[2]<<" "<<torque_y[2]<<" "<<torque_z[2]<<"\n";
}

void remove_spin_projection(int n_atoms,
			    std::vector<double> &sx,
			    std::vector<double> &sy,
			    std::vector<double> &sz,
			    std::vector<double> &grad_total_x,
			    std::vector<double> &grad_total_y,
			    std::vector<double> &grad_total_z,
			    std::vector<double> &s_dot_grad
			    )
{
 std::fill(s_dot_grad.begin(), s_dot_grad.end(),0);

  //std::cout<<"RMPROJgrad_total_x"<<grad_total_x[1]<<"\n";
  //std::cout<<"RMPROJgrad_total_y"<<grad_total_y[1]<<"\n";
  //std::cout<<"RMPROJgrad_total_z"<<grad_total_z[1]<<"\n";

 // std::cout<<"RMPROJ sx"<<sx[1]<<"\n";
 // std::cout<<"RMPROJ sy"<<sy[1]<<"\n";
 // std::cout<<"RMPROJ sz"<<sz[1]<<"\n";


 for(int i=0; i<n_atoms; i++)
 {
  s_dot_grad[i]=sx[i]*grad_total_x[i]+
  	        sy[i]*grad_total_y[i]+
	        sz[i]*grad_total_z[i];
 // std::cout<<"RMPROJ s_dot_grad"<<s_dot_grad[i]<<"\n";
 // std::cout<<"RMPROJ s_dot_grad*sx"<<s_dot_grad[i]*sx[i]<<"\n";
  //std::cout<<"RMPROJ s_dot_grad*sy"<<s_dot_grad[i]*sy[i]<<"\n";
  //std::cout<<"RMPROJ s_dot_grad*sz"<<s_dot_grad[i]*sz[i]<<"\n";



  grad_total_x[i]=grad_total_x[i]-(s_dot_grad[i]*sx[i]);
  grad_total_y[i] -= s_dot_grad[i]*sy[i];
  grad_total_z[i] -= s_dot_grad[i]*sz[i];
 }

//  std::cout<<"RMPROJgrad_total_x"<<grad_total_x[1]<<"\n";
//  std::cout<<"RMPROJgrad_total_y"<<grad_total_y[1]<<"\n";
//  std::cout<<"RMPROJgrad_total_z"<<grad_total_z[1]<<"\n\n";

}

void total_f(int n_atoms,
	     std::vector<double> &total_x,
	     std::vector<double> &total_y,
	     std::vector<double> &total_z,
	     std::vector<double> &real_x,
	     std::vector<double> &real_y,
	     std::vector<double> &real_z
	     )
{
 std::fill(total_x.begin(), total_x.end(),0.0);     
 std::fill(total_y.begin(), total_y.end(),0.0);
 std::fill(total_z.begin(), total_z.end(),0.0);
 std::fill(real_x.begin(), real_x.end(),0.0);
 std::fill(real_y.begin(), real_y.end(),0.0);
 std::fill(real_z.begin(), real_z.end(),0.0);

 // call gradient computing functions
 gradient::zeeman_f(input::no_of_atoms,
                    input::hx,
		    input::hy,
		    input::hz,
		    input::H,
		    input::spin_moment_list,
		    gradient::zeeman_x,
		    gradient::zeeman_y,
		    gradient::zeeman_z,
		    atom::material_id);

 gradient::exchange_f(input::no_of_atoms,
                      spin::sx,
                      spin::sy,
                      spin::sz,
                      gradient::exchange_x,
                      gradient::exchange_y,
                      gradient::exchange_z,
                      atom::material_id,
                      utils::i_list,
		      utils::start,
		      utils::end,
		      utils::Jmatrix);

 gradient::uniaxial_anisotropy_f(input::no_of_atoms,
                                 input::anis_const_list,
                                 input::anis_dir_x,
                                 input::anis_dir_y,
                                 input::anis_dir_z,
                                 gradient::anisotropy_x,
                                 gradient::anisotropy_y,
                                 gradient::anisotropy_z,
                                 spin::sx,
                                 spin::sy,
                                 spin::sz,
                                 atom::material_id);

 gradient::constraints_wrt_spins_f(input::no_of_atoms,
                                   input::v0_constraint_z,
                                   spin::sx,
                                   spin::sy,
                                   spin::sz,
                                   spin::lambda,
                                   spin::global,
                                   gradient::lagrange_x,
                                   gradient::lagrange_y,
                                   gradient::lagrange_z);

 gradient::constraints_wrt_lambda_f(input::no_of_atoms,
                                    input::v0_constraint_z,
     	                            spin::sx,
           	                    spin::sy,
                                    spin::sz,
                       	            spin::global,
                             	    gradient::lambda_param);

 //compute the total gradient
 for(int i = 0; i<n_atoms;i++)
 {
  gradient::total_x[i] = gradient::zeeman_x[i] +
	       	         gradient::exchange_x[i] +
                         gradient::anisotropy_x[i] +
                         gradient::lagrange_x[i];

  gradient::total_y[i] = gradient::zeeman_y[i] +
	         	 gradient::exchange_y[i] +
                         gradient::anisotropy_y[i] +
                         gradient::lagrange_y[i];

  gradient::total_z[i] = gradient::zeeman_z[i] +
                         gradient::exchange_z[i] +
                         gradient::anisotropy_z[i] +
                         gradient::lagrange_z[i];

  gradient::real_x[i] = gradient::zeeman_x[i] +
                        gradient::exchange_x[i] +
                        gradient::anisotropy_x[i];

  gradient::real_y[i] = gradient::exchange_y[i] +
                        gradient::anisotropy_y[i] +
                        gradient::zeeman_y[i];

  gradient::real_z[i] = gradient::exchange_z[i] +
                        gradient::anisotropy_z[i] +
                        gradient::zeeman_z[i];

 }
}
}
