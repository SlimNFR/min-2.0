//OUTPUT MODULE - HEADER FILE//
///////////////////////////////////



#ifndef OUTPUT_H
#define OUTPUT_H

//C++ libraries
#include <vector>
#include <math.h>
#include <iostream>

//Code defined libraries

namespace output{

//namespace variables


//namespace functions

//function to print the spin # and the spin components
void spin(int n_atoms,
          std::vector<double>&x,
          std::vector<double>&y,
	  std::vector<double>&z);

//---function to print the atom #,real coordinates and material id
void atom(int n_atoms,
          std::vector<double> &x,
          std::vector<double> &y,
          std::vector<double> &z,
          std::vector<int> &mat_id);

//---function to output the various energy components in the system
void energy(double zeeman,
            double exchange,
            double uniaxial_anisotropy,
            double lagrange_term,
            double en_total,
            double en_real);

//---function to output the gradient terms
void gradient(int n_atoms,
              std::vector<double> &grad_zeeman_x,
              std::vector<double> &grad_zeeman_y,
              std::vector<double> &grad_zeeman_z,
              std::vector<double> &grad_exch_x,
              std::vector<double> &grad_exch_y,
              std::vector<double> &grad_exch_z,
              std::vector<double> &grad_anis_x,
              std::vector<double> &grad_anis_y,
              std::vector<double> &grad_anis_z,
              std::vector<double> &grad_lagrange_x,
              std::vector<double> &grad_lagrange_y,
              std::vector<double> &grad_lagrange_z,
	      std::vector<double> &grad_l_param,
	      std::vector<double> &grad_total_x,
              std::vector<double> &grad_total_y,
              std::vector<double> &grad_total_z,
              std::vector<double> &grad_real_x,
              std::vector<double> &grad_real_y,
              std::vector<double> &grad_real_z
              );

//---function to output the torque components plus its modulus
void torque(int n_atoms,
            std::vector<double> &torque_x,
            std::vector<double> &torque_y,
            std::vector<double> &torque_z,
            std::vector<double> &torque_mod);
//---function to ouput the interaction list, the matrix exchange and others..
void lists(int n_mat,
           std::vector<int>&int_list,
           std::vector<int>&start,
           std::vector<int>&end,
           std::vector<std::vector<double> >&Jmat);

}

#endif // INPUT_H
                                                                                                                                                                                  
