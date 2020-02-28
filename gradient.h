//GRADIENT MODULE - SOURCE FILE//
////////////////////////////////

#ifndef GRADIENT_H
#define GRADIENT_H

//C++ libraries
#include <vector>
#include <iostream>

//Code defined libraries

//namespace for the gradient terms
namespace gradient {

//namespace variables

extern std::vector<double> zeeman_x;
extern std::vector<double> zeeman_y; //arrays to store the gradient of the zeeman energy with respect to the spins' compontents on x,y,z
extern std::vector<double> zeeman_z;

extern std::vector<double> exchange_x; //
extern std::vector<double> exchange_y; // arrays to store the gradient of the exchange energy term with respect to the spins' components on x,y,z
extern std::vector<double> exchange_z; //

extern std::vector<double> anisotropy_x;
extern std::vector<double> anisotropy_y; //arrays to store the gradient of the anis.energy term with respect to the spins' components on x, y, z
extern std::vector<double> anisotropy_z;

extern std::vector<double> lagrange_x;
extern std::vector<double> lagrange_y; //arrays to store the gradient of the lagrangian term with respect to the spins' components on x,y,z
extern std::vector<double> lagrange_z;

extern std::vector<double> lambda_param; //array to store the gradient of the energy with respect to the lambda parameters (the constraints)

extern std::vector<double> total_x;
extern std::vector<double> total_y; //arrays to store the total gradient (+ Lagrange term)
extern std::vector<double> total_z;

extern std::vector<double> real_x;
extern std::vector<double> real_y; //arrays to store the real gradient (- Lagrange term)
extern std::vector<double> real_z;

extern std::vector<double> torque_x; //array to store the modulus of the torque acting upon each of the spins
extern std::vector<double> torque_y; //array to store the modulus of the torque acting upon each of the spins
extern std::vector<double> torque_z; //array to store the modulus of the torque acting upon each of the spins

extern std::vector<double> torque_mod; //array to store the modulus of the torque acting upon each of the spins

extern std::vector<double>spin_dot_gradient; //array to store the projection of the spin on the gradient's direction

extern double magnitude_spin_term; //array to store the modulus of the spin related gradient
extern double magnitude_lambda_term; //array to store the modulus of the lambda related gradient

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
              std::vector<int> &mat_id);

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
                std::vector<std::vector<double> > &Jmat);

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
                           std::vector<int> &mat_id);

void constraints_wrt_spins_f(int n_atoms,
                             double v0z_constraint,
                             std::vector<double> &sx,
                             std::vector<double> &sy,
                             std::vector<double> &sz,
                             std::vector<double> &constraint,
			     std::vector<double> &global,
                             std::vector<double> &grad_lagrange_x,
                             std::vector<double> &grad_lagrange_y,
                             std::vector<double> &grad_lagrange_z);

void constraints_wrt_lambda_f(int n_atoms,
                              double v0z_constraint,
                              std::vector<double> &sx,
                              std::vector<double> &sy,
                              std::vector<double> &sz,
                              std::vector<double> &global,
                              std::vector<double> &grad_l_param);

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
              std::vector<double> &torque_mod);

void remove_spin_projection(int n_atoms,
                            std::vector<double> &sx,
                            std::vector<double> &sy,
                            std::vector<double> &sz,
                            std::vector<double> &grad_total_x,
                            std::vector<double> &grad_total_y,
                            std::vector<double> &grad_total_z,
                            std::vector<double> &s_dot_grad);
void total_f(int n_atoms,
             std::vector<double> &total_x,
             std::vector<double> &total_y,
             std::vector<double> &total_z,
             std::vector<double> &real_x,
             std::vector<double> &real_y,
             std::vector<double> &real_z);

}




#endif // GRADIENT_H
