//UPDATE MODULE - HEADER FILE///
////////////////////////////////


#ifndef UPDATE_H
#define UPDATE_H

//C++ libraries
#include <vector>
#include <iostream>

//Code defined libraries


namespace update{

//namespace variables

extern std::vector<double> lambda_dir;
extern std::vector<double> spin_dir_x;
extern std::vector<double> spin_dir_y;
extern std::vector<double> spin_dir_z;

//namespace functions

void lambda_dir_f(std::vector<double> &lambda,
                  std::vector<double> &update_lambda,
                  std::vector<double> &grad_lambda_param,
                  double magnitude_lambda_grad,
                  double step);

void spins_dir_f(int n_atoms,
                 std::vector<double> &sx,
                 std::vector<double> &sy,
                 std::vector<double> &sz,
                 std::vector<double> &update_sx,
                 std::vector<double> &update_sy,
                 std::vector<double> &update_sz,
                 std::vector<double> &grad_tot_x,
                 std::vector<double> &grad_tot_y,
                 std::vector<double> &grad_tot_z,
                 double &magnitude_spin_grad,
                 double step);

}
#endif // UPDATE_H
