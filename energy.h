//ENERGY MODULE - HEADER FILE//
//////////////////////////////



#ifndef ENERGY_H
#define ENERGY_H



//C++ libraries
#include<vector>

//Code defined libraries


//namespace for the energy terms which need to be minimised
namespace energy {


//namespace variables
extern double zeeman;
extern double exchange; //variable to save the exchange energy
extern double uniaxial_anisotropy; //variable to save the anisotropy energy
extern double lagrange_term; //variable to save the energy corresponding to the lagrangian terms
extern double total; //variable to save the total energy of the system(+ the lagrange term)
extern double real; //variable to save the real energy of the system(- the Lagrange term)

//namespace functions
void compute_f(); //func which allows energy calculations

//func to compute the zeeman energy
double zeeman_f(int n_atoms,
                int field_x,
                int field_y,
                int field_z,
                int field_amp,
                std::vector<double> &spin_mom,
                std::vector<double> &sx,
                std::vector<double> &sy,
                std::vector<double> &sz,
	        std::vector<int> &mat_id);

//function to compute the exchange energy
double exchange_f(int n_atoms,
                  std::vector<double> &sx,
                  std::vector<double> &sy,
                  std::vector<double> &sz,
                  std::vector<int> &mat_id,
                  std::vector<int> &int_list,
                  std::vector<int> &start,
                  std::vector<int> &end,
                  std::vector<std::vector<double> > &Jmat);

//function to calculate the uniaxial anisotropy energy
double uniaxial_anisotropy_f(int natoms,
                             std::vector<double> &anis_cst,
                             std::vector<double> &anis_dir_x,
                             std::vector<double> &anis_dir_y,
                             std::vector<double> &anis_dir_z,
                             std::vector<double> &sx,
                             std::vector<double> &sy,
                             std::vector<double> &sz,
                             std::vector<int> &mat_id);

//function to compute the energy related to the lagrangian terms
double lagrange_term_f(int n_atoms,
                         double v0z_constraint,
                         std::vector<double> &sx,
                         std::vector<double> &sy,
                         std::vector<double> &sz,
                         std::vector<double> &constraint,
                         std::vector<double> &spin_global);

//function which computes the total energy
double total_f(double zeeman,
               double exchange,
               double anisotropy,
               double lagrange);

}


#endif // ENERGY_H
