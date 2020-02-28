//UTILS MODULE - HEADER FILE//
//////////////////////////////

#ifndef UTILS_H
#define UTILS_H


//C++ libraries
#include <vector>
#include <fstream>


//Code defined libraries
namespace utils {


//namespace variables

extern std::vector<std::vector<double> >Jmatrix;
extern std::vector<double>global_magnetisation; // array to save the global magnetisation vector's components
extern std::vector<double> sums; //array to save the sums of the spins' components on x,y,z and the sum of the x components squared + sum of the y components squared + sum of the z components squared
extern std::vector<int>i_list;//array to store the list of interactions(including pairs)
extern std::vector<int>start;//used in conjunction with the i_list array
extern std::vector<int>end;//idem above


extern std::vector<double>temp_x;
extern std::vector<double>temp_y;
extern std::vector<double>temp_z;
extern std::vector<double>temp_lambda;       //temp variables to store the current state
extern std::vector<double>temp_gradient;
extern double temp_energy;


//namespace functions

//function to generate the exchange matrix
void generate_exchange_mat(int n_materials,
                           std::vector<double>&exch_list,
                           std::vector<std::vector<double> >&Jmat);

//function to generate the list of interactions
void generate_interac_list(int n_atoms,
                           std::vector<double> &x,
                           std::vector<double> &y,
                           std::vector<double> &z,
                           std::vector<int> &int_list,
                           std::vector<int> &start,
                           std::vector<int> &end);

//function to compute the global magnetisation
void compute_global_magnetisation(int n_atoms,
                                  std::vector<double> &sx,
                                  std::vector<double> &sy,
                                  std::vector<double> &sz,
                                  std::vector<double> &global); //function to compute the global magnetisation

double compute_XYZROW_vec_magnitude(std::vector<double> &x,
                                    std::vector<double> &y,
                                    std::vector<double> &z);

double compute_XYZ_vec_magnitude(double &x,
                                 double &y,
                                 double &z);

double compute_ROW_vec_magnitude(std::vector<double> &v);


double max_val_array(std::vector<double> &vec,
	             double& max_val);//function to extract the maximum value in one array

void print_all(std::ofstream& f1, std::ofstream& f2,std::ofstream& f3,int step);

void save_state(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                std::vector<double>& lambda,
                std::vector<double>& gradient,
                double& energy);


double compute_modul(std::vector<double> vec);


void compute_magnetisation_sums(std::vector<double>x, std::vector<double> y, std::vector<double>z,
                                std::vector<double>& sums); //function to compute the sums of the spin components on x,y,z and the  (m_ix)^2 + (m_iy)^2 + (m_iz)^2  term


}



#endif // UTILS_H
