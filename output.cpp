//OUTPUT MODULE - SOURCE FILE/////
//////////////////////////////////


//C++ libraries
#include <vector>
#include <math.h>
#include <iostream>

//Code defined libraries
#include "gradient.h"
namespace output{


void spin(int n_atoms,
          std::vector<double>&x,
          std::vector<double>&y,
          std::vector<double>&z)
{
    //std::cout<<"Printing the spins' components:"<<"\n";
    std::cout<<"#|sx|sy|sz\n";
    for(int count_spin = 0; count_spin<n_atoms;count_spin++)
    {
        std::cout<<count_spin<<" "
                 <<x[count_spin]<<" "
                 <<y[count_spin]<<" "
                 <<z[count_spin]<<"\n";
    }
    std::cout<<"\n";

}


void atom(int n_atoms,
          std::vector<double>&x,
          std::vector<double>&y,
          std::vector<double>&z,
          std::vector<int>&mat_id)
{
 std::cout<<"#|x|y|z|mat\n";
 for(int count_atom = 0; count_atom<n_atoms; count_atom++)
 {
 std::cout<<count_atom<<" " //print atom #
          <<x[count_atom]<<" " //print x position
          <<y[count_atom]<<" " //print y position
          <<z[count_atom]<<" " //print z position
          <<mat_id[count_atom]<<"\n"; //print material id of atom
 }
 std::cout<<"\n";
}

void energy(double zeeman,
	    double exchange,
	    double uniaxial_anisotropy,
	    double lagrange_term,
	    double en_total,
	    double en_real)
{
    std::cout<<"Zeeman energy: " <<zeeman<<" J\n";
    std::cout<<"Exchange energy: " <<exchange<<" J\n";
    std::cout<<"Anisotropy energy: "<<uniaxial_anisotropy<<" J\n";
    std::cout<<"Lagrange energy: " <<lagrange_term<<" J\n";
    std::cout<<"Total energy(+ Lagrange term): "<<en_total<<" J\n";
    std::cout<<"Real energy(- Lagrange term): "<<en_real<<" J"<<"\n\n";

}

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
	      )
{
 std::cout<<"Zeeman gradient:\n";
 for(int i=0; i<n_atoms; i++)
 {
  std::cout<<grad_zeeman_x[i]<<" "
           <<grad_zeeman_y[i]<<" "
           <<grad_zeeman_z[i]<<" "<<"\n";
 }
    
 std::cout<<"Exchange gradient:\n";
 for(int i=0; i<n_atoms; i++)
 {
  std::cout<<grad_exch_x[i]<<" "
	   <<grad_exch_y[i]<<" "
           <<grad_exch_z[i]<<" "<<"\n";
 }

 std::cout<<"Anisotropy gradient:\n";
 for(int i=0; i<n_atoms; i++)
 {
  std::cout<<grad_anis_x[i]<<" "
           <<grad_anis_y[i]<<" "
           <<grad_anis_z[i]<<" "<<"\n";
 }

 std::cout<<"Lagrange-term gradient:\n";
 for(int i=0; i<n_atoms; i++)
 {
  std::cout<<grad_lagrange_x[i]<<" "
           <<grad_lagrange_y[i]<<" "
           <<grad_lagrange_z[i]<<" "<<"\n";
 }
 std::cout<<"Gradient wrt. to lambda params:\n";
 for(int i=0; i<grad_l_param.size();i++)
 {
  std::cout<<grad_l_param[i]<<" ";
 }
 std::cout<<"\nTotal gradient(+ Lagrange term):\n";
 for(int i=0; i<n_atoms;i++)
 {
  std::cout<<grad_total_x[i]<<" "
	   <<grad_total_y[i]<<" "
	   <<grad_total_z[i]<<"\n";
 }
 std::cout<<"Total gradient(- Lagrange term):\n";
 for(int i=0; i<n_atoms;i++)
 {
  std::cout<<grad_real_x[i]<<" "
           <<grad_real_y[i]<<" "
           <<grad_real_z[i]<<"\n";
 }
 

}

void torque(int n_atoms,
	    std::vector<double> &torque_x,
	    std::vector<double> &torque_y,
	    std::vector<double> &torque_z,
	    std::vector<double> &torque_mod)
{
 std::cout<<"Torque: x|y|z|modulus\n";
 for(int i=0; i<n_atoms;i++)
 {
  std::cout<<torque_x[i]<<" "
	   <<torque_y[i]<<" "
	   <<torque_z[i]<<" "
	   <<torque_mod[i]<<"\n"; 
 }

}

void lists(int n_mat,
	   std::vector<int>&int_list,
	   std::vector<int>&start,
	   std::vector<int>&end,
	   std::vector<std::vector<double> >&Jmat)
{
//print Jmat
 std::cout<<"Jmatrix\n";
 for(int i = 0; i<n_mat;i++)
 {       
  for(int j = 0; j<n_mat;j++)
  {
   std::cout<<Jmat[i][j]<<" ";
  }
   std::cout<<"\n";
 }
 std::cout<<"\n";
//print interaction list
 std::cout<<"interaction list: \n";
 for(int i=0;i<int_list.size();i++)
 {
  std::cout<<int_list[i]<<" ";
 }
  std::cout<<"\n\n";
//print start
 std::cout<<"start list: \n";
 for(int i=0;i<start.size();i++)
 {
  std::cout<<start[i]<<" ";
 }
  std::cout<<"\n\n";
//print end
 std::cout<<"end list: \n";
 for(int i=0;i<end.size();i++)
 {
  std::cout<<end[i]<<" ";
 }
}

}

