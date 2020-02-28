//SOLVER MODULE - SOURCE FILE//
///////////////////////////////


//C++ libraries
#include <vector>
#include <iostream>
#include <fstream>

//Code defined libraries
#include "input.h"
#include "energy.h"
#include "gradient.h"
#include "utils.h"
#include "update.h"
#include "test.h"
//namespace for the gradient terms
namespace solver{

//namespace variables

//namespace functions

void sdescent(int n_atoms,
	      double v0z_constraint,
	      double &magnitude_spin_grad,
	      double &magnitude_lambda_grad,
	      std::vector<double> &grad_total_x,
	      std::vector<double> &grad_total_y,
	      std::vector<double> &grad_total_z,
	      std::vector<double> &grad_real_x,
	      std::vector<double> &grad_real_y,
	      std::vector<double> &grad_real_z,
	      std::vector<double> &grad_lambda_param,
	      std::vector<double> &spin_dot_grad,
	      std::vector<double> &torque_x,
	      std::vector<double> &torque_y,
	      std::vector<double> &torque_z,
	      std::vector<double> &torque_mod,
	      std::vector<double> &spin_dir_x,
	      std::vector<double> &spin_dir_y,
	      std::vector<double> &spin_dir_z,
	      std::vector<double> &lambda_dir,
	      std::vector<double> &sx,
	      std::vector<double> &sy,
	      std::vector<double> &sz,
	      std::vector<double> &global,
	      std::vector<double> &lambda)
{
    std::ofstream f1;
    std::ofstream f2;
    std::ofstream f3;


    f1.open("output.txt", std::ofstream::out);
    f2.open("energy_and_gradient.txt",std::ofstream::out);
    f3.open("all.txt",std::ofstream::out);


int iter_spin, iter_constr, total_iter = 0;

double alpha_spin=0.1;
double alpha_constr=1;
double tol_torque=1e-5;
double max_iter_spin=5e3;
double max_iter_constr=250;
double max_torque_x=0.0;
double max_torque_y=0.0;
double max_torque_z=0.0;
double max_torque_mod=0.0;
double vz_dif;

iter_constr=0;
do
{
 iter_spin=0;
 do
 {
  //std::cout<<"\n\n";
  iter_spin++;
  //std::cout<<iter_spin<<"\n\n";

 //energy::compute_f();

//  test::spin();
 //test::gradient();
  //std::cout<<"Computing gradient.."<<"\n";
  gradient::total_f(n_atoms,
                    grad_total_x,
                    grad_total_y,
                    grad_total_z,
                    grad_real_x,
                    grad_real_y,
                    grad_real_z);
  //std::cout<<"SOLVERgrad_total_x"<<grad_total_x[1]<<"\n";
  //std::cout<<"SOLVERgrad_total_y"<<grad_total_y[1]<<"\n";
  //std::cout<<"SOLVERgrad_total_z"<<grad_total_z[1]<<"\n";
 
// test::gradient();
  //std::cout<<"Removing spin projection.."<<"\n";
  gradient::remove_spin_projection(n_atoms,
                                   sx,sy,sz,
                                   grad_total_x,
                                   grad_total_y,
                                   grad_total_z,
                                   spin_dot_grad); 
  //std::cout<<"SOLVERgrad_total_x"<<grad_total_x[1]<<"\n";
  //std::cout<<"SOLVERgrad_total_y"<<grad_total_y[1]<<"\n";
  //std::cout<<"SOLVERgrad_total_z"<<grad_total_z[1]<<"\n";

 
  magnitude_spin_grad=utils::compute_XYZROW_vec_magnitude(grad_total_x,
   		 		  	                  grad_total_y,
						 	  grad_total_z);

// std::cout<<"MAGNITUDE SPIN GRAD"<<magnitude_spin_grad<<"\n";
// test::gradient();

  if(magnitude_spin_grad<0.01)magnitude_spin_grad=0.1;


  //std::cout<<"Computing torque.."<<"\n";
/*  gradient::torque_f(n_atoms,
                     sx,sy,sz,
                     grad_total_x,
                     grad_total_y,
                     grad_total_z,
                     torque_x,
                     torque_y,
                     torque_z,
                     torque_mod);

  utils::max_val_array(torque_x,max_torque_x);
  utils::max_val_array(torque_y,max_torque_y);
  utils::max_val_array(torque_z,max_torque_z);

  utils::max_val_array(torque_mod,max_torque_mod);
  */

  for(int i =0;i<n_atoms;i++)
  {
   double Heff_x = alpha_spin*grad_total_x[i]/magnitude_spin_grad;
   double Heff_y = alpha_spin*grad_total_y[i]/magnitude_spin_grad;
   double Heff_z = alpha_spin*grad_total_z[i]/magnitude_spin_grad;

   torque_x[i] = sy[i]*Heff_z - sz[i]*Heff_y;
   torque_y[i] = sz[i]*Heff_x - sx[i]*Heff_z;
   torque_z[i] = sx[i]*Heff_y - sy[i]*Heff_x;

   torque_mod[i] = sqrt(torque_x[i]*torque_x[i] +
		        torque_y[i]*torque_y[i] +
			torque_z[i]*torque_z[i]);
 
  }
  utils::max_val_array(torque_x,max_torque_x);
  utils::max_val_array(torque_y,max_torque_y);
  utils::max_val_array(torque_z,max_torque_z); 
  
  utils::max_val_array(torque_mod,max_torque_mod);
  utils::print_all(f1,f2,f3,iter_spin+total_iter);

  //std::cout<<"Update spins.."<<"\n"; 
  update::spins_dir_f(n_atoms,
                      sx,sy,sz,
                      spin_dir_x,
		      spin_dir_y,
		      spin_dir_z,
		      grad_total_x,
 		      grad_total_y,
 		      grad_total_z,
		      magnitude_spin_grad,
		      alpha_spin);
 
  //test::gradient();
  /*test::spin();
  std::cout<<"lx"<<lambda[0]<<"\n";
  std::cout<<"ly"<<lambda[1]<<"\n";
  std::cout<<"lz"<<lambda[2]<<"\n";
  std::cout<<"lz"<<lambda[3]<<"\n";
  test::gradient();
  test::torque();
 
  //std::cout<<"Update_sx"<<update::spin_dir_x[1]<<"\n";
  //std::cout<<"Update_sy"<<update::spin_dir_y[1]<<"\n";
  //std::cout<<"Update_sz"<<update::spin_dir_z[1]<<"\n";
   
   std::cout<<"LAMBDA"<<lambda[0]<<"\n";
   std::cout<<"LAMBDA"<<lambda[1]<<"\n";
   std::cout<<"LAMBDA"<<lambda[2]<<"\n";
   std::cout<<"LAMBDA"<<lambda[3]<<"\n";

   std::cout<<"max_torque_x"<<max_torque_x<<"\n";
   std::cout<<"max_torque_y"<<max_torque_y<<"\n";
   std::cout<<"max_torque_z"<<max_torque_z<<"\n";*/
   }while(iter_spin<=max_iter_spin &&
          max_torque_mod>tol_torque);
         /* max_torque_x>tol_torque &&
          max_torque_y>tol_torque &&
          max_torque_z>tol_torque);//end do..while*/

 utils::compute_global_magnetisation(n_atoms,
                                    sx,sy,sz,
                                    global);

 std::cout<<"VOZ_CONSTRAINT:"<<v0z_constraint<<"\n";
 std::cout<<"GLOBAL"<<global[2]<<"\n";
 double vz_dif = fabs(v0z_constraint - global[2]);
 if(vz_dif>=0.001)
 {
  magnitude_lambda_grad=utils::compute_ROW_vec_magnitude(grad_lambda_param);
  if(magnitude_lambda_grad<0.01)magnitude_lambda_grad=1.0;
  update::lambda_dir_f(lambda,
 		       lambda_dir,
		       grad_lambda_param,
		       magnitude_lambda_grad,
		       alpha_constr);

 
  if(fabs(lambda[0]) > 100.0*input::lambda_x ) {lambda[0]=input::lambda_x; }
  if(fabs(lambda[1]) > 100.0*input::lambda_y ) {lambda[1]=input::lambda_y; }
  if(fabs(lambda[2]) > 100.0*input::lambda_z ) {lambda[2]=input::lambda_z; }
  if(fabs(lambda[3]) > 100.0*input::l1 ) {lambda[3]=input::l1; }

  iter_constr++;
  std::cout<<"Vz_dif: "<<vz_dif<<"\n";
  //test::spin();
  //test::gradient();
  //test::torque();

  total_iter = iter_spin + iter_constr;
 }
 else break;
}while(iter_constr<max_iter_constr);//end do..while
f1.close();


}//end sdescent function

}//end numerical_methods namespace
          
