//TEST MODULE - SOURCE FILE///
//////////////////////////////



//C++ libraries
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>

//Code defined libraries
#include "input.h"
#include "energy.h"
#include "utils.h"
#include "structure.h"
#include "output.h"
#include "init.h"
#include "gradient.h"
#include "update.h"
#include "solver.h"
namespace test{


//namespace variables



//namespace functions
void initialisation()
{
init::sys();
energy::compute_f();
output::energy(energy::zeeman,
               energy::exchange,
               energy::uniaxial_anisotropy,
               energy::lagrange_term,
               energy::total,
               energy::real);

output::spin(input::no_of_atoms,
                 spin::sx,
                 spin::sy,
                 spin::sz);

output::atom(input::no_of_atoms,
             atom::pos_x,
             atom::pos_y,
             atom::pos_z,
             atom::material_id
             );
output::lists(input::n_materials,
              utils::i_list,
              utils::start,
              utils::end,
              utils::Jmatrix);

}

void spin()
{
 output::spin(input::no_of_atoms,
              spin::sx,
	      spin::sy,
	      spin::sz);

}

void energy()
{
energy::compute_f();
output::energy(energy::zeeman,
               energy::exchange,
               energy::uniaxial_anisotropy,
               energy::lagrange_term,
               energy::total,
               energy::real);


}


void gradient()
{
 gradient::total_f(input::no_of_atoms,
                   gradient::total_x,
       	           gradient::total_y,
            	   gradient::total_z,
                   gradient::real_x,
      	           gradient::real_y,
           	   gradient::real_z);

 output::gradient(input::no_of_atoms,
            	  gradient::zeeman_x,
      		  gradient::zeeman_y,
                  gradient::zeeman_z,
     		  gradient::exchange_x,
    		  gradient::exchange_y,
     		  gradient::exchange_z,
       		  gradient::anisotropy_x,
		  gradient::anisotropy_y,
     		  gradient::anisotropy_z,
		  gradient::lagrange_x,
     		  gradient::lagrange_y,
	  	  gradient::lagrange_z,
	       	  gradient::lambda_param,
		  gradient::total_x,
		  gradient::total_y,
		  gradient::total_z,
		  gradient::real_x,
		  gradient::real_y,
		  gradient::real_z);

}

void torque()
{

/* gradient::torque_f(input::no_of_atoms,
                    spin::sx,
                    spin::sy,
               	    spin::sz,
                    gradient::total_x,
              	    gradient::total_y,
                    gradient::total_z,
                    gradient::torque_x,
     	            gradient::torque_y,
                    gradient::torque_z,
                    gradient::torque_mod);
*/
 output::torque(input::no_of_atoms,
                gradient::torque_x,
                gradient::torque_y,
                gradient::torque_z,
                gradient::torque_mod);

}

void gradient_no_parallel_comp()
{
 gradient::remove_spin_projection(input::no_of_atoms,
  	                          spin::sx,
           	                  spin::sy,
                    	          spin::sz,
                       	          gradient::total_x,
	                          gradient::total_y,
 	    		          gradient::total_z,
	 		          gradient::spin_dot_gradient);

 output::gradient(input::no_of_atoms,
                  gradient::zeeman_x,
                  gradient::zeeman_y,
                  gradient::zeeman_z,
                  gradient::exchange_x,
                  gradient::exchange_y,
                  gradient::exchange_z,
                  gradient::anisotropy_x,
                  gradient::anisotropy_y,
                  gradient::anisotropy_z,
                  gradient::lagrange_x,
                  gradient::lagrange_y,
                  gradient::lagrange_z,
                  gradient::lambda_param,
                  gradient::total_x,
                  gradient::total_y,
                  gradient::total_z,
                  gradient::real_x,
                  gradient::real_y,
                  gradient::real_z);
}

void steepest_descent()
{
 solver::sdescent(input::no_of_atoms,
                    input::v0_constraint_z,
                    gradient::magnitude_spin_term,
                    gradient::magnitude_lambda_term,
              	    gradient::total_x,
                    gradient::total_y,
                    gradient::total_z,
                    gradient::real_x,
                    gradient::real_y,
                    gradient::real_z,
		    gradient::lambda_param,
		    gradient::spin_dot_gradient,
                    gradient::torque_x,
                    gradient::torque_y,
                    gradient::torque_z,
                    gradient::torque_mod,
                    update::spin_dir_x,
              	    update::spin_dir_y,
                    update::spin_dir_z,
              	    update::lambda_dir,
              	    spin::sx,
              	    spin::sy,
              	    spin::sz,
		    spin::global,
              	    spin::lambda);

}

}
