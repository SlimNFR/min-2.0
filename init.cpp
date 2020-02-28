//INITIALISE MODULE - SOURCE FILE//
//////////////////////////////////


//C++ libraries
#include <vector>
#include <iostream>

//Code defined libraries
#include "init.h"
#include "structure.h"
#include "input.h"
#include "utils.h"
#include "gradient.h"
#include "update.h"

namespace init{


void alloc_memory()
{

spin::sx.resize(input::no_of_atoms,0.0);
spin::sy.resize(input::no_of_atoms,0.0);
spin::sz.resize(input::no_of_atoms,0.0);
spin::lambda.resize(4,0.0);

atom::pos_x.resize(input::no_of_atoms,0.0);
atom::pos_y.resize(input::no_of_atoms,0.0);
atom::pos_z.resize(input::no_of_atoms,0.0);
atom::material_id.resize(input::no_of_atoms,0.0);

gradient::zeeman_x.resize(input::no_of_atoms,0.0);
gradient::zeeman_y.resize(input::no_of_atoms,0.0);
gradient::zeeman_z.resize(input::no_of_atoms,0.0);
gradient::exchange_x.resize(input::no_of_atoms,0.0);
gradient::exchange_y.resize(input::no_of_atoms,0.0);
gradient::exchange_z.resize(input::no_of_atoms,0.0);
gradient::anisotropy_x.resize(input::no_of_atoms,0.0);
gradient::anisotropy_y.resize(input::no_of_atoms,0.0);
gradient::anisotropy_z.resize(input::no_of_atoms,0.0);
gradient::lagrange_x.resize(input::no_of_atoms,0.0);
gradient::lagrange_y.resize(input::no_of_atoms,0.0);
gradient::lagrange_z.resize(input::no_of_atoms,0.0);
gradient::lambda_param.resize(4,0.0);
gradient::total_x.resize(input::no_of_atoms,0.0);
gradient::total_y.resize(input::no_of_atoms,0.0);
gradient::total_z.resize(input::no_of_atoms,0.0);
gradient::real_x.resize(input::no_of_atoms,0.0);
gradient::real_y.resize(input::no_of_atoms,0.0);
gradient::real_z.resize(input::no_of_atoms,0.0);
gradient::torque_x.resize(input::no_of_atoms,0.0);
gradient::torque_y.resize(input::no_of_atoms,0.0);
gradient::torque_z.resize(input::no_of_atoms,0.0);
gradient::torque_mod.resize(input::no_of_atoms,0.0);
gradient::spin_dot_gradient.resize(input::no_of_atoms,0.0);

update::lambda_dir.resize(4,0.0);
update::spin_dir_x.resize(input::no_of_atoms,0.0);
update::spin_dir_y.resize(input::no_of_atoms,0.0);
update::spin_dir_z.resize(input::no_of_atoms,0.0);

}


void sys()
{
//generate lattice
atom::generate_f(input::no_of_atoms,
		input::nx,
		input::ny,
		input::nz,
		input::natoms_material,
		atom::pos_x,
		atom::pos_y,
		atom::pos_z,
		atom::material_id);
//initial spin values
spin::init_config_f(input::no_of_atoms,
   		    input::s0x,
                    input::s0y,
                    input::s0z,
		    spin::sx,
		    spin::sy,
		    spin::sz,
		    atom::material_id);
//initial constraints
spin::init_constr_f(input::lambda_x,
                    input::lambda_y,
                    input::lambda_z,
                    input::l1,
                    spin::lambda); 
//generate the exchange matrix [n_mat]X[n_mat]
utils::generate_exchange_mat(input::n_materials,
                             input::exchange_list,
                             utils::Jmatrix);

utils::generate_interac_list(input::no_of_atoms,
                             atom::pos_x,
                             atom::pos_y,
                             atom::pos_z,
                             utils::i_list,
                             utils::start,
                             utils::end);

}



}
