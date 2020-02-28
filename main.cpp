//MAIN MODULE - SOURCE FILE///
//////////////////////////////


//----------------------------------------------------------------APPEALING NAME FOR ME CODE: BLABLA V.1.0-------------------------------------------------------//


//C++ libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <fenv.h>
#include <chrono>

//Code defined libraries
#include "input.h"
#include "output.h"
#include "init.h"
#include "energy.h"
#include "utils.h"
#include "structure.h"
#include "test.h"
#include "gradient.h"
#include "solver.h"

int main()
{
//feenableexcept(FE_INVALID | FE_OVERFLOW);
    auto start = std::chrono::high_resolution_clock::now();

    std::ofstream f1;
    std::ofstream f2;
    std::ofstream f3;


    f1.open("spins.txt", std::ofstream::out);
    f2.open("energy_and_gradient.txt",std::ofstream::out);
    f3.open("all.txt",std::ofstream::out);
    init::alloc_memory();
    init::sys();
    test::gradient();
  //  test::energy();
 //   test::spin();
//    test::gradient_no_parallel_comp();
    test::steepest_descent();
    //test::spin();
    //test::torque();
    //std::cout<<"SPINI: "<<spin::x[0]<<","<<spin::y[0]<<","<<spin::z[0]<<"\n";

//      utils::save_state(utils::temp_x,utils::temp_y,utils::temp_z,utils::temp_lambda,utils::temp_gradient,utils::temp_energy);



      //---------------------------------------compute total energy and total gradient --------------------//

//      energy::total_f();
//      gradient::total_f();
//      std::cout<<utils::compute_modul(gradient::total);


/*

     // ---------------------------------------------------test steepest descent--------------------------------------//
            double temp_en;
            double mod_grad;
            double max_torque_before;            
            double ax=0, ay=0, az=1;   
            double vz_dif;         

            double lx ,ly, lz, ll1;            

            double max_iter = 5e3;
            double tol = 1e-5;
            double alpha = 0.01;
            double max_torque; 

            int j;//pas total
            int k; //step for each energy minimum calculation; resets with each new global magnetisation constraint
            int itot=10; //this value shall be used later to implement loop counters in the energy barrier calculations 

		
           //make sure spins point initially on z
            std::fill (spin::x.begin(),spin::x.end(),0.001);
            std::fill (spin::y.begin(),spin::y.end(),0.001);
            std::fill (spin::z.begin(),spin::z.end(),0.999);
	    
	    //main for loop: the l parameter sets the global value of the magnetisation
            for(int l = 0; l>=0 ; l=l-1)
            //  for(int l = 0; l>=-1 ; l=l-1)
            {

            v0_constrain_z = (double)l/(itot+0.0);//set the global magnetisation
            //v0_constrain_z = 0.0;
	    
        //            if (l == 0) {
//            std::fill (spin::x.begin(),spin::x.end(),0);
//            std::fill (spin::y.begin(),spin::y.end(),1);
//            std::fill (spin::z.begin:q(),spin::z.end(),0);
//            }

	    //make sure the first spin points in the 0,0,1 direction
          // spin::x[0]=ax;
          // spin::y[0]=ay;
          // spin::z[0]=az;
	   //std::fill (spin::x.begin(),spin::x.end(),0.001);
          // std::fill (spin::y.begin(),spin::y.end(),0.001);
           //std::fill (spin::z.begin(),spin::z.end(),0.999);
	    
	   //if(v0_constrain_z<0.0){spin::lambda[3]=-l1;}
    
	    //make sure the last spin points in the 0,0,-1 direction//this can help acquiring the DW configuration
           spin::x[no_of_atoms-1]= ax;
           spin::y[no_of_atoms-1]= ay;
           spin::z[no_of_atoms-1]=-az;

	    //save the values of the lagrangian parameters before you do any updates using the search vector
            lx = spin::lambda[0];
	    ly = spin::lambda[1];
	    lz = spin::lambda[2];
	    ll1 = spin::lambda[3];
	
	    j = 1;
	    do
	    {

	    k = 1;

	    //do..while statement which searches for the minimum energy configuration at a given global magnetisation state
            do
            {
		//compute the global magn, the energy and the gradient    
                utils::compute_global_magnetisation();
                energy::total_f();
                temp_en =energy::total;
                gradient::total_f();

			    //1.compute the modulus of the gradient vector and reset its value if it's smaller than a certain lower margin
			    //2. update the search vector with the newly found directions
                            for(unsigned int i = 0; i<3*no_of_atoms;i++)    
			    {
                               mod_grad = utils::compute_modul(gradient::total);
                               if (mod_grad<0.01) mod_grad=0.1;
                               s_vec_spins[i] = gradient::total[i]/mod_grad;
                           //    std::cout<< "i= " <<i <<"  " <<gradient::total[i]<<" "<<mod_grad<<std::endl<<std::endl;
                            }

		   //update the spin variables using the newly found search vector	    
                   update_spins(spin::x,spin::y,spin::z,s_vec_spins,alpha);


	           //check at a certain number of steps whether the max_torque_before is getting smaller than the imposed tolerance. if not, reset the lagrange parameters with
		   //the ones prior to having any problems, see below.
		   //
		   //shouldn't this be k%1000 == 0 ????????????????????????????????
                    if ( (k%500 == 0) && max_torque_before > tol ){ 
                            spin::lambda[0]=lx;
                            spin::lambda[1]=ly;
                            spin::lambda[2]=lz;
                            spin::lambda[3]=ll1;
                	    }
			  

		   //compute the energy, global magn, gradient and the absolute difference of the z components for the actual and theoretical global magnetisation
                   energy::total_f();
                   utils::compute_global_magnetisation();
                   gradient::total_f();
                   vz_dif= fabs(spin::global[2]-v0_constrain_z);


	       //compute the new maximum value of the gradient array   
               utils::max_val_array(gradient::torque,max_torque);
	       //update k step
               k++;
              }while(k<= max_iter && max_torque > tol);


	    	if(vz_dif>=0.01 || max_torque>tol)
		{

		   for(unsigned int i = -3; i<4;i++)
                            {
                              s_vec_lambda[i] = gradient::total[i+3*no_of_atoms]/mod_grad;
			       mod_grad = utils::compute_modul(gradient::total);
                               if (mod_grad<0.01) mod_grad=1;
                               
                           //    std::cout<< "i= " <<i <<"  " <<gradient::total[i]<<" "<<mod_grad<<std::endl<<std::endl;
                            }

		   //if the updated lagrange parameters are bigger by a certain margin compared to the initial ones (in the globals.h file), then reset them
                   if(fabs(spin::lambda[0]) > 10.0*lambda_x ) {  spin::lambda[0]=lambda_x; }
                   if(fabs(spin::lambda[1]) > 10.0*lambda_y ) {  spin::lambda[1]=lambda_y; }
                   if(fabs(spin::lambda[2]) > 10.0*lambda_z ) {  spin::lambda[2]=lambda_z; }
                   if(fabs(spin::lambda[3]) > 10.0*l1 )       {  spin::lambda[3]=l1; }



                   //update the spin variables using the newly found search vector
                   update_lambda(spin::lambda,s_vec_lambda,alpha);
		j = j+1;

		}

		else
		{
		   //if everything was smooth, save the current values of the lagrange parameters as checkpoints
                   lx = spin::lambda[0];
                   ly = spin::lambda[1];
                   lz = spin::lambda[2];
                   ll1 = spin::lambda[3];

		  break;
		} 

                }while(j<250);
			


     //       std::cout<<"!!!!!!!!!!!!!!!!!!    max_torque "<< max_torque << "  vz_dif=" << vz_dif <<"   "<< k << std::endl;
//            system("echo in code");
  //          std::string cp = "cp " + file2 + "energy_and_gradient0.dat; cp "+ file3 + "all0.dat; gnuplot ../plot_num_analytic.gnu;";
	      std::cout<<"Zeeman energy is: "<<energy::zeeman<<"\n";         

	    //print at the end of each minimum energy calculation and update the global step counter, j.
            utils::print_all(f1,f2,f3,k);
        

           //screen output
            std::cout<<"!!!!!!!!!!!!!!!!!!    mz= "<< v0_constrain_z << "  <sz>=" << spin::global[2] <<"  no. of lambda updates= "<< j <<
                       "  max_torque= "<< max_torque << "  energy= "<<energy::real_energy<<"\n";

//	    //more screen outpu
//            for(int i = 0; i<no_of_atoms;i++)
///              {
//	         // std::cout<<"gradient["<<i<<"]="<<gradient::total[i+0*no_of_atoms]<< "  "<<gradient::total[i+1*no_of_atoms]<< "  "<<gradient::total[i+2*no_of_atoms] <<"  mod_grad=  " << mod_grad <<std::endl;
//              }
///

            }
	 auto finish = std::chrono::high_resolution_clock::now();
     	 std::chrono::duration<double> elapsed = finish - start;
	 std::cout<<" Elapsed time:"<<elapsed.count() <<" s"<<"\n";

 
	



            f1.close();
            f2.close();
            f3.close();
*/	
    return 0;
}
