 //UTILS MODULE - SOURCE FILE//
//////////////////////////////


//C++ libraries
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>//std::distance

//Code defined libraries
#include "utils.h"
#include "input.h"
#include "structure.h"
#include "gradient.h"



//DEFINING THE MODULES' FUNCTIONS & INITIALISING THE EXTERNAL VARIABLES


namespace utils {


//namespace variables

std::vector<std::vector<double> >Jmatrix;
std::vector<double>global_magnetisation = {0,0,0};
std::vector<double> sums = {0,0,0,0};
std::vector<int>i_list;
std::vector<int>start;
std::vector<int>end;




std::vector<double>temp_x;
std::vector<double>temp_y;
std::vector<double>temp_z;
std::vector<double>temp_lambda;
std::vector<double>temp_gradient;
double temp_energy;





//namespace functions
//
//
void generate_exchange_mat(int n_materials,
			   std::vector<double>&exch_list,
			   std::vector<std::vector<double> >&Jmat)
{
 Jmat.resize(n_materials, std::vector<double>(n_materials));

 for(int i=0;i<n_materials;i++)
 {
  for(int j=0;j<n_materials;j++)
  {
   if(abs(j-i)>1.5)//needs to be only greater than 1, but for safety..
   { 
    Jmat[i][j]=0;
    continue;
   }
   else
   {
    Jmat[i][j] = exch_list[i+j];
   }    
  }
 }


}

void generate_interac_list(int n_atoms,
			   std::vector<double> &x,
			   std::vector<double> &y,
			   std::vector<double> &z,
			   std::vector<int> &int_list,
			   std::vector<int> &start,
			   std::vector<int> &end)
{
start.resize(n_atoms);
end.resize(n_atoms);
double r0=1.0;//range of interactions(currently equal with the unit cell step)
const double rtoll=r0*1.0e-5;

for(int i=0; i<n_atoms; i++)
{
 int count_interactions=0;//count the total interactions an atom has
 for(int j=0; j<n_atoms; j++)
  {
    if(i==j)continue;//exclude self interactions
    else
    {
    //compute the distance between atoms
    double rij=sqrt( pow((x[i]-x[j]),2.0) +
		     pow((y[i]-y[j]),2.0) +
		     pow((z[i]-z[j]),2.0) );
    if(rij-r0<rtoll)//if the distance agrees with the interaction range..
    {
     count_interactions++;
     int_list.push_back(j);//update interaction list
     int index_j=std::distance(int_list.begin(),int_list.end());
     end[i]=index_j;//update end list
    }
    } 
  }
  start[i]=end[i]-count_interactions;
}

}


void compute_global_magnetisation(int n_atoms,
				  std::vector<double> &sx,
				  std::vector<double> &sy,
				  std::vector<double> &sz,
				  std::vector<double> &global)
{
 double modulus = 0, term1 = 0, term2 = 0, term3 = 0;
 
 for(int count_atom=0;count_atom<n_atoms;count_atom++)
 {
  term1 = term1 + sx[count_atom];
  term2 = term2 + sy[count_atom];
  term3 = term3 + sz[count_atom];

 }

  modulus = sqrt( pow(term1,2.0) + pow(term2,2.0) + pow(term3,2.0)  );
  global[0]= term1/modulus;
  global[1] = term2/modulus;
  global[2] = term3/modulus;
 


}


double compute_XYZROW_vec_magnitude(std::vector<double> &x,
	           	            std::vector<double> &y,
			            std::vector<double> &z)
{

 double sum=0.0;	
 for(int i=0;i<x.size();i++)
 {
  sum += x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
 }
 return sqrt(sum);
 
}

double compute_XYZ_vec_magnitude(double &x,
				 double &y,
				 double &z)
{

return sqrt(x*x + y*y + z*z);
}

double compute_ROW_vec_magnitude(std::vector<double> &vec)
{
double sum=0.0;
for(int i=0; i<vec.size();i++)
{
 sum +=vec[i]*vec[i];
}
 return sqrt(sum);

}




double max_val_array(std::vector<double> &vec, double& max_val)
{

max_val = fabs(vec[0]);

for (int i = 0; i<vec.size(); i++)
{
    if(max_val<fabs(vec[i]))max_val=fabs(vec[i]);
}

return max_val;

}




void print_all(std::ofstream& f1, std::ofstream& f2,std::ofstream& f3,int step)
{




for(int i = 0;i<input::no_of_atoms;i++)
{


   f1<<std::setprecision(10);
   f1<<i<< " " <<step<<" "<<spin::sx[i]<<" "<<spin::sy[i]<<" "<<spin::sz[i]<<" "<<sqrt(spin::sx[i]*spin::sx[i]+spin::sy[i]*spin::sy[i]+spin::sz[i]*spin::sz[i]) //1-6
        << " " << gradient::torque_x[i]<< " " << gradient::torque_y[i]<< " " << gradient::torque_z[i] <<" "<<gradient::torque_mod[i] //7-10
	<< " " <<gradient::total_x[i]<<" "<<gradient::total_y[i]<<" "<<gradient::total_z[i]<<" "<<gradient::lambda_param[0]<<" "<<gradient::lambda_param[1]<<" "
											        <<gradient::lambda_param[2]<<" "<<gradient::lambda_param[3]<<" "
												<<gradient::magnitude_spin_term<<" "<<gradient::magnitude_lambda_term// 11-19
	<< " " <<spin::lambda[0]<<" "<<spin::lambda[1]<<" "<<spin::lambda[2]<<" "<<spin::lambda[3]<<"\n";//20-23


   /*
   f3 <<std::setprecision(10)  <<
     i <<" "<<spin::x[i]<<" "<<spin::y[i]<<" "<<spin::z[i]<<" " <<v0_constrain_z<<" "<<  // 1-5
	" "<<gradient::total[i]<<" "<<gradient::total[i+no_of_atoms]<<" "<<gradient::total[i+2*no_of_atoms]<< //6-8
	" "<< gradient::exchange_x[i]   <<" "<< gradient::exchange_y[i] <<" "<< gradient::exchange_z[i] << //9-11
	" "<< gradient::anisotropy_x[i] <<" "<< gradient::anisotropy_y[i] <<" "<< gradient::anisotropy_z[i] << //12-14
	" "<< gradient::zeeman_x[i]     <<" "<< gradient::zeeman_y[i] <<" "<< gradient::zeeman_z[i] << //15-17
	" "<< gradient::lagrangian_x[i] <<" "<< gradient::lagrangian_y[i] <<" "<< gradient::lagrangian_z[i] << //18-20
    " " << spin::lambda[0] << " "<< spin::lambda[1] << " "<< spin::lambda[2] << " "<< spin::lambda[3] << //21-24
    " " << gradient::torque_X[i] << " "<< gradient::torque_Y[i] << " "<< gradient::torque_Z[i] <<  //25-27
    " "<< step<<"\n";
*/




// f1<<spin::lambda[i]<<" "<<spin::lambda[i+1]<<" "<<spin::lambda[i+2]<<" "<<" "<<spin::lambda[i+3] << "  "<<gradient::total[i]<<" "<<gradient::total[i+no_of_atoms]<<" "<<gradient::total[i+2*no_of_atoms]<<" "<<gradient::total[i+3*no_of_atoms + 3]<<"\n";

}
//f1<<"\n\n";


//f2<<std::setprecision(10);
//f2<<energy::real_energy<<" "<<utils::compute_modul(gradient::real_energy)<<" "<<step<<" "<<v0_constrain_z<<"\n";
//f2<<(energy::real_energy + J*(no_of_atoms-1))/(5)<<" "<<utils::compute_modul(gradient::real_energy)<<" "<<step<<" "<<v0_constrain_z<<"\n";
//f2<<(energy::real_energy + J*(no_of_atoms-1))/(no_of_atoms+0.0)<<" "<<utils::compute_modul(gradient::real_energy)<<" "<<step<<" "<<v0_constrain_z<<"\n";
//f2<<(energy::real_energy/(no_of_atoms))<<" "<<utils::compute_modul(gradient::real_energy)<<" "<<step<<" "<<v0_constrain_z<<"\n";



 //   <<gradient::total[*no_of_atoms]<<" "<<gradient::total[4*no_of_atoms + 1]<<gradient::total[4*no_of_atoms + 2];



}

/*

void save_state(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                std::vector<double>& lambda,
                std::vector<double>& gradient,
                double& energy)
{

x.resize(no_of_atoms);
y.resize(no_of_atoms);
z.resize(no_of_atoms);
lambda.resize(no_of_atoms+3);
gradient.resize(4*no_of_atoms+3);


x = spin::x;
y = spin::y;
z = spin::z;
lambda = spin::lambda;
gradient = gradient::total;
energy = energy::total;




}


double compute_modul(std::vector<double> vec)
{
    double sum = 0;

    for(unsigned int i = 0;i<vec.size();i++)
    {
        sum = sum + vec[i]*vec[i];

    }

    return sqrt(sum);
}

void compute_magnetisation_sums(std::vector<double> x, std::vector<double> y, std::vector<double> z,
                                std::vector<double>& sums)
{


    double modulus = 0, term1 = 0, term2 = 0, term3 = 0;

    for(int i = 0;i < no_of_atoms;i++)
    {

    term1 = term1 + x[i];
    term2 = term2 + y[i];
    term3 = term3 + z[i];

    }

    modulus =  pow(term1,2.0) + pow(term2,2.0) + pow(term3,2.0); // sum of the x components squared + sum of the y components squared + sum of the z components squared

    sums[0] = term1;
    sums[1] = term2;
    sums[2] = term3;
    sums[3] = modulus;



}


*/

}
