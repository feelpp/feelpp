#include <iostream>
#include <fstream>
#include <cmath>

//using namespace std;

double modulo (double a, double b);

int main (int argc, char *argv[]){
    std::string file_name = "velocity.dat";
	
	double dt = 1.e-3;
	double t = 0.;
	double cycle_time = 0.505;
	
	int nbp = (cycle_time / dt) * 10.;
	
	double *time = new double[nbp];
	double *velocity = new double[nbp];

	//Compute
	t = 0.;
	for (int i = 0; i < nbp; i++){
		time[i] = t;
		velocity[i] =	  28937551.3654412*pow(modulo(t, cycle_time), 10)- 200716534.75997*pow(modulo(t, cycle_time), 9)+ 376821586.480469*pow(modulo(t, cycle_time), 8)- 337967328.956058*pow(modulo(t, cycle_time), 7)+ 170259854.833669*pow(modulo(t, cycle_time), 6)- 50632071.0861583*pow(modulo(t, cycle_time), 5)+ 8764788.99269973*pow(modulo(t, cycle_time), 4)- 820265.547108161*pow(modulo(t, cycle_time), 3)+ 36825.5149442343*pow(modulo(t, cycle_time), 2)- 699.3336471009*modulo(t, cycle_time)+ 44.1999236907;
		t += dt;
	}
	
	//Save
    std::ofstream file(file_name.c_str(), std::ios::out | std::ios::trunc);
	if (file){
		for (int i = 0; i < nbp; i++){
			file << time[i] << "\t" << velocity[i] << std::endl;
		}
		file.close();
	}
	else{
		std::cout << "!!! ERROR !!!" << std::endl;
	}
	
	return 0;
}

double modulo (double a, double b)
{
	double mod = a;
	while (mod >= b)
    {
		mod -= b;
	}
	return mod;
}



