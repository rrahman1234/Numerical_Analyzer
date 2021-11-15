#include"numerical_solver.h"
#include"quadratic_solve.h"
#include"newton_raphson_solve.h"
#include<iostream>
#include<vector>


using namespace std;


int main()
{
	double X = 3.0;
	int order = 2;
        //Quadratic Solve
	//double coff[] = {1.0, 2.0};
	//quadratic_solve quad_Sol(X,order,coff);
	//quad_Sol.solve();
	//quad_Sol.print_result();

	//Newton Raphson
	vector<double> coff;
	vector<double> diff_coff;
	coff.push_back(1.0);
	coff.push_back(3.0);
	coff.push_back(1.0);
	diff_coff.push_back(3.0);
	diff_coff.push_back(2.0);

	newton_raphson_solve NR_Sol(X, 2, 2.0, 50, 0.001, coff, diff_coff);		
	NR_Sol.solve();
	//NR_Sol.function(X, coff);

	return 0;
}
