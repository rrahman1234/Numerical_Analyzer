#include"numerical_solver.h"
#include"quadratic_solve.h"
#include"bisection_solve.h"
#include"linear_solve.h"
#include"newton_raphson_solve.h"
#include<iostream>
#include<vector>
#include"armadillo"


using namespace std;
using namespace arma;

int main()
{
	double X = 3.0;
	int order = 1;

    //Quadratic Solve
	//double coff[] = {1.0, 2.0};
	//quadratic_solve quad_Sol(X,order,coff);
	//quad_Sol.solve();
	//quad_Sol.print_result();
    
    int n_rows = 2;
    int n_cols = 2;
	mat Eq(n_rows, n_cols);
    Eq.fill(0.0);
    linear_solve LinSolv(Eq, order, n_rows, n_cols);
    LinSolv.solve();

    ////Newton Raphson
	//vector<double> coff;
	////vector<double> diff_coff;
	//coff.push_back(3.0);
	//coff.push_back(0.0);
	//coff.push_back(-2.0);
	//coff.push_back(1.0);
    //
	////diff_coff.push_back(3.0);
	////diff_coff.push_back(2.0);

	////newton_raphson_solve NR_Sol(X, 2, 2.0, 50, 0.001, coff, diff_coff);		
	////NR_Sol.solve();
	////NR_Sol.function(X, coff);

	//bisection_solve BSolv(-10, 20, 3, 100, 0.0001, coff); 
	//BSolv.solve();
    //BSolv.get_list(2, coff);
    //BSolv.get_list(3, coff);

    

	return 0;
}
