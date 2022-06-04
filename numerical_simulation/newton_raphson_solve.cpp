#include "numerical_solver.h"
#include "newton_raphson_solve.h"
#include<math.h>
#include<iostream>
#include<vector>
#include<utility>

using namespace std;


newton_raphson_solve::newton_raphson_solve(double X, int eq_order, double X_init, int Max_Iter, double error_tol, vector<double>& coeffs, vector<double>& diff_coeffs): m_X_init(X_init), m_Max_Iter(Max_Iter), m_error_tol(error_tol), coefficient(coeffs), diff_coefficients(diff_coeffs)
{
	x = X;
	order  = eq_order;
	cout << "Derived Class: Newton Raphson Class \n";
}


void newton_raphson_solve::solve()
{
	cout << "Solve \n";
   	
	double xn, xn_plus_1;
	double f1, f0;
	double df0;
	int iter = 0;
 	
	xn = m_X_init;
	f0 = function(xn, coefficient);
	df0 = function(xn, diff_coefficients); 	
	f1 = function(xn, coefficient);
	do
	{
		xn_plus_1 = xn - (f0/df0);
		cout<<"Iteration: "<< iter <<" x = " << xn_plus_1 << " and f(x) = "<< f1 << "\n";
	        xn = xn_plus_1;	
		iter++;
		f0 = function(xn,coefficient); 	
		f1 = function(xn_plus_1, coefficient);
		df0 = function(xn,diff_coefficients); 
	} while(fabs(f1) > m_error_tol);
	
    solution_vector.push_back(xn_plus_1); 	

}

vector<double> newton_raphson_solve::get_solution()
{
    cout << "Solution: Newton-Raphson" << endl;
    return solution_vector;
}

	
