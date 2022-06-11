#include "numerical_solver.h"
#include "newton_raphson_solve.h"
#include "bisection_solve.h"
#include<math.h>
#include<iostream>
#include<vector>
#include<utility>

using namespace std;

bisection_solve::bisection_solve(double X_left, double X_right, int eq_order, int Max_Iter, double error_tol, vector<double>& coeffs):x_left(X_left), x_right(X_right), m_Max_Iter(Max_Iter), m_error_tol(error_tol), order(eq_order), coefficient(coeffs)
{
	cout << "Derived Class: Bisection Method Class \n";
}

void bisection_solve::solve()
{
	cout << "Solve \n";

	double x_mid;
	double f1, f2;
    int Iter = 0;

	f1 = function(x_left, coefficient);
	f2 = function(x_right, coefficient);
    if (f1 * f2 >= 0) 
	{
        cout << "Choose new interval" << endl;
        exit(0);
    }
    x_mid = x_left;
    do 
    {
        cout << "Iteration:" << Iter << " Left: " << x_left << " " << " Right: " << x_right << " Mid: " << x_mid << " Function: " << function(x_mid, coefficient) << "  diff: " << (x_right - x_left) << "---" << endl;
        x_mid = (x_right + x_left)/2;
        if (function(x_mid, coefficient)*function(x_left, coefficient) < 0)
          x_right = x_mid;
        else
          x_left = x_mid;
            Iter++;
    } while ((x_right-x_left) >= m_error_tol);
    cout << "The value of root is : " << x_mid << "\n";
    solution_vector.push_back(x_mid);
}
	
vector<double> bisection_solve::get_solution()
{
    cout << "Solution: Bisection Method" << endl;
    return solution_vector;
}

