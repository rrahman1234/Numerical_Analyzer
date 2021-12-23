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



vector<double> bisection_solve::solve()
{
	cout << "Solve \n";

	double x_mid;
	double f1, f2;
    int Iter = 0;


	//try
    //{
    //    if (f1 * f2 >= 0) 
	//    {
    //        throw "You have not assumed right a and b\n";
    //    }
    // }
    //
    //catch (const char* msg) 
    //{
    //    cerr << msg << endl;
    //}
    do 
    {
	    f1 = function(x_left, coefficient);
	    f2 = function(x_right, coefficient);
        if (f1 * f2 >= 0) 
	    {
            cout << "Choose new interval" << endl;
            break;
        }
        else
        {
            // Find middle point
            x_mid = (x_right + x_left)/2;
            // Check if middle point is root
            //if (function(x_mid, coefficient) == 0.0)
                //break;
            // Decide the side to repeat the steps
            if (function(x_mid, coefficient)*function(x_left, coefficient) < 0)
              x_left = x_mid;
            else
              x_right = x_mid;
            cout << "Iteration:" << Iter << "\n";
            Iter++;
        }
    } while ((x_right-x_left) >= m_error_tol && function(x_mid, coefficient) != 0.0);
    cout << "The value of root is : " << x_mid << "\n";
    solution.push_back(x_mid);
    return solution;
}
	

