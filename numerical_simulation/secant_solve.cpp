#include "numerical_solver.h"
#include "newton_raphson_solve.h"
#include "secant_solve.h"
#include<math.h>
#include<iostream>
#include<vector>
#include<utility>
#include<cmath>

using namespace std;

secant_solve::secant_solve(double X_left, double X_right, int eq_order, int Max_Iter, double error_tol, vector<double>& coeffs):x_left(X_left), x_right(X_right), m_Max_Iter(Max_Iter), m_error_tol(error_tol), order(eq_order), coefficient(coeffs)
{
	cout << "Derived Class: Secant Method Class \n";
}

void secant_solve::solve()
{
	cout << "Solve \n";

	double x_mid;
	double f1, f2, dx, xl, rts;
    int Iter = 0;

	f1 = function(x_left, coefficient);
	f2 = function(x_right, coefficient);
    
    if (f1 * f2 >= 0) 
	{
        cout << "Choose new interval" << endl;
        exit(0);
    }
    
    //if (fabs(f1) < fabs(f2)) {
	//	rts=x_left;
	//	xl=x_right;
	//	swap(f1,f2);
	//} else {
	//	xl=x_left;
	//	rts=x_right;
	//}
    
    if (fabs(f1) < fabs(f2)) {
		rts=x_right;
		xl=x_left;
		swap(f1,f2);
	} else {
		xl=x_right;
		rts=x_left;
	}
	
    for (int j=0; j<m_Max_Iter; j++) {
		dx=(xl-rts)*f2/(f2-f1);
		xl=rts;
		f1=f2;
		rts += dx;
		f2=function(rts, coefficient);
		if (fabs(dx) < m_error_tol || f2 == 0.0)
        {
            x_mid = rts;
            break;
        }
        cout << "Iteration:" << j << " Left: " << xl << " " << " Right: " << xl+dx << " Mid: " << x_mid << " Function: " << f2 << "  diff: " << dx << "---" << endl;
	}

    cout << "The value of root is : " << x_mid << "\n";
    solution_vector.push_back(x_mid);
}
	
vector<double> secant_solve::get_solution()
{
    cout << "Solution: Bisection Method" << endl;
    return solution_vector;
}

