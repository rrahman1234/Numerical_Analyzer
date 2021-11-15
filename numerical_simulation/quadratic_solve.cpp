#include"quadratic_solve.h"
#include<math.h>
#include<iostream>
#include<vector>
#include <exception>

using namespace std;


quadratic_solve::quadratic_solve(double X, int eq_order, double* coeffs)
{
	
	x = X;
	order = eq_order;
	coefficient = coeffs; 
	cout << "Derived Class: Quadratic Class \n";
}


vector<double> quadratic_solve::solve()
{
	double a = coefficient[0];
	double b = coefficient[1];
	double c = coefficient[2];
	double delta = pow(b,2) - 4.0*a*c; 
        double y1, y2;	
	
	char err_MSG;	
	//catch Exception
	try
	{
		if (delta < 0)
		{
			throw runtime_error("Error with imaginary root \n");
		}
	}

	catch(exception& e)
	{
		cout << e.what() ;
	}

	catch(...)
	{
		cout << "(b^2-4ac) must be greater than zero \n";
	}

	
	y1 = (-b + sqrt(delta))/2.0;
	y2 = (-b - sqrt(delta))/2.0;

	solution.push_back(y1);
	solution.push_back(y2);

	//cout << solution << "\n";
	return solution;
}
