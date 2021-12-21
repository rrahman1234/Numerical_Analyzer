#include "numerical_solver.h"
#include <iostream>
#include<vector>
#include <utility>
#include<math.h>
#include<numeric>


using namespace std;

numerical_solver::numerical_solver()
{
	x = 0.0;
	order = 1;

}


numerical_solver::numerical_solver(double x, int order): x(x), order(order)
{
	cout << "BaseClass \n";
}


void numerical_solver:: get_list(double x1, vector<double> poly_coeff)
{
    cout << "Printing Equation \n" << endl;

    int i = 0;
    for(auto it: poly_coeff)
	{
		equation_components.push_back(make_pair(make_pair(poly_coeff[i], poly_coeff[i]*pow(x1, i)), i));	
		i += 1;

	}
    
    cout << " ****** Printing Terms ******" << endl; 
    for (const auto& it : equation_components)
    {
        cout << "Coefficient:" << it.first.first << " Value of different terms:" << it.first.second << " Power of term: " << it.second << endl;
    }
    cout << " ****** Printing Terms: Done ******" << endl; 
    
} 

void numerical_solver::print_result()
{
	cout << "Printing Solutions: \n";
	for(auto s: solution)
	{
		cout << s << "\n";
	}
}

double numerical_solver::function(double x1, vector<double> poly_coeff)
{
	int i=0;
	
	vector<double> fn;
 	double func;

	for(auto it: poly_coeff)
	{
		fn.push_back(poly_coeff[i]*pow(x1, i));	
		i += 1;

	}

	func = accumulate(fn.begin(), fn.end(), 0.0);
	

	return func;
}





/*
vector<pair<double, double>>  numerical_solver::function(double x1, vector<double> poly_coeff, vector<double> diff_poly_coeff)
{
	
	int i=0;   	

	//x1 = x;
	vector<pair<double, double>> polynomial;
	vector<double> fn;
	vector<double> dfn;

	
	for(int j = 0; j < poly_coeff.size(); j++)
	{
		polynomial.push_back(make_pair(poly_coeff[i]*pow(x1, i), diff_poly_coeff[i]*pow(x1, i)));
		i += 1;
	}

	for(int i=0; i<polynomial.size(); ++i)
	{
  		fn.push_back(polynomial[i].first);
  		dfn.push_back(polynomial[i].second);
	}

	
	for(int i=0; i<polynomial.size(); ++i)
	{
  		std::cout << polynomial[i].first << '\n';
	}

	return polynomial;

}*/


//numerical_solver::numerical_solver(double variable, int eq_order)
//{
//	x = variable;
//	order = eq_order;
//	cout << "Base Class \n";
//}
