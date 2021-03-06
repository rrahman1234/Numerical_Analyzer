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


void numerical_solver:: get_list(double x1, vector<double>& poly_coeff)
{
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

void numerical_solver::print_result(vector<double>& solution)
{
	cout << "Printing Solutions: \n";
	for(auto s: solution)
	{
		cout << s << "\n";
	}
}

double numerical_solver::function(double x1, vector<double>& poly_coeff)
{
	
	vector<double>* fn = new vector<double>();
 	double func;

	for(int i = 0; i < poly_coeff.size(); i++)
	{
        fn->insert(fn->begin()+i, poly_coeff[i]*pow(x1, poly_coeff.size()-1-i));	
    }

	func = accumulate(fn->begin(), fn->end(), 0.0);
	
	return func;
}

vector<double>* numerical_solver::function_poly_terms(double x1, vector<double>& poly_coeff)
{
	int i=0;
	
	vector<double>* fn;
 	double func;

	for(int i = 0; i < poly_coeff.size(); i++)
	{
        fn->insert(fn->begin()+i, poly_coeff[i]*pow(x1, poly_coeff.size()-1-i));	
        i += 1;
	}
	
	return fn;
}

vector<double>* numerical_solver::function_poly_diff_terms(double x1, vector<double>& poly_coeff)
{
	int i=0;
	
	vector<double>* fn;
 	double func;

	for(size_t i=0; i<poly_coeff.size(); i++)
	{
        if(i > 0)
        {
            fn->insert(fn->begin()+i, poly_coeff[i]*pow(x1, poly_coeff.size()-1-i-1));	
        }
    }
	
	return fn;
}

