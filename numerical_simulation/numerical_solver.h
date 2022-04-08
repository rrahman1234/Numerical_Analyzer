#include <iostream>
#include <vector>
#include <utility>

#ifndef NUMERICAL_SOLVER_H
#define NUMERICAL_SOLVER_H

using namespace std;

class numerical_solver
{
	protected:
		double x, y;
		int i;
		int order;
		vector<double> solution;
        vector< pair< pair <double, double>, int> > equation_components;

	public:
		numerical_solver();
		numerical_solver(double x, int order);
		//virtual vector<double> solve() = 0;
		void print_result();
        void get_list(double x1, vector<double> poly_coeff);
		double function(double x1, vector<double> poly_coeff);
};


#endif 
