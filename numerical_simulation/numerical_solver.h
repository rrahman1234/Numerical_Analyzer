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

	public:
		numerical_solver();
		numerical_solver(double x, int order);
		virtual vector<double> solve() = 0;
		void print_result();
		//vector<pair<double, double>> function(double x1, vector<double> poly_coeff, vector<double> diff_poly_coeff);
		double function(double x1, vector<double> poly_coeff);
};


#endif 
