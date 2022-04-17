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
        vector< pair< pair <double, double>, int> > equation_components;
        vector<double> solution_vector;
        vector<double> solution;

	public:
		numerical_solver();
		numerical_solver(double x, int order);
		void print_result(vector<double>& solution);
        void get_list(double x1, vector<double>& poly_coeff);
		double function(double x1, vector<double>& poly_coeff);
        virtual vector<double> get_solution() = 0;
		virtual void solve() = 0;
};


#endif 
