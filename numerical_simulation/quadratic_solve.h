#include"numerical_solver.h"
#include<iostream>
#include<vector>

#ifndef QUADRATIC_SOLVE_H
#define QUADRATIC_SOLVE_H

using namespace std;


class quadratic_solve: public numerical_solver
{
	private:
		double* coefficient;

	public:
		quadratic_solve(double X, int eq_order, double* coeffs);
		vector<double> solve();
};



#endif

