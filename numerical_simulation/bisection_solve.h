#include"numerical_solver.h"
#include<iostream>
#include<vector>

#ifndef BISECTION_SOLVE_H
#define BISECTION_SOLVE_H


using namespace std;


class bisection_solve: public numerical_solver
{
    private:
		double x_left;
		double x_right;
		double m_error_tol;
		int m_Max_Iter;
		int order;
		vector<double>& coefficient;
	    
    public:
		bisection_solve(double X_left, double X_right, int eq_order, int Max_Iter, double error_tol, vector<double>& coeffs);
		void solve();
        vector<double> get_solution();
};

#endif
