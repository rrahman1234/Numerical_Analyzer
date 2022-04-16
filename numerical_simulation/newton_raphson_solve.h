#include"numerical_solver.h"
#include<iostream>
#include<vector>

#ifndef NEWTON_RAPHSON_SOLVE_H
#define NEWTON_RAPHSON_SOLVE_H


using namespace std;


class newton_raphson_solve: public numerical_solver
{
	private:
		double m_X_init;
		int m_Max_Iter;
		double m_error_tol;
		vector<double>& coefficient;
		vector<double>& diff_coefficients;
	    vector<double> solution_vector;
    public:
		newton_raphson_solve(double X, int eq_order, double X_init, int Max_Iter, double error_tol, vector<double>& coeffs, vector<double>& diff_coefficients);
		void solve();
        vector<double> get_solution();
};


#endif
