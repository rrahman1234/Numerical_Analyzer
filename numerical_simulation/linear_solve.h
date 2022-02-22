#include"numerical_solver.h"
#include<iostream>
#include<vector>
#include"armadillo"


#ifndef LINEAR_SOLVE_H
#define LINEAR_SOLVE_H

using namespace std;
using namespace arma;


class linear_solve: public numerical_solver
{
    private:
        int num_rows, num_cols;
        int order;
        mat EqnMat;
        vec b_right_side;

    public:
        typedef std::vector<double> stdvec;
        linear_solve(mat LinEqs, vec b_eq, int eq_order, int num_rows, int num_cols);
		vector<double> solve();
};

#endif
