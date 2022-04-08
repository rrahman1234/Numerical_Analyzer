#include"numerical_solver.h"
#include<iostream>
#include<vector>
#include"armadillo"
#include <Eigen/Dense>
#include <string>  

#ifndef LINEAR_SOLVE_H
#define LINEAR_SOLVE_H

using namespace std;
using namespace arma;
using namespace Eigen;

class linear_solve: public numerical_solver
{
    private:
        int num_rows, num_cols;
        int order;
        string solver_type; 
        mat EqnMat;
        vec b_right_side;
        MatrixXd EigenEqMat;
        VectorXd Eigen_b_right_side;

    public:
        typedef std::vector<double> stdvec;
        linear_solve(mat LinEqs, vec b_eq, int eq_order, int num_rows, int num_cols);
        linear_solve(MatrixXd LinEqs, VectorXd b_eq, int eq_order, int num_rows, int num_cols);
		//vector<double> solve();
        //VectorXd solve_Eigen();
        //VectorXd solve();
        void solve(vector<double>& solution); 
        void solve(VectorXd& solution); 

};

#endif
