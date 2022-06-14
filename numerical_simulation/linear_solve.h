#include"numerical_solver.h"
#include<iostream>
#include<vector>
#include"armadillo"
#include <Eigen/Dense>
#include <string>  

#ifndef LINEAR_SOLVE_H
#define LINEAR_SOLVE_H

using namespace std;

class linear_solve: public numerical_solver
{
    private:
        int num_rows, num_cols;
        int order;
        int num_iterations;
        string solver_type; 
        arma::mat EqnMat;
        arma::vec b_right_side;
        Eigen::MatrixXd EigenEqMat;
        Eigen::VectorXd Eigen_b_right_side;
        Eigen::VectorXd solution;

    public:
        typedef std::vector<double> stdvec;
        linear_solve(arma::mat& LinEqs, arma::vec& b_eq, int num_rows, int num_cols);
        linear_solve(Eigen::MatrixXd& LinEqs, Eigen::VectorXd& b_eq, int num_rows, int num_cols, string solver_name, int num_iter = 100);
        void solve(); 
        vector<double> get_solution();
        void split_matrix();
};

#endif
