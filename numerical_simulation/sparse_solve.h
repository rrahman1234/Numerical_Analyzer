#include"numerical_solver.h"
#include<iostream>
#include<vector>
#include"armadillo"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <string>  

#ifndef SPARSE_SOLVE_H
#define SPARSE_SOLVE_H

using namespace std;
using namespace arma;
using namespace Eigen;


typedef Eigen::SparseMatrix<double> SpMatrx;


class sparse_solve: public numerical_solver
{
    private:
        int num_rows, num_cols;
        int order;
        string solver_type; 
        sp_mat EqnMat;
        vec b_right_side;
        SpMatrx EigenEqMat;
        VectorXd Eigen_b_right_side;
        VectorXd solution;

    public:
        typedef std::vector<double> stdvec;
        sparse_solve(sp_mat LinEqs, vec b_eq, int num_rows, int num_cols);
        sparse_solve(SpMatrx LinEqs, VectorXd b_eq, int num_rows, int num_cols, string solver_name);
        void solve(); 
        vector<double> get_solution();
        //void solve(VectorXd& solution, string solver_type); 

};

#endif

