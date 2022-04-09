#include"numerical_solver.h"
#include "sparse_solve.h"
#include<math.h>
#include<iostream>
#include<vector>
#include<utility>
#include "armadillo"
#include <string> 
#include <Eigen/Sparse>


using namespace Eigen;
using namespace std;
using namespace arma;


sparse_solve::sparse_solve(sp_mat LinEqs, vec b_eq, int eq_order, int n_rows, int n_cols): EqnMat(LinEqs), b_right_side(b_eq), order(eq_order), num_rows(n_rows), num_cols(n_cols)
{
    cout << "Linear Solver Class: Armadillo Library" << endl;
}

sparse_solve::sparse_solve(SpMatrx LinEqs, VectorXd b_eq, int eq_order, int n_rows, int n_cols): EigenEqMat(LinEqs), Eigen_b_right_side(b_eq), order(eq_order), num_rows(n_rows), num_cols(n_cols)
{
    cout << "Linear Solver Class: Eigen Library" << endl;
}



void sparse_solve::solve(vector<double>& solution)
{
    cout << "Solving" << endl;
    
    EqnMat.print("EqnMat:");
    b_right_side.print("b_right:");
    vec EqSol = arma::spsolve(EqnMat, b_right_side);
    EqSol.print("Solution:");

    for (size_t i = 0; i < EqSol.n_rows; ++i) {
        double var = EqSol(i);
        solution.push_back(var);
    }
}


void sparse_solve::solve(VectorXd& solution)
{
    cout << "Solving" << endl;
    cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
    cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;

    Eigen::SimplicialCholesky<SpMatrx> chol(EigenEqMat);
    solution = chol.solve(Eigen_b_right_side);
}



