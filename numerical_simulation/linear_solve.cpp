#include"numerical_solver.h"
#include "linear_solve.h"
#include<math.h>
#include<iostream>
#include<vector>
#include<utility>
#include "armadillo"
#include <string> 
#include <Eigen/Dense>


using namespace Eigen;
using namespace std;
using namespace arma;


linear_solve::linear_solve(mat LinEqs, vec b_eq, int eq_order, int n_rows, int n_cols): EqnMat(LinEqs), b_right_side(b_eq), order(eq_order), num_rows(n_rows), num_cols(n_cols)
{
    cout << "Linear Solveer Class: Armadillo Library" << endl;
}

linear_solve::linear_solve(MatrixXd LinEqs, VectorXd b_eq, int eq_order, int n_rows, int n_cols): EigenEqMat(LinEqs), Eigen_b_right_side(b_eq), order(eq_order), num_rows(n_rows), num_cols(n_cols)
{
    cout << "Linear Solveer Class: Eigen Library" << endl;
}

//vector<double> linear_solve::solve()
//{
//    cout << "Solving" << endl;
//    
//    EqnMat.print("EqnMat:");
//    b_right_side.print("b_right:");
//    vec EqSol = arma::solve(EqnMat, b_right_side);
//    EqSol.print("Solution:");
//
//    for (size_t i = 0; i < EqSol.n_rows; ++i) {
//        double var = EqSol(i);
//        solution.push_back(var);
//    }
//    
//    return solution;
//}
//
//
//VectorXd linear_solve::solve()
//{
//    cout << "Solving" << endl;
//    cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
//    cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;
//
//    VectorXd solution = EigenEqMat.lu().solve(Eigen_b_right_side);
//
//    return solution;
//}

void linear_solve::solve(vector<double>& solution)
{
    cout << "Solving" << endl;
    
    EqnMat.print("EqnMat:");
    b_right_side.print("b_right:");
    vec EqSol = arma::solve(EqnMat, b_right_side);
    EqSol.print("Solution:");

    for (size_t i = 0; i < EqSol.n_rows; ++i) {
        double var = EqSol(i);
        solution.push_back(var);
    }
    
    //return solution;
}


void linear_solve::solve(VectorXd& solution)
{
    cout << "Solving" << endl;
    cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
    cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;

    solution = EigenEqMat.lu().solve(Eigen_b_right_side);

    //return solution;
}


