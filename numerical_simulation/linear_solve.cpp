#include"numerical_solver.h"
#include "linear_solve.h"
#include<math.h>
#include<iostream>
#include<vector>
#include<utility>
#include "armadillo"

using namespace std;
using namespace arma;


linear_solve::linear_solve(mat LinEqs, vec b_eq, int eq_order, int n_rows, int n_cols): EqnMat(LinEqs), b_right_side(b_eq), order(eq_order), num_rows(n_rows), num_cols(n_cols)
{
    cout << "Linear Solveer Class" << endl;
}

vector<double> linear_solve::solve()
{
    cout << "Solving" << endl;
   
    EqnMat.print("EqnMat:");
    b_right_side.print("b_right:");

    vec EqSol = arma::solve(EqnMat, b_right_side);
    EqSol.print("Solution:");

    //vector<double> LinSolution(num_rows);

    for (size_t i = 0; i < EqSol.n_rows; ++i) {
        double var = EqSol(i);
        solution.push_back(var);
    }
    

    return solution;
}
