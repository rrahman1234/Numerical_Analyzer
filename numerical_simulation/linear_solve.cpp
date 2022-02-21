#include"numerical_solver.h"
#include "linear_solve.h"
#include<math.h>
#include<iostream>
#include<vector>
#include<utility>
#include "armadillo"

using namespace std;
using namespace arma;


linear_solve::linear_solve(mat LinEqs, int eq_order, int n_rows, int n_cols):A(LinEqs), order(eq_order), num_rows(n_rows), num_cols(n_cols)
{
    cout << "Linear Solveer Class" << endl;
}

vector<double> linear_solve::solve()
{
    cout << "Solving" << endl;
    vector<double> v = {1.0};

    return v;
}
