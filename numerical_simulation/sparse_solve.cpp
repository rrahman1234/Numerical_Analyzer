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


sparse_solve::sparse_solve(sp_mat LinEqs, vec b_eq, int n_rows, int n_cols): EqnMat(LinEqs), b_right_side(b_eq), num_rows(n_rows), num_cols(n_cols)
{
    cout << "Linear Solver Class: Armadillo Library" << endl;
    order = 1;
}

sparse_solve::sparse_solve(SpMatrx LinEqs, VectorXd b_eq, int n_rows, int n_cols, string solver_name): EigenEqMat(LinEqs), Eigen_b_right_side(b_eq), num_rows(n_rows), num_cols(n_cols), solver_type(solver_name)
{
    cout << "Linear Solver Class: Eigen Library" << endl;
    order = 2;
}


void sparse_solve::solve()
{
    cout << "Solving" << endl;
    cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
    cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;
    if (order == 1)
    { 
        cout << "Solving" << endl;
        EqnMat.print("EqnMat:");
        b_right_side.print("b_right:");
        vec EqSol = arma::spsolve(EqnMat, b_right_side);
        EqSol.print("Solution:");
        for (size_t i = 0; i < EqSol.n_rows; ++i) {
            double var = EqSol(i);
            solution_vector.push_back(var);
        }
    }
    else if ((order == 2) && (solver_type == "QR"))
    {
        Eigen::SparseQR<SpMatrx, Eigen::COLAMDOrdering<int>> solver(EigenEqMat);
        solver.compute(EigenEqMat);
        solution = solver.solve(Eigen_b_right_side);
    }
    else if ((order == 2) && (solver_type == "LU"))
    {
        bool is_square = true;
        if (num_rows != num_cols) {
                is_square = false;
                throw std::runtime_error("Not Square Matrix!");
        }

        if(is_square)
        {
            Eigen::SparseLU<SpMatrx, Eigen::COLAMDOrdering<int>> solver(EigenEqMat);
            solver.compute(EigenEqMat);
            solution = solver.solve(Eigen_b_right_side);
        }
    }
    else if ((order == 2) && (solver_type == "LDLT"))
    {
        Eigen::SimplicialLDLT<SpMatrx, Eigen::Lower, Eigen::NaturalOrdering<int>> solver(EigenEqMat);
       
        bool is_PosDef;
        is_PosDef = is_positive_semi_definite<SpMatrx, Eigen::SimplicialLDLT<SpMatrx, Eigen::Lower, Eigen::NaturalOrdering<int>>&> (EigenEqMat, solver);

        if(is_PosDef)
        {
            solver.compute(EigenEqMat);
            solution = solver.solve(Eigen_b_right_side);
        }
    }
    solution_vector.insert(solution_vector.end(), std::make_move_iterator(solution.data()), std::make_move_iterator(solution.data() + solution.size()));
}


vector<double> sparse_solve::get_solution()
{
    if (order == 1)
    {
        cout << "Solution: Linear Solver with Armadillo" << " " << endl;
    }
    else if (order == 2)
    {
        cout << "Solution: Linear Solver with Eigen: " << solver_type << " " << endl;
    }
    return solution_vector;
}

