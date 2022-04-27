#include"numerical_solver.h"
#include "linear_solve.h"
#include<math.h>
#include<iostream>
#include<vector>
#include<utility>
#include "armadillo"
#include <string> 
#include <Eigen/Dense>
#include <typeinfo>
#include <stdexcept>


using namespace Eigen;
using namespace std;
using namespace arma;


linear_solve::linear_solve(mat& LinEqs, vec& b_eq, int n_rows, int n_cols): EqnMat(LinEqs), b_right_side(b_eq),  num_rows(n_rows), num_cols(n_cols)
{
    cout << "Linear Solver Class: Armadillo Library" << endl;
    order = 1;
}

//linear_solve::linear_solve(MatrixXd LinEqs, VectorXd b_eq, int n_rows, int n_cols): EigenEqMat(LinEqs), Eigen_b_right_side(b_eq), num_rows(n_rows), num_cols(n_cols)
//{
//    cout << "Linear Solver Class: Eigen Library" << endl;
//    order = 2;
//}

linear_solve::linear_solve(MatrixXd& LinEqs, VectorXd& b_eq, int n_rows, int n_cols, string solver_name): EigenEqMat(LinEqs), Eigen_b_right_side(b_eq), num_rows(n_rows), num_cols(n_cols), solver_type(solver_name)
{
    cout << "Linear Solver Class: Eigen Library" << endl;
    order = 2;
}


void linear_solve::solve()
{
    if (order == 1)
    {
        cout << "Solving" << endl;
        
        EqnMat.print("EqnMat:");
        b_right_side.print("b_right:");
        vec EqSol = arma::solve(EqnMat, b_right_side);
        EqSol.print("Solution:");

        for (size_t i = 0; i < EqSol.n_rows; ++i) {
            double var = EqSol(i);
            solution_vector.push_back(var);
        }
    }
    else if ((order == 2) && (solver_type == "LU"))
    {
        cout << "Solving" << endl;
        cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
        cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;

        FullPivLU<MatrixXd> lu(EigenEqMat);
        //solution = EigenEqMat.lu().solve(Eigen_b_right_side);
        solution = lu.solve(Eigen_b_right_side);
        solution_vector.insert(solution_vector.end(), std::make_move_iterator(solution.data()), std::make_move_iterator(solution.data() + solution.size()));
    }
    else if ((order == 2) && (solver_type == "QR"))
    {
        cout << "Solving" << endl;
        cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
        cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;

        FullPivHouseholderQR<MatrixXd> qr(EigenEqMat);
        solution = qr.solve(Eigen_b_right_side);
        solution_vector.insert(solution_vector.end(), std::make_move_iterator(solution.data()), std::make_move_iterator(solution.data() + solution.size()));
    }
    else if ((order == 2) && (solver_type == "LLT"))
    {
        cout << "Solving" << endl;
        cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
        cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;

        bool is_PosDef = true;

        Eigen::LLT<Eigen::MatrixXd>  llt(EigenEqMat);
        //check for positive definite condition
        if (!EigenEqMat.isApprox(EigenEqMat.transpose()) || llt.info() == Eigen::NumericalIssue) {
            is_PosDef = false;
            throw std::runtime_error("Possibly non semi-positive definitie matrix!");
        }         

        if(is_PosDef)
        {
            solution = llt.solve(Eigen_b_right_side);
            solution_vector.insert(solution_vector.end(), std::make_move_iterator(solution.data()), std::make_move_iterator(solution.data() + solution.size()));
        }
    }
    else if ((order == 2) && (solver_type == "LDLT"))
    {
        cout << "Solving" << endl;
        cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
        cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;

        bool is_PosSemDef = true;;

        Eigen::LDLT<Eigen::MatrixXd>  ldlt(EigenEqMat);
        //check for positive definite condition
        if (!EigenEqMat.isApprox(EigenEqMat.transpose()) || ldlt.info() == Eigen::NumericalIssue || !ldlt.isNegative() || !ldlt.isPositive()) {
            is_PosSemDef = false;
            throw std::runtime_error("Not Positive or negative semidefinite!");
        }         

        if(is_PosSemDef)
        {
            solution = ldlt.solve(Eigen_b_right_side);
            solution_vector.insert(solution_vector.end(), std::make_move_iterator(solution.data()), std::make_move_iterator(solution.data() + solution.size()));
        }
    }
}


vector<double> linear_solve::get_solution()
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

