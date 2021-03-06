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


using namespace std;
using namespace arma;
using namespace Eigen;

linear_solve::linear_solve(arma::mat& LinEqs, arma::vec& b_eq, int n_rows, int n_cols): EqnMat(LinEqs), b_right_side(b_eq),  num_rows(n_rows), num_cols(n_cols)
{
    cout << "Linear Solver Class: Armadillo Library" << endl;
    order = 1;
}

linear_solve::linear_solve(Eigen::MatrixXd& LinEqs, Eigen::VectorXd& b_eq, int n_rows, int n_cols, string solver_name, int num_iter, std::vector<double>* solution_init, double error_tolerance): EigenEqMat(LinEqs), Eigen_b_right_side(b_eq), num_rows(n_rows), num_cols(n_cols), solver_type(solver_name), num_iterations(num_iter), solution_initial(solution_init), err_tol(error_tolerance) 
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

        bool is_PosDef;

        Eigen::LLT<Eigen::MatrixXd>  llt(EigenEqMat);
        is_PosDef = is_positive_semi_definite<Eigen::MatrixXd, Eigen::LLT<Eigen::MatrixXd>> (EigenEqMat, llt);

        if(is_PosDef)
        {
            solution = llt.solve(Eigen_b_right_side);
            solution_vector.insert(solution_vector.end(), std::make_move_iterator(solution.data()), std::make_move_iterator(solution.data() + solution.size()));
        }
    }
    else if ((order == 2) && (solver_type == "LLT_Upper"))
    {
        cout << "Solving" << endl;
        cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
        cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;

        bool is_PosDef = true;

        Eigen::LLT<Eigen::MatrixXd, Eigen::Upper>  llt(EigenEqMat);
        is_PosDef = is_positive_semi_definite<Eigen::MatrixXd, Eigen::LLT<Eigen::MatrixXd, Eigen::Upper>&>(EigenEqMat, llt);
                 
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

        bool is_PosSemDef;

        Eigen::LDLT<Eigen::MatrixXd>  ldlt(EigenEqMat);
        
        if (ldlt.isPositive() == true)
        {
            std::cout << "Positive" << endl;
        }
        else if (ldlt.isNegative() == true)
        {
            std::cout << "Negative" << endl;
        }

        is_PosSemDef = is_positive_negative_semi_definite<Eigen::MatrixXd, Eigen::LDLT<Eigen::MatrixXd>&>(EigenEqMat, ldlt);

        if(is_PosSemDef)
        {
            solution = ldlt.solve(Eigen_b_right_side);
            solution_vector.insert(solution_vector.end(), std::make_move_iterator(solution.data()), std::make_move_iterator(solution.data() + solution.size()));
        }
    }
    else if ((order == 2) && (solver_type == "LDLT_Upper"))
    {
        cout << "Solving" << endl;
        cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
        cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;
        bool is_PosSemDef = true;

        Eigen::LDLT<Eigen::MatrixXd, Eigen::Upper>  ldlt(EigenEqMat);
        is_PosSemDef = is_positive_negative_semi_definite<Eigen::MatrixXd, Eigen::LDLT<Eigen::MatrixXd, Eigen::Upper>&>(EigenEqMat, ldlt);

        if(is_PosSemDef)
        {
            solution = ldlt.solve(Eigen_b_right_side);
            solution_vector.insert(solution_vector.end(), std::make_move_iterator(solution.data()), std::make_move_iterator(solution.data() + solution.size()));
        }
    }
    else if ((order == 2) && (solver_type == "SVD"))
    {
        cout << "Solving" << endl;
        cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
        cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;

        Eigen::JacobiSVD<Eigen::MatrixXd, ComputeThinU | ComputeThinV> jacob_svd(EigenEqMat);
        solution = jacob_svd.solve(Eigen_b_right_side);
        solution_vector.insert(solution_vector.end(), std::make_move_iterator(solution.data()), std::make_move_iterator(solution.data() + solution.size()));
    }
    else if ((order == 2) && (solver_type == "Gauss-Siedel"))
    {
        cout << "Solving" << endl;
        cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
        cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;
       
        int iter_count = 0;

        solution.resize(num_rows, 1);
        Eigen::VectorXd solution_x0 = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(solution_initial->data(), solution_initial->size());
      
        int N_mat = EigenEqMat.rows();
        bool converged_done = false;
        bool not_applicable = false;

        if(isApplicable<Eigen::MatrixXd>(EigenEqMat))
        {
            while(num_iterations > 0)
            {
                for(int i=0; i<EigenEqMat.rows(); i++) {
                    solution(i) = Eigen_b_right_side(i)/EigenEqMat(i,i);
                    for(int j=0; j<EigenEqMat.cols(); j++)
                    {
                        if(j == i) continue;
                        solution(i) = solution(i) - (EigenEqMat(i, j)/EigenEqMat(i, i))*solution_x0(j);
                    }
                    if(fabs(solution_x0(i)-solution(i)) < err_tol)
                    {
                        converged_done = true;
                    }
                    solution_x0(i) = solution(i);
                }

            if(converged_done == true)
            {
               cout << "Solution converged" << endl;
                break;
            }
            num_iterations--;
            iter_count++;
            cout << "End of Iteration#: " << iter_count << endl;
            }
        }
        else
        {
            cout << "Gauss-Siedel not applicable" << endl;
            not_applicable = true;
        }

        if(converged_done && !not_applicable)
        {
            solution_vector.insert(solution_vector.end(), std::make_move_iterator(solution.data()), std::make_move_iterator(solution.data() + solution.size()));
        }
        else
        {
            cout << "------------------------ Either Not Converged or Gauss-Siedel not Applicable ------------------------" << endl;
            abort();
        }
    }
    else if ((order == 2) && (solver_type == "Gauss-Jordan"))
    {
        cout << "Solving" << endl;
        cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
        cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;
        
        solution.resize(num_rows, 1);

        cout << "****" << endl;
        cout << EigenEqMat << endl;
      
        //If Matrix is invertble
        bool Is_invertible = EigenEqMat.fullPivLu().isInvertible();

        if(Is_invertible)
        {
            //Combine A and b (from Ax=b) 
            Eigen::MatrixXd EigenEqMat_Comb(EigenEqMat.rows(), EigenEqMat.cols()+Eigen_b_right_side.cols());
            EigenEqMat_Comb << EigenEqMat, Eigen_b_right_side;

            cout << "**** After Combining ****" << endl;
            cout << EigenEqMat_Comb << endl;

            //normalize and exchange rows
            for(int i=0; i < EigenEqMat_Comb.rows(); i++)
            {
                if(!pivotZero<int, Eigen::MatrixXd&>(i, EigenEqMat_Comb))
                {
                    normalize<int, Eigen::MatrixXd&>(i, EigenEqMat_Comb);
                }
                else
                {
                    exchangeRow<int, Eigen::MatrixXd&>(i, EigenEqMat_Comb);
                    normalize<int, Eigen::MatrixXd&>(i, EigenEqMat_Comb);
                }

                for(int k = i+1; k < EigenEqMat_Comb.rows(); k++)
                {
                    double zeroFactor = EigenEqMat_Comb(k, i);
                    for(int j = 0; j < EigenEqMat_Comb.cols(); j++)
                    {
                        EigenEqMat_Comb(k, j) -= zeroFactor*EigenEqMat_Comb(i, j);
                    }
                }
            }

            cout << "**** Combined A and b ****" << endl;
            cout << EigenEqMat_Comb << endl;

            //Put back EigenEqMat and Eigen_b_right_side
            EigenEqMat = EigenEqMat_Comb.block(0, 0 ,EigenEqMat.rows(), EigenEqMat.cols());
            Eigen_b_right_side = EigenEqMat_Comb.block(0, EigenEqMat_Comb.cols()-1, EigenEqMat_Comb.rows(), Eigen_b_right_side.cols());

            //Back-substitution
            for (int i = EigenEqMat.rows()-1; i >= 0; i--)
            {
                solution(i) = Eigen_b_right_side(i);   //EigenEqMat(i, EigenEqMat.cols()-1);
                for (int j=i+1; j<EigenEqMat.cols(); j++)
                {
                    solution(i) -= EigenEqMat(i,j)*solution(j);
                }
            }
            solution_vector.insert(solution_vector.end(), std::make_move_iterator(solution.data()), std::make_move_iterator(solution.data() + solution.size()));
        }
        else
        {
            cout << "------------------------ Either Not Converged or Gauss-Jordan method not Applicable ------------------------" << endl;
            abort();
        }
    }
    else if ((order == 2) && (solver_type == "Gauss-Elimination"))
    {
        cout << "Solving" << endl;
        cout << "EigenEqnMat:" << endl << EigenEqMat << endl;
        cout << "Right side of the equation" << endl << Eigen_b_right_side << endl;
        solution.resize(num_rows, 1);
        
        int N_mat = EigenEqMat.rows();
        bool converged_done = false;
        bool not_applicable = false;
        
        int i, j;

        if(!isApplicable_gaussElimination<Eigen::MatrixXd>(EigenEqMat))
        {
            cout << "Not Applicable: Diagonal element can not be zero" << endl;
            abort();
        }
        else
        {
            for(j=0; j < EigenEqMat.cols(); j++)
            {
                for(i=j+1; i < EigenEqMat.rows(); i++)
                {
                    double ratio = EigenEqMat(i,j)/EigenEqMat(j,j);
                    for(int k=0; k < EigenEqMat.cols(); k++)
                    {
                        EigenEqMat(i,k) = EigenEqMat(i,k) - ratio*EigenEqMat(j,k);
                    }
                    Eigen_b_right_side(i) = Eigen_b_right_side(i) - ratio*Eigen_b_right_side(j);
                }
            }
        }
        

        //Solution: Back Substitution
        solution(EigenEqMat.rows()-1) = Eigen_b_right_side(EigenEqMat.rows()-1)/EigenEqMat(EigenEqMat.rows()-1, EigenEqMat.cols()-1); 
        for(int i=EigenEqMat.rows()-2; i >= 0; i--)
        {
            int sum = 0;
            
            for(int j=i+1; j < EigenEqMat.cols(); j++)
            {
                sum = sum + EigenEqMat(i,j)*solution(j);
            }
            solution(i) = (Eigen_b_right_side(i) - sum)/EigenEqMat(i,i);
        }
        solution_vector.insert(solution_vector.end(), std::make_move_iterator(solution.data()), std::make_move_iterator(solution.data() + solution.size()));
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

