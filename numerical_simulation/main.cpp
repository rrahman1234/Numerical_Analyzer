#include"numerical_solver.h"
#include"quadratic_solve.h"
#include"bisection_solve.h"
#include"linear_solve.h"
#include"sparse_solve.h"
#include"newton_raphson_solve.h"
#include<iostream>
#include<vector>
#include"armadillo"
#include <Eigen/Dense>


using namespace std;
using namespace arma;
using namespace Eigen;

int main()
{
	double X = 3.0;
	int order = 0;

    //Quadratic Solve
	//double coff[] = {1.0, 2.0};
	//vector<double> quad;
    //quadratic_solve quad_Sol(X,order,coff);
    //quad_Sol.solve();
    //quad = quad_Sol.get_solution();
	//quad_Sol.print_result(quad);
    
    //Newton Raphson
	//vector<double> nr;
	//vector<double> coff;
	//vector<double> diff_coff;
	//coff.push_back(3.0);
	//coff.push_back(0.0);
	//coff.push_back(-2.0);
	//coff.push_back(1.0);
	//diff_coff.push_back(3.0);
	//diff_coff.push_back(2.0);
	//newton_raphson_solve NR_Sol(X, 2, 2.0, 50, 0.001, coff, diff_coff);		
    //NR_Sol.solve();
	//nr = NR_Sol.get_solution();
	//NR_Sol.print_result(nr);
    

    //Bisection
	//vector<double> bisec;
	//vector<double> coff;
	//coff.push_back(3.0);
	//coff.push_back(0.0);
	//coff.push_back(-2.0);
	//coff.push_back(1.0);
	//bisection_solve BSolv(-10, 20, 3, 100, 0.0001, coff); 
    //BSolv.solve();
	//bisec = BSolv.get_solution();
	//BSolv.print_result(bisec);


    cout << "*********" << endl;
    
    int n_rows = 2;
    int n_cols = 2;
	mat Eq_arma(n_rows, n_cols);
    Eq_arma.fill(0.0);
    Eq_arma = {{3, 2}, {7, 1}};
    vec b_arma(n_rows);
    b_arma.fill(0.0);
    b_arma = {16, 19};
    linear_solve LinSolv_arma(Eq_arma, b_arma, n_rows, n_cols);
    vector<double> LinSolvSol_arma; 
    LinSolv_arma.solve(); 
    LinSolvSol_arma = LinSolv_arma.get_solution();
    LinSolv_arma.print_result(LinSolvSol_arma);

    cout << "*********" << endl;

	MatrixXd Eq(n_rows, n_cols);
    VectorXd b(n_rows);
    Eq << 3, 2, 
          7, 1;
    b << 16, 19;
    cout << Eq << endl;
    cout << b<< endl;
    cout << "LU decomposition" << endl;
    linear_solve LinSolv_Eigen(Eq, b, n_rows, n_cols, "LU");
    vector<double> LinSolvSol_Eigen;
    LinSolv_Eigen.solve();
    LinSolvSol_Eigen = LinSolv_Eigen.get_solution();
    LinSolv_Eigen.print_result(LinSolvSol_Eigen);

    cout << "QR decomposition" << endl;
    linear_solve LinSolv_Eigen_2(Eq, b, n_rows, n_cols, "QR");
    vector<double> LinSolvSol_Eigen_2;
    LinSolv_Eigen_2.solve();
    LinSolvSol_Eigen_2 = LinSolv_Eigen_2.get_solution();
    LinSolv_Eigen_2.print_result(LinSolvSol_Eigen_2);
    
    cout << "LDLT decomposition" << endl;
    linear_solve LinSolv_Eigen_3(Eq, b, n_rows, n_cols, "LDLT");
    vector<double> LinSolvSol_Eigen_3;
    LinSolv_Eigen_3.solve();
    LinSolvSol_Eigen_3 = LinSolv_Eigen_3.get_solution();
    LinSolv_Eigen_3.print_result(LinSolvSol_Eigen_3);

    //cout << "*********" << endl;

    //int n_rows_sp = 1000;
    //int n_cols_sp = 1000;
	//sp_mat Eq_arma_sp = sprandu<sp_mat>(n_rows_sp, n_cols_sp, 0.1);
    //vec b_arma_sp(n_rows_sp);
    //b_arma_sp.fill(1.0);
    //sparse_solve LinSolv_arma_sp(Eq_arma_sp, b_arma_sp, n_rows_sp, n_cols_sp);
    //vector<double> LinSolvSol_arma_sp; 
    //LinSolv_arma_sp.solve();
    //LinSolvSol_arma_sp = LinSolv_arma_sp.get_solution();
    //LinSolv_arma_sp.print_result(LinSolvSol_arma_sp);

    //cout << "*********" << endl;
	//
    //SpMatrx Eq_2(5, 5);
    //VectorXd b_2(5);
    //
    //Eq_2.insert(0,1) = 3;
    //Eq_2.insert(1,0) = 22;
    //Eq_2.insert(1,4) = 17;
    //Eq_2.insert(2,0) = 17;
    //Eq_2.insert(2,1) = 5;
    //Eq_2.insert(2,3) = 1;
    //Eq_2.insert(4,2) = 14;
    //Eq_2.insert(3,4) = 8;

    //b_2 << 16, 19, 20, 12, 1;
    //cout << Eq_2 << endl;
    //cout << b_2 << endl;
    //linear_solve LinSolv_2(Eq_2, b_2, n_rows, n_cols, "LU");
    //vector<double> LinSolvSol_2;
    //LinSolv_2.solve();
    //LinSolvSol_2 = LinSolv_2.get_solution();
    //LinSolv_2.print_result(LinSolvSol_2);
    //
    //cout << "*********" << endl;
    

   
	return 0;
}
