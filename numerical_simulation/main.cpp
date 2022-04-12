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
	int order = 1;

    //Quadratic Solve
	//double coff[] = {1.0, 2.0};
	//quadratic_solve quad_Sol(X,order,coff);
	//quad_Sol.solve();
	//quad_Sol.print_result();
    
    //int n_rows = 2;
    //int n_cols = 2;
	//mat Eq_arma(n_rows, n_cols);
    //Eq_arma.fill(0.0);
    //Eq_arma = {{3, 2}, {7, 1}};
    //vec b_arma(n_rows);
    //b_arma.fill(0.0);
    //b_arma = {16, 19};
    //linear_solve LinSolv_arma(Eq_arma, b_arma, order, n_rows, n_cols);
    //vector<double> LinSolvSol_arma; 
    //LinSolv_arma.solve(LinSolvSol_arma);
    //LinSolv_arma.print_result();

    //cout << "*********" << endl;

    //////int n_rows = 2;
    //////int n_cols = 2;
	//MatrixXd Eq(n_rows, n_cols);
    //VectorXd b(n_rows);
    //Eq << 3, 2, 
    //      7, 1;
    //b << 16, 19;
    //cout << Eq << endl;
    //cout << b<< endl;
    //linear_solve LinSolv(Eq, b, order, n_rows, n_cols);
    //VectorXd LinSolvSol;
    //LinSolv.solve(LinSolvSol);
    //cout << "Solution:" << endl << LinSolvSol << endl;

    cout << "*********" << endl;

    int n_rows = 1000;
    int n_cols = 1000;
	sp_mat Eq_arma = sprandu<sp_mat>(1000, 1000, 0.1);
    vec b_arma(n_rows);
    b_arma.fill(1.0);
    sparse_solve LinSolv_arma(Eq_arma, b_arma, order, n_rows, n_cols);
    vector<double> LinSolvSol_arma; 
    LinSolv_arma.solve(LinSolvSol_arma);
    LinSolv_arma.print_result();

    cout << "*********" << endl;
	
    SpMatrx Eq(5, 5);
    VectorXd b(5);
    
    Eq.insert(0,1) = 3;
    Eq.insert(1,0) = 22;
    Eq.insert(1,4) = 17;
    Eq.insert(2,0) = 17;
    Eq.insert(2,1) = 5;
    Eq.insert(2,3) = 1;
    Eq.insert(4,2) = 14;
    Eq.insert(3,4) = 8;

    b << 16, 19, 20, 12, 1;
    cout << Eq << endl;
    cout << b<< endl;
    linear_solve LinSolv(Eq, b, order, n_rows, n_cols);
    VectorXd LinSolvSol;
    LinSolv.solve(LinSolvSol);
    cout << "Solution:" << endl << LinSolvSol << endl;

    cout << "*********" << endl;
    
    sparse_solve SprsSolv(Eq, b, order, n_rows, n_cols); 
    VectorXd SprsSolvSol;
    SprsSolv.solve(SprsSolvSol, "QR");
    cout << "Solution:" << endl << SprsSolvSol << endl;

   
    //////Newton Raphson
	//vector<double> coff;
	//////vector<double> diff_coff;
	//coff.push_back(3.0);
	//coff.push_back(0.0);
	//coff.push_back(-2.0);
	//coff.push_back(1.0);
    ////
	//////diff_coff.push_back(3.0);
	//////diff_coff.push_back(2.0);

	//////newton_raphson_solve NR_Sol(X, 2, 2.0, 50, 0.001, coff, diff_coff);		
	//////NR_Sol.solve();
	//////NR_Sol.function(X, coff);

	//bisection_solve BSolv(-10, 20, 3, 100, 0.0001, coff); 
	//BSolv.solve();
    ////BSolv.get_list(2, coff);
    ////BSolv.get_list(3, coff);

    

	return 0;
}
