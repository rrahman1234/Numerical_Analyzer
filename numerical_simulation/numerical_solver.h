#include <iostream>
#include <vector>
#include <utility>
#include"armadillo"
#include <Eigen/Dense>


#ifndef NUMERICAL_SOLVER_H
#define NUMERICAL_SOLVER_H

using namespace std;

class numerical_solver
{
	protected:
		double x, y;
		int i;
		int order;
        vector< pair< pair <double, double>, int> > equation_components;
        vector<double> solution_vector;
        //vector<double> solution;

	public:
		numerical_solver();
		numerical_solver(double x, int order);
		void print_result(vector<double>& solution);
        void get_list(double x1, vector<double>& poly_coeff);
		double function(double x1, vector<double>& poly_coeff);
		vector<double>* function_poly_terms(double x1, vector<double>& poly_coeff);
		vector<double>* function_poly_diff_terms(double x1, vector<double>& poly_coeff);
        virtual vector<double> get_solution() = 0;
		virtual void solve() = 0;

        template<typename U, typename T>
        bool is_positive_semi_definite(U EqnMat, T obj)
        {
            bool is_PosDef = true;
            if (!EqnMat.isApprox(EqnMat.transpose()) || obj.info() == Eigen::NumericalIssue) {
                    is_PosDef = false;
                    throw std::runtime_error("Possibly non semi-positive definitie matrix!");
            }
            return is_PosDef;
        }

        template<typename U, typename T>
        bool is_positive_negative_semi_definite(U EqnMat, T obj)
        {
            bool is_PosSemDef = true;
            if ((!EqnMat.isApprox(EqnMat.transpose())) || (obj.info() == Eigen::NumericalIssue) || (obj.isPositive() == false)) {
            is_PosSemDef = false;
            throw std::runtime_error("Not Positive or negative semidefinite!");
            }         
            return is_PosSemDef;
        }

        template<typename U = int, typename T>
        double findSum(int i, T Mat)
        {
            double sum = 0;
            for(int j=0; j<Mat.rows(); j++)
            {
              if(i!=j)
                sum+=Mat(i,j);
            }
            return sum;
        }

        template<typename T>
        bool isApplicable(T Mat)
        {
           for(int i=0;i<Mat.rows();i++)
           {
                for(int j=0;j<Mat.cols();j++)
                {
                   if(fabs(Mat(i,i))>findSum(i, Mat))
                      break;
                   else  
                      return false;   
                }
            }
            return true;
        }

        template<typename U = int, typename T>
        bool pivotZero(int row_indx, T& Mat)
        {
            if(Mat(row_indx, row_indx) == 0)
            {
                return true;
            }
            else
                return false;
        }

        template<typename U = int, typename T>
        void exchangeRow(int row_indx, T& Mat)
        {
            int k = row_indx + 1;
            int max_indx = row_indx + 1;

            while(k < Mat.rows())
            {
                if(Mat(max_indx, row_indx) < Mat(k, row_indx))
                    max_indx  = k;
                k++;
            }

            for (int j = 0; j < Mat.cols(); j++)
            {
                double temp;
                temp = Mat(row_indx, j);
                Mat(row_indx, j) = Mat(max_indx, j);
                Mat(max_indx, j) = temp;
            }
        }

        template<typename U = int, typename T>
        void normalize(int row_indx, T& Mat)
        {
            double  normal = Mat(row_indx, row_indx);
            for(int j = 0; j < Mat.cols(); j++)
            {
                Mat(row_indx, j) /= normal;
            }
        }  


        template<typename T>
        bool isApplicable_gaussElimination(T Mat)
        {
           for(int i=0;i<Mat.rows();i++)
           {
                for(int j=0;j<Mat.cols();j++)
                {
                   if(fabs(Mat(i,i)) == 0.0)
                   {
                      return false;   
                      break;
                   }
                }
           }
        return true;
        }
};


#endif 
