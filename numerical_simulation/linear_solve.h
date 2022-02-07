#include"numerical_solver.h"
#include<iostream>
#include<vector>


#ifndef LINEAR_SOLVE_H
#define LINEAR_SOLVE_H

using namespace std;



class linear_solve: public numerical_solver
{
    private:
        int n_rows, n_cols;
        
