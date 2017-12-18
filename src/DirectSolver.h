/*
* This file is part of SimpleNumericalLirary, an open-source cross-platform
* Numerical Library
* 2007  Rafael Sachetto
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
* Contact e-mail: Rafael Sachetto <rsachetto@gmail.com>
* Program URL   : http://lib-matrix.sourceforge.net/
*
*/

#ifndef SOLVER_H_
#define SOLVER_H_
#include "Matrix.h"
#include "Vector.h"
#include "cmath"
#include<iostream>

using namespace std;

//! Direct Solver.
/*!
Implements Direct linear systems solvers.
(Gauss Elimination, Pivot Gauss Elimination, LU Factorazation,
 Doolittle Factorization and Cholesky Factorization)
*/

class DirectSolver
{
public:

	//! Constructor.
	/*!
		Construct a empty solver;
	*/
	DirectSolver(){};

	//! Constructor.
	/*!
		 Construct a solver that solves the linear system Ax = b;
		\param A a Matrix that represents the system.
		\param b the right-hand side of the system.
	*/
	DirectSolver(Matrix &A, Vector &b);
	
	//! Destructor.
	/*!
		Destroy the objects;
	*/
	virtual ~DirectSolver();
	
	//! Set the system Matrix.
	/*!
		Set the Matrix that will represent the system
		\param A a Matrix that represents the system.
	*/
	void setMatrixA(Matrix &A);
	
	//! Set the right-hand side of the system.
	/*!
	 	Set the Vector that will represents right-hand side of the system.
		\param b the right-hand side of the system.
	*/
	void setVectorB(Vector &b);
	
	//! Solve the system
	/*!
	 	Solve the system using the Gaussian Elimination;
		\return the Vector containing the result
	*/
	Vector gaussElimination();
	
	//! Solve the system
	/*!
	 	Solve the system using the Gaussian Elimination with partial pivoting
		\return the Vector containing the result
	*/
	Vector pivotGaussElimination();
	
	//! Solve the system
	/*!
	 	Solve the system using the Gaussian Elimination with total pivoting
		\return the Vector containing the result
	*/	
	Vector totalPivotGaussElimination();
	
	//! Solve the system
	/*!
	 	Solve the system using the LU factorization with Gaussian Elimination;
		\return the Vector containing the result
	*/
	Vector LUGauss();
	
	//! Solve the system
	/*!
	 	Solve the system using the LU factorization with Gaussian Elimination and pivoting
		\return the Vector containing the result
	*/
	Vector LUpivotGauss();
	
	//! Solve the system
	/*!
	 	Solve the system using the LU factorization with Gaussian Elimination and pivoting
	 	and refine the result
	 	\param maxIterations maximum number of refinements
	 	\param tolerance error tolarance
		\return the Vector containing the result
	*/
	Vector refinedLUpivotGauss(int maxIterations, double tolerance);
	
	//! Solve the system
	/*!
	 	Solve the system using the Doolittle factorization
		\return the Vector containing the result
	*/
	Vector doolittleFactorization();
	
	//! Solve the system
	/*!
	 	Solve the system using the Cholesky factorization
		\return the Vector containing the result
	*/
	Vector choleskyFactorization();
	
	//! Print the system
	/*!
	 	Print a human readable linear system
	*/
	void printSystem();
	
protected:
	Matrix A;	
	Vector b;
	Vector retroSubstitution(Matrix A, Vector b);

};

#endif /*SOLVER_H_*/
