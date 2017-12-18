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

#ifndef ITERATIVESOLVER_H_
#define ITERATIVESOLVER_H_
#include "Matrix.h"
#include "Vector.h"
#include "cmath"
#include<iostream>

using namespace std;

//! Iterative Solver.
/*!
Implements Iterative linear systems solvers.
(Jacob, Gauss-Seidel, SOR, Gradient Method, Conjugated Gradient Method)
*/

class IterativeSolver
{
public:
	//! Constructor.
	/*!
		Construct a empty solver;
	*/
	IterativeSolver();
	
	//! Constructor.
	/*!
		 Construct a solver that solves the linear system Ax = b;
		\param A a Matrix that represents the system.
		\param b the right-hand side of the system.
	*/
	IterativeSolver(Matrix &A, Vector &b);
	
	//! Destructor.
	/*!
		Destroy the objects;
	*/
	virtual ~IterativeSolver();
	
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
	 	Solve the system using the Jacob Method
	 	\param x_0 Initial guess
	 	\param maxIteration maximum number of iterations
	 	\param tolerance error tolerance
		\return the Vector containing the result
	*/
	Vector jacob(Vector x_0, int maxIteration, double tolerance);
	
	//! Solve the system
	/*!
	 	Solve the system using the Jacob Method (Initial guess setted to 0)
	 	\param maxIteration maximum number of iterations
	 	\param tolerance error tolerance
		\return the Vector containing the result
	*/
	Vector jacob(int maxIteration, double tolerance);
	
	//! Solve the system
	/*!
	 	Solve the system using the Gauss-Seidel Method
	 	\param x_0 Initial guess
	 	\param maxIteration maximum number of iterations
	 	\param tolerance error tolerance
		\return the Vector containing the result
	*/
	Vector gaussSeidel(Vector x_0, int maxIteration, double tolerance);
	
	//! Solve the system
	/*!
	 	Solve the system using the Gauss-Seidel Method (Initial guess setted to 0)
	 	\param maxIteration maximum number of iterations
	 	\param tolerance error tolerance
		\return the Vector containing the result
	*/
	Vector gaussSeidel(int maxIteration, double tolerance);
	
	//! Solve the system
	/*!
	 	Solve the system using the SOR Method
	 	\param x_0 Initial guess
	 	\param maxIteration maximum number of iterations
	 	\param tolerance error tolerance
	 	\param omega SOR omega parameter
		\return the Vector containing the result
	*/
	Vector SOR(Vector x_0, int maxIteration, double tolerance, double omega);
	
	//! Solve the system
	/*!
	 	Solve the system using the SOR Method (Initial guess setted to 0)
	 	\param maxIteration maximum number of iterations
	 	\param tolerance error tolerance
	 	\param omega SOR omega parameter
		\return the Vector containing the result
	*/
	Vector SOR(int maxIteration, double tolerance, double omega);
	
	//! Solve the system
	/*!
	 	Solve the system using the Gradient Method
	 	\param v_0 Initial guess
	 	\param maxIteration maximum number of iterations
	 	\param tolerance error tolerance
		\return the Vector containing the result
	*/
	Vector gradientMethod(Vector v_0, int maxIteration, double tolerance);
	
	//! Solve the system
	/*!
	 	Solve the system using the Gradient Method (Initial guess setted to 0)
	 	\param maxIteration maximum number of iterations
	 	\param tolerance error tolerance
		\return the Vector containing the result
	*/
	Vector gradientMethod(int maxIteration, double tolerance);
	
	//! Solve the system
	/*!
	 	Solve the system using the Conjugated Gradient Method
	 	\param v_0 Initial guess
	 	\param maxIteration maximum number of iterations
	 	\param tolerance error tolerance
		\return the Vector containing the result
	*/
	Vector conjugatedGradientMethod(Vector v_0, int maxIteration, double tolerance);
	
	//! Solve the system
	/*!
	 	Solve the system using the Conjugated Gradient Method (Initial guess setted to 0)
	 	\param maxIteration maximum number of iterations
	 	\param tolerance error tolerance
		\return the Vector containing the result
	*/
	Vector conjugatedGradientMethod(int maxIteration, double tolerance);
	Vector conjugatedGradientMethod2(Vector v_0, int maxIteration, double tolerance);
	
	//! Print the system
	/*!
	 	Print a human readable linear system
	*/
	void printSystem();
	
protected:
	Matrix A;	
	Vector b;
	double stopCriteriaJacob(Vector x_0, Vector x_1);
	bool stopCriteriaGradient(Vector v_0, Vector v_1, Vector r_1, double tolerance);
};

#endif /*ITERATIVESOLVER_H_*/
