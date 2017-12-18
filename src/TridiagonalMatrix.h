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

#ifndef TRIDIAGONALMATRIX_H_
#define TRIDIAGONALMATRIX_H_
#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include <string>
#include "MatrixFileHandler.h"
#include "Vector.h"


using namespace std;

//! An nxn Tridiagonal Matrix.
/*!
A class that represents a Tridiagonal Matrix of doubles
*/
class TridiagonalMatrix
{
public:

	//! Constructor.
	/*!
		Construct a empty matrix;
	*/
	TridiagonalMatrix();
	
	//! Constructor.
	/*!
		 Construct a matrix with the gived diagonal size
		\param diagonalSize the principal diagonal size of the matrix
	*/
	TridiagonalMatrix(int diagonalSize);
	

	//! Destructor.
	/*!
		Destroy the objects;	*/
	
	virtual ~TridiagonalMatrix();	

	//! Constructor.
	/*!
		 Construct a copy of the matrix m
		\param m the matrix to be copied
	*/
	TridiagonalMatrix(const TridiagonalMatrix &m);
	
	//! Get Method
	/*!
	 	Get a element of the matrix 
	 	\param i the line of the element
	 	\param j the column of the element
		\return the element located at line i and column j if 
		the elements belongs a one of the diagonal, 0 otherwise
	*/
	double getIJ( int i, int j);
	
	//! Set Method
	/*!
	 	Set a element of the matrix if this elemente belongs to one of tree diagonals 
		\param i the line to put the element
		\param j the column to put the element
		\param value the value to set
	*/
	void setIJ( int i, int j, double value);
	
	//! Overloaded operator =
	/*!
	 	Attributes a matrix to another
	 	\param rhs the matrix to be attributed
	*/	
	const TridiagonalMatrix& operator = (const TridiagonalMatrix &rhs);
	
	//! Overloaded operator ()
	/*!
	 	Get a element of the matrix 
	 	\param i the line of the element
	 	\param j the column of the element
		\return the element located at line i and column j if 
		the elements belongs a one of tree diagonals, 0 otherwise
		
		ex: double value = matrix(1,1);
	*/
	double operator () (int i, int j) const;
	
	
	double getUpperDiagonalValue(int index);
	double getPrincipalDiagonalValue(int index);
	double getLowerDiagonalValue(int index);
	void printMatrix();
	

private:
	double *lowerDiagonal, *principalDiagonal, *upperDiagonal;
	int diagonalSize;
	void alocateMatrix(int diagonalSize);
	void copyMatrix(double *destLower,double *destPrinc,double *destUpper , double *srcLower, double *srcPrinc, double *srcUpper);
	void freeDiagonal(double *diagonal);
	
};

#endif /*TRIDIAGONALMATRIX_H_*/
