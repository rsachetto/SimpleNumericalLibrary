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

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include <string>
#include <iomanip>
#include "MatrixFileHandler.h"
#include "Vector.h"

using namespace std;

class MatrixException : public runtime_error
{
    public:
        MatrixException(const string &msg) : runtime_error(msg){}
};

//! An nxn Matrix.
/*!
A class that represents a nxn Matrix of doubles
*/
class Matrix
{
public:
	//! Constructor.
	/*!
		Construct a empty matrix;
	*/
	Matrix() {matrix = NULL;}
	
	//! Constructor.
	/*!
		 Construct a matrix numLinesXnumCols
		\param numLines lines number of the matrix
		\param numCols columns number of the matrix
	*/
	Matrix( int numLines , int numCols );
	
	//! Constructor.
	/*!
		 Construct a copy of the matrix m
		\param m the matrix to be copied
	*/
	Matrix(const Matrix &m);
	
	//! Constructor.
	/*!
		 Construct a matrix with the values in mat
		\param mat source two dimensional array
		\param numLines lines number of the matrix
		\param numCols columns number of the matrix
		
	*/
	Matrix(double **mat, int numLines , int numCols );
	
	//! Destructor.
	/*!
		Destroy the objects;
	*/
	virtual ~Matrix();
	
	//! Get Method
	/*!
	 	Get the number of lines of the matrix 
		\return the number of lines of the matrix
	*/
	int getNumLines() const;
	
	//! Get Method
	/*!
	 	Get the number of columns of the matrix 
		\return the number of columns of the matrix
	*/
	int getNumCols() const;
	
	//! Get Method
	/*!
	 	Get a element of the matrix 
	 	\param i the line of the element
	 	\param j the column of the element
		\return the element located at line i and column j
	*/
	double getIJ( int i, int j) const;
	
	//! Set Method
	/*!
	 	Set a element of the matrix 
		\param i the line to put the element
		\param j the column to put the element
		\param value the value to set
	*/
	void setIJ( int i, int j, double value);
	
	//! Print the matrix
	/*!
	 	Print the matrix in a nxn form
	*/	
	void printMatrix() const;
	
	//! Print the matrix
	/*!
	 	Print the matrix in a nxn form with a precision
	 	\param precision the precision to print
	*/
	void printMatrix(int precision) const;
	
	//! Overloaded operator =
	/*!
	 	Attributes a matrix to another
	 	\param rhs the matrix to be attributed
	*/
	const Matrix& operator = (const Matrix &rhs);
	
	//! Overloaded operator +
	/*!
	 	Sums two matrices
	 	\param rhs the right-hand side of the operation
	 	\return the result of operation
	*/
	Matrix operator + (const Matrix &rhs) const;
	
	//! Overloaded operator -
	/*!
	 	Substracts two matrices
	 	\param rhs the right-hand side of the operation
	 	\return the result of operation
	*/
	Matrix operator - (const Matrix &rhs) const;
	
	//! Overloaded operator *
	/*!
	 	Multiplies two matrices
	 	\param rhs the right-hand side of the operation
	 	\return the result of operation
	*/
	Matrix operator * (const Matrix &rhs) const;
	
	//! Overloaded operator *
	/*!
	 	Multiplies a matrix by a vector
	 	\param rhs the right-hand side of the operation
	 	\return the result of operation
	*/
	Vector operator * (Vector &rhs) const;
	
	//! Overloaded operator *
	/*!
	 	Multiplies a matrix by a scalar
	 	\param escalar the right-hand side of the operation
	 	\return the result of operation
	*/
	Matrix operator * (double escalar);
	
	//! Overloaded operator /
	/*!
	 	Divides a matrix by a scalar
	 	\param escalar the right-hand side of the operation
	 	\return the result of operation
	*/
	Matrix operator / (double escalar);
	
	//! Overloaded operator ()
	/*!
	 	Set a element of the matrix 
		\param i the line to put the element
		\param j the column to put the element
		
		ex: matrix(1,1) = 5.0;
	*/
	double &operator () (int i, int j);
	
	//! Overloaded operator ()
	/*!
	 	Get a element of the matrix 
	 	\param i the line of the element
	 	\param j the column of the element
		\return the element located at line i and column j
		
		ex: double value = matrix(1,1);
	*/
	double operator () (int i, int j) const;
	
	//! Verify the allocation of the matrix
	/*!
	 	Verify if the matrix was allocated	 	
		\return true if the matrix was allocates, false otherwise
	*/
	bool isNull();
	
	//! Get the transposed matrix
	/*!
	 	Calculate the transposed matrix	 	
		\return the transposed matrix		
	*/
	Matrix getTransposed() const;
	
	//! Get the inverse matrix
	/*!
	 	Calculate the inverse matrix	 	
		\return the inverse matrix		
	*/
	Matrix getInverseMatrix();
	
	//! Get the principal diagonal
	/*!
	 	Get the principal diagonal of the matrix	 	
		\return the principal diagonal of the matrix		
	*/	
	Vector getPrincipalDiagonal() const;
	
	//! Get the secundary diagonal
	/*!
	 	Get the secundary diagonal of the matrix	 	
		\return the secundary diagonal of the matrix		
	*/	
	Vector getSecundaryDiagonal() const;
	
	//! Get a line of the matrix
	/*!
	 	Get a specific line of the matrix
	 	\param line line number to get
		\return a line of the matrix
	*/	
	Vector getLine(int line) const;
	
	//! Get a column of the matrix
	/*!
	 	Get a specific column of the matrix
	 	\param column column number to get
		\return a column of the matrix
	*/	
	Vector getColumn(int column) const;	
	
	//! Get the determinant of the matrix
	/*!
	 	Calculate the determinant of the matrix
		\return the determinant
	*/	
	double getDeterminant() const;
	
	void singleValueDecomposition(Matrix &U, Vector &w, Matrix &v);
	
	void multiplyLineEscalar(int line, double escalar);
	void sumLines(int destLine,int line1, int line2);
	void sumLineVector(int destLine,int line, Vector v);
	void addColumnRight(Vector &b);
	void addColumnLeft(Vector &b);	
	Matrix getDiagonalMatrix() const;
	void changeLines(int line1, int line2);
	void readMatrixFromFile(const string fileName);
	void writeMatrixToFile(const string fileName) const;
	
	//! Get the norm  of the matrix
	/*!
	 	Calculate the norm of the matrix
		\return the norm 
	*/	
	double getNorm() const;
	
	//! Get the conditioning  of the matrix
	/*!
	 	Calculate the conditioning of the matrix
		\return the conditioning 
	*/	
	double getCond();	
	
	
protected:
	double **matrix;
	int numLines, numCols;
	void alocateMatrix(int numLines, int numCols);
	void copyMatrix(double **dest, double **src);
	void freeMatrix();
	Matrix reducedMatrix(Matrix &, int colunm);
	Vector retroSubstitution(Matrix A, Vector b);
};

#endif //_MATRIX_H_
