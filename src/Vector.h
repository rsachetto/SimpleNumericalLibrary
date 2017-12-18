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

#ifndef VECTOR_H_
#define VECTOR_H_

#include <string>
#include <cmath>
#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<iomanip>
#include <stdexcept>

using namespace std;


class VectorException : public runtime_error
{
    public:
        VectorException(const string &msg) : runtime_error(msg){}
};

//! An n Vector.
/*!
A class that represents a n Vector of doubles
*/
class Vector
{
public:
	//! Constructor.
	/*!
		Construct a empty vector;
	*/
	Vector(){vector = NULL;};
	
	//! Constructor.
	/*!
		 Construct a vector with n elements
		\param size number of elements of the vector
	*/
	Vector(int size);
	
	//! Constructor.
	/*!
		 Construct a copy of the Vector v
		\param v the Vector to be copied
	*/
	Vector(const Vector &v);
	
	//! Constructor.
	/*!
		 Construct a vector with the values in vec
		\param vec source array
		\param size the size of source array
	*/
	Vector(double *vec, int size);
	
	//! Destructor.
	/*!
		Destroy the objects;
	*/
	virtual ~Vector();	
	
	//! Get Method
	/*!
	 	Get the size of the Vector
		\return the size of the Vector
	*/
	int getSize() const;
	
	//! Get Method
	/*!
	 	Get a element of the vector 
	 	\param index the position of the element
		\return the element located at position index
	*/
	double getI( int index ) const;
	
	//! Set Method
	/*!
	 	Set a element of the vector 
		\param index the position of the element
		\param value the value to set
	*/
	void setI( int index, double value );
	
	//! Print the vector
	/*!
	 	Print the vector in a linear form
	*/	
	void printVector() const;
	
	//! Overloaded operator =
	/*!
	 	Attributes a vector to another
	 	\param rhs the vector to be attributed
	*/
	const Vector& operator = (const Vector & rhs);
	
	//! Overloaded operator +
	/*!
	 	Sums two vectors
	 	\param rhs the right-hand side of the operation
	 	\return the result of operation
	*/
	Vector operator + (const Vector & rhs) const;
	
	//! Overloaded operator -
	/*!
	 	Substracts two vectors
	 	\param rhs the right-hand side of the operation
	 	\return the result of operation
	*/
	Vector operator - (const Vector & rhs) const;
	
	//! Overloaded operator -
	/*!
	 	Invert the signal of a vector
	 	\return the result of operation
	*/
	Vector operator - () const;
	
	//! Overloaded operator *
	/*!
	 	Multiplies two vectors
	 	\param rhs the right-hand side of the operation
	 	\return the result of operation
	*/
	Vector operator * (const Vector & rhs) const;
	
	//! Overloaded operator *
	/*!
	 	Multiplies a vector by a scalar
	 	\param escalar the right-hand side of the operation
	 	\return the result of operation
	*/
	Vector operator * (double escalar);
	
		//! Overloaded operator /
	/*!
	 	Divides a vector by a scalar
	 	\param escalar the right-hand side of the operation
	 	\return the result of operation
	*/
	Vector operator / (double escalar);
	
	//! Overloaded operator []
	/*!
	 	Set a element of the matrix 
		\param i the position to put the element
		
		ex: vector[0] = 5.0;
	*/
	double &operator [] (int i);
	
	//! Overloaded operator []
	/*!
	 	Get a element of the vector 
	 	\param i the position to put the element
		\return the element located at position i
		
		ex: double value = vector[0];
	*/
	double operator [] (int i) const;
	
	//! Calculates inner product of two vectors
	/*!
	 	Calculates inner product of two vectors
	 	\param a the first vector
	 	\param b the second vector
		\return the inner product
		
		ex: double inner = Vector::innerProduct(v1,v2);
	*/
	static double innerProduct(Vector a, Vector b);
	
	//! Get the norm2  of the vector
	/*!
	 	Calculate the norm2 of the vector
		\return the norm 
	*/	
	double norm2() const;
	
	//! Verify the allocation of the vector
	/*!
	 	Verify if the vector was allocated	 	
		\return true if the vector was allocates, false otherwise
	*/
	bool isNull();
	
	double multipilyComponents();
	void readVectorFromFile(const string);
	void writeVectorToFile(const string) const;
	friend ostream &operator<<( ostream &, const Vector & );

private:
	double *vector;
	int size;
	void alocateVector();
	void freeVector();
	void copyVector(double *dest, double *src);
	
};

#endif /*VECTOR_H_*/
