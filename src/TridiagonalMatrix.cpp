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

#include "TridiagonalMatrix.h"

TridiagonalMatrix::TridiagonalMatrix(int diagonalSize) {	
	this->diagonalSize = diagonalSize;
	alocateMatrix(diagonalSize);
	lowerDiagonal[0] = 0;
	upperDiagonal[diagonalSize-1] = 0;
}

TridiagonalMatrix::TridiagonalMatrix() {	
	lowerDiagonal = NULL;
	principalDiagonal = NULL;
	upperDiagonal = NULL;	
}

TridiagonalMatrix::~TridiagonalMatrix(){
	
	free(lowerDiagonal);
	free(principalDiagonal);
	free(upperDiagonal);
	
}

TridiagonalMatrix::TridiagonalMatrix(const TridiagonalMatrix &m) {
	
	diagonalSize = m.diagonalSize;
	alocateMatrix(diagonalSize);
	copyMatrix(lowerDiagonal,principalDiagonal,upperDiagonal, m.lowerDiagonal, m.principalDiagonal, m.upperDiagonal);

}

void TridiagonalMatrix::alocateMatrix(int diagonalSize) {

	lowerDiagonal = (double*)calloc(diagonalSize,sizeof(double));
	principalDiagonal = (double*)calloc(diagonalSize,sizeof(double));
	upperDiagonal = (double*)calloc(diagonalSize,sizeof(double));

	if(!lowerDiagonal || !principalDiagonal || !upperDiagonal) {
		cerr << "Error on memory allocation";
		exit(0);	
	}	
}

void TridiagonalMatrix::copyMatrix(double *destLower, double *destPrinc, double *destUpper, double *srcLower, double *srcPrinc, double *srcUpper) {

	for(int i = 0; i < diagonalSize; i++) {	
		destLower[i] = srcLower[i];
		destPrinc[i] = srcPrinc[i];
		destUpper[i] = srcUpper[i];		
	}

}

void TridiagonalMatrix::freeDiagonal(double *diagonal) {

	free(diagonal);

}

const TridiagonalMatrix& TridiagonalMatrix::operator =(const TridiagonalMatrix &rhs) {

	if(this == &rhs)
		return *this;
	
	if(lowerDiagonal != NULL)		
		freeDiagonal(lowerDiagonal);
		
	if(principalDiagonal != NULL)		
		freeDiagonal(principalDiagonal);	
	
	if(upperDiagonal != NULL)		
		freeDiagonal(upperDiagonal);
		
	this->diagonalSize = rhs.diagonalSize;
	alocateMatrix(rhs.diagonalSize);
	copyMatrix(lowerDiagonal, principalDiagonal, upperDiagonal, rhs.lowerDiagonal, rhs.principalDiagonal, rhs.upperDiagonal);
	
	return *this;

}

double TridiagonalMatrix::operator () (int i, int j) const {
	
	if((i < diagonalSize) && (j < diagonalSize)) {
		if(fabs(i-j)>1)
			return 0.0;
		else if((i == j) && (i < diagonalSize))
			return principalDiagonal[i];
		else if ((i+1 == j) && (i < diagonalSize-1)) 
			return upperDiagonal[i];
		else
			return lowerDiagonal[i]; 
	}
	else {
		cout << "Matrix Index out of bounds: "<< "i = " << i << ", j = " << j; 
		exit(0);
	}
}

double TridiagonalMatrix::getIJ(int i, int j) {
	
	if((i < diagonalSize) && (j < diagonalSize)) {
		if(fabs(i-j)>1)
			return 0;
		if((i == j) && (i < diagonalSize))
			return principalDiagonal[i];
		else if ((i+1 == j) && (i < diagonalSize-1)) 
			return upperDiagonal[i];
		else 
			return lowerDiagonal[i]; 
	}
	else {
		cout << "Matrix Index out of bounds: "<< "i = " << i << ", j = " << j; 
		exit(0);
	}
}

void TridiagonalMatrix::setIJ(int i, int j, double value) {	
	
	if((i < diagonalSize) && (j < diagonalSize)) {
		if((i == j) && (i < diagonalSize))
			principalDiagonal[i] = value;
		else if ((i+1 == j) && (i < diagonalSize-1)) 
			upperDiagonal[i] = value;
		else if ((j == i-1) && (i > 0))
			lowerDiagonal[i] = value; 
	}
	else {
		cout << "Matrix Index out of bounds: "<< "i = " << i << ", j = " << j; 
		exit(0);
	}
}

void TridiagonalMatrix::printMatrix() {

	for(int i = 0; i < diagonalSize; i++) {
		for(int j = 0; j < diagonalSize; j++) {			
			cout << this->getIJ(i,j) << "\t";	
		}
		cout << endl;
	}

}

double TridiagonalMatrix::getUpperDiagonalValue(int index) {	
	return upperDiagonal[index];
}

double TridiagonalMatrix::getPrincipalDiagonalValue(int index) {
	return principalDiagonal[index];
}

double TridiagonalMatrix::getLowerDiagonalValue(int index) {
	return lowerDiagonal[index];
}
