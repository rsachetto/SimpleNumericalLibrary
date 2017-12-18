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

#include "TridiagonalSolver.h"

TridiagonalSolver::TridiagonalSolver(TridiagonalMatrix &A, Vector &b) {
	
	this->A = A;
	this->b = b;
}

TridiagonalSolver::~TridiagonalSolver()
{
}

void TridiagonalSolver::setMatrixA(TridiagonalMatrix &A) {

	this->A = A;

}

void TridiagonalSolver::setVectorB(Vector &b) {

	this->b = b;

}

Vector TridiagonalSolver::gaussElimination() {

	int n = b.getSize()-1;
	Vector lower(n+1);
	Vector principal(n+1);
	Vector upper(n+1);
	Vector x(n+1);
	
	for (int i = 0; i <= n; ++i) {
		upper[i] = A.getUpperDiagonalValue(i);
		principal[i] = A.getPrincipalDiagonalValue(i);
		lower[i] = A.getLowerDiagonalValue(i);
	}
	
	double m;
	
	for(int k = 1; k <= n; k++) {
	
		m = lower[k]/principal[k-1];
		principal[k] = principal[k] - m*upper[k-1];
		b[k] = b[k] - m*b[k-1];
	
	}
	
	x[n] = b[n]/principal[n];
	
	for(int k = n-1; k >=0; k--) {
	
		x[k] = (b[k] - upper[k]*x[k+1])/principal[k];
	
	}
	
	return x;
}
