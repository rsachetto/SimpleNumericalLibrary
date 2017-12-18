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

#ifndef NONLINEARSOLVER_H_
#define NONLINEARSOLVER_H_
#include "Vector.h"
#include "Matrix.h"
#include "DirectSolver.h"
#include "Function.h"
#include<cstdio>
#include<iostream>
#include<cstdlib>

class NonLinearSolver
{
public:
	NonLinearSolver(function **jacobi, function *f, Vector x0);
	virtual ~NonLinearSolver();	
	Vector solve(int maxIter, double tolerance);
	int getNIterations();
	static function ** alocateJacobi(int n);
	static function * alocateF(int n);
	
private:
	
	int n;
	int nIter;
	function **jacobi;
	function *f;
	Vector x;	
	void copyJacoby(function **dest, function **src);
	void copyF(function *dest, function *src);
	bool stopCriteria(Vector x_0, Vector x_1, double tolerance);
	
};

#endif /*NONLINEARSOLVER_H_*/
