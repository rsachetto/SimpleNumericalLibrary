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

#include "NonLinearSolver.h"

NonLinearSolver::NonLinearSolver(function **jacobi, function *f, Vector x0) {
	
	this->n = x0.getSize();
	this->x = x0;
	this->jacobi = jacobi;
	this->f = f;
	//alocateJacobi();
	//alocateF();
	//copyJacoby(this->jacobi, jacobi);
	//copyF(this->f,f); 
}

NonLinearSolver::~NonLinearSolver() {
	
	if(jacobi) {
		for(int i = 0; i < n; i++) {
			if(jacobi[i])
				free(jacobi[i]);
		}
		if(jacobi)
			free(jacobi);
	}
	
	free(f);	
	
}

function ** NonLinearSolver::alocateJacobi(int n) {

	function ** jacobi;
	jacobi = (function**)malloc(n*sizeof(function));
	
	if(!jacobi){
		cerr << "Error on memory allocation";
		exit(0);
	}
	
	for(int i = 0; i < n; i++) {
		jacobi[i] = (function*)malloc(n*sizeof(function));
		if(!jacobi[i]){
			 cerr << "Error on memory allocation";
			 exit(0);
		}
	}
	return jacobi;
}

function * NonLinearSolver::alocateF(int n) {

	function *f;
	f = (function*)malloc(n*sizeof(function));
	if(!f) {
		 cerr << "Error on memory allocation";
		 exit(0);
	}
	return f;
}

void NonLinearSolver::copyJacoby(function **dest, function **src) {
	
	for(int i = 0; i < n; i++) {	
		for(int j = 0; j < n; j++) {
			dest[i][j] = src[i][j];
		}		
	}
}

void NonLinearSolver::copyF(function *dest, function *src) {
	
	for(int i = 0; i < n; i++) {	
		dest[i] = src[i];				
	}
}

Vector NonLinearSolver::solve(int maxIter, double tolerance) {
	
	Matrix J(n,n);
	Vector F(n);
	Vector s(n);
	Vector x1(n);
	DirectSolver *ds;
	int k = 0;	
	
		
	while(k < maxIter) {
	
		for(int i = 0; i < n; i++) {
			
			F[i] = -f[i](x); 
			for(int j = 0; j < n; j++) {
				
				J(i,j) = jacobi[i][j](x);				
				
			}				
			
		}
	
		ds = new DirectSolver(J,F);
		
		s = ds->pivotGaussElimination();
		x1 = x + s;	 
		
		if(stopCriteria(x, x1, tolerance)) {
			nIter = k;
			return x1;
		}
		
		x = x1;
		
		k++;
	
	}
	
	nIter = k;
	
	return x;
	
}

int NonLinearSolver::getNIterations() {	
	return nIter;
}

bool NonLinearSolver::stopCriteria(Vector x_0, Vector x_1, double tolerance) {

	return ((x_1-x_0).norm2() < tolerance);

}
