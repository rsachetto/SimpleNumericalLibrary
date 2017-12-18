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

#include "DirectSolver.h"

DirectSolver::DirectSolver(Matrix &A, Vector &b) {

	this->A = A;
	this->b = b;	

}

DirectSolver::~DirectSolver() {
	
}

void DirectSolver::setMatrixA(Matrix &A) {

	this->A = A;

}

void DirectSolver::setVectorB(Vector &b) {

	this->b = b;

}

Vector DirectSolver::gaussElimination() {

	if(A.isNull() || b.isNull()) {	
		cout << "Erro solve()!!"<< endl;	
		exit(0);
	}
	
	int n = A.getNumCols()-1;
	double m;
	Vector x(A.getNumCols());

	for(int k = 0; k <= n-1; k++) {
		for(int i = k+1; i <= n; i++) {
			m = A(i,k)/A(k,k);
			A(i,k) = 0;
			for(int j = k+1; j <= n; j++) {
				A(i,j) = A(i,j)-m*A(k,j);
			}			
			b[i] = b[i]-m*b[k];			
		}
	} 

	x = retroSubstitution(A,b);
	
	return x;
	
}

Vector DirectSolver::retroSubstitution(Matrix A, Vector b) {

	Vector x(b.getSize());
	int n = A.getNumLines() - 1;
	double s;
	x[n] = b[n]/A(n,n);
	
	for(int i = n-1; i >= 0; i--) {
		s = 0;
		for(int k = i+1; k <= n; k++) {
			s = s + A(i,k)*x[k];
			x[i] = (b[i]-s)/A(i,i);		
		}
	}	
	return x;
}

Vector DirectSolver::pivotGaussElimination() {
	
	int n = A.getNumCols()-1;
	double pv,aux,m;
	int r;

	for(int k = 0; k <= n-1; k++) {
	
		pv = fabs(A(k,k));
		r = k;
	
		for(int i = k+1; i <= n; i++) {
		
			if(fabs(A(i,k)) > pv) {
			
				pv = fabs(A(i,k));
				r = i;
			
			}		
		}
		if(pv == 0) {
		
			cout << "Singular Matrix" << endl;
			return Vector(n+1);
		
		}
		if(r != k) {
		
			for(int j = 0; j <= n; j++) {
			
				aux = A(k,j);
				A(k,j) = A(r,j);
				A(r,j) = aux;
			
			}
			aux = b[k];
			b[k] = b[r];
			b[r] = aux;
		
		}
		
		for(int i = k+1; i <= n; i++) {
			m = A(i,k)/A(k,k);
			A(i,k) = 0;
			for(int j = k+1; j <= n; j++) {
				A(i,j) = A(i,j)-m*A(k,j);
			}			
			b[i] = b[i]-m*b[k];			
		}	
	}	
	
	return retroSubstitution(A,b);

}

Vector DirectSolver::totalPivotGaussElimination() {
	
	int n = A.getNumCols()-1, aux2;
	double pv,aux,m;
	int r,s;
	Vector x(n+1);
	Vector p(n+1);
	
	for(int i = 0; i <= n; i++) {
		
		p[i] = i;
	
	}
	
	for(int k = 0; k <= n-1; k++) {
	
		pv = fabs(A(k,k));
		r = k;
		s = k;
	
		for(int i = k+1; i <= n; i++) {
		
			if(fabs(A(i,k)) > pv) {
			
				pv = fabs(A(i,k));
				r = i;
			
			}		
		}
		
		if(pv == 0) {
		
			cout << "Singular Matrix" << endl;
			return Vector(n+1);
		
		}
		if(r != k) {
		
			for(int j = 0; j <= n; j++) {
			
				aux = A(k,j);
				A(k,j) = A(r,j);
				A(r,j) = aux;
			
			}
			aux = b[k];
			b[k] = b[r];
			b[r] = aux;
		
		}
				
		for(int i = k; i <= n; i++) {
			for(int j = k; j <= n; j++) {
				if(fabs(A(i,j)) > pv) {				
					pv = fabs(A(i,s));
					s = j;
				
				}		
			}
		}
		
		if(s != k) {
		
			for(int j = 0; j <= n; j++) {
			
				aux = A(j,k);
				A(j,k) = A(j,s);
				A(j,s) = aux;
			
			}
			
			aux = p[k];
			p[k] = p[r];
			p[r] = aux;
		
		}
		
		for(int i = k+1; i <= n; i++) {
			m = A(i,k)/A(k,k);
			A(i,k) = 0;
			for(int j = k+1; j <= n; j++) {
				A(i,j) = A(i,j)-m*A(k,j);
			}			
			b[i] = b[i]-m*b[k];			
		}	
	}
		
	x = retroSubstitution(A,b);

	for(int i = 0; i <=n; i++) {
		
		aux2 = (int)p[i];
		x[i] = x[aux2];
	
	}
	
	return x;

	
}

Vector DirectSolver::LUGauss() {

	if(A.isNull() || b.isNull()) {	
		cout << "Erro solve()!!"<< endl;	
		exit(0);
	}
	
	int n = A.getNumCols()-1;
	double m,s,aux;
	Vector x(A.getNumCols()), y(A.getNumCols());

	for(int k = 0; k <= n-1; k++) {
		for(int i = k+1; i <= n; i++) {
			m = A(i,k)/A(k,k);
			A(i,k) = m;
			for(int j = k+1; j <= n; j++) {
				A(i,j) = A(i,j)-m*A(k,j);
			}		
		}
	}
	
	y[0] = b[0];
	
	for(int i = 1; i <= n; i++) {
		s = 0;
		for(int k = 0; k <= (i-1); k++) {
			if(i==k) {
				aux = 1;
			}
			else {
				aux = A(i,k);
			}
			s = s + aux * y[k];
			y[i] = b[i] - s;		
		}		
	
	}

	x = retroSubstitution(A,y);
	
	return x;

}

Vector DirectSolver::LUpivotGauss() {
	
	int n = A.getNumCols()-1;
	double pv,aux,m,s;
	int aux2;
	Vector y(A.getNumLines()), p(A.getNumLines()), c(A.getNumLines());
	int r;

	for(int i = 0; i <= n; i++) {
		
		p[i] = i;
	
	}

	for(int k = 0; k <= n-1; k++) {
	
		pv = fabs(A(k,k));
		r = k;
	
		for(int i = k+1; i <= n; i++) {
		
			if(fabs(A(i,k)) > pv) {
			
				pv = fabs(A(i,k));
				r = i;
			
			}		
		}
		if(pv == 0) {
		
			cout << "Singular Matrix" << endl;
			return Vector(n+1);
		
		}
		if(r != k) {
		
			for(int j = 0; j <= n; j++) {
			
				aux = A(k,j);
				A(k,j) = A(r,j);
				A(r,j) = aux;
			
			}
			aux = p[k];
			p[k] = p[r];
			p[r] = aux;
		
		}
		
		for(int i = k+1; i <= n; i++) {
			m = A(i,k)/A(k,k);
			A(i,k) = m;
			for(int j = k+1; j <= n; j++) {
				A(i,j) = A(i,j)-m*A(k,j);
			}						
		}	
	}	
	
	for(int i = 0; i <=n; i++) {
		
		aux2 = (int)p[i];
		c[i] = b[aux2];
	
	}
	
	y[0] = c[0];
	
	for(int i = 1; i <= n; i++) {
		s = 0;
		for(int k = 0; k <= (i-1); k++) {
			if(i==k) {
				aux = 1;
			}
			else {
				aux = A(i,k);
			}
			s = s + aux * y[k];
			y[i] = c[i] - s;		
		}		
	
	}
	
	return retroSubstitution(A,y);

}

Vector DirectSolver::refinedLUpivotGauss(int maxIterarions, double tolerance) {

	int n = A.getNumCols()-1;
	Matrix A_orig(A);
	double pv,aux,m,s;
	int aux2;
	Vector y(A.getNumLines()), p(A.getNumLines()), c(A.getNumLines());
	int r;
	Vector x(n+1);
	Vector e_0(n+1);
	Vector e_1(n+1);
	Vector res(n+1);
	
	for(int i = 0; i <= n; i++) {
		
		p[i] = i;
	
	}

	for(int k = 0; k <= n-1; k++) {
	
		pv = fabs(A(k,k));
		r = k;
	
		for(int i = k+1; i <= n; i++) {
		
			if(fabs(A(i,k)) > pv) {
			
				pv = fabs(A(i,k));
				r = i;
			
			}		
		}
		if(pv == 0) {
		
			cout << "Singular Matrix" << endl;
			return Vector(n+1);
		
		}
		if(r != k) {
		
			for(int j = 0; j <= n; j++) {
			
				aux = A(k,j);
				A(k,j) = A(r,j);
				A(r,j) = aux;
			
			}
			aux = p[k];
			p[k] = p[r];
			p[r] = aux;
		
		}
		
		for(int i = k+1; i <= n; i++) {
			m = A(i,k)/A(k,k);
			A(i,k) = m;
			for(int j = k+1; j <= n; j++) {
				A(i,j) = A(i,j)-m*A(k,j);
			}						
		}	
	}	
	
	for(int i = 0; i <=n; i++) {
		
		aux2 = (int)p[i];
		c[i] = b[aux2];
	
	}
	
	y[0] = c[0];
	
	for(int i = 1; i <= n; i++) {
		s = 0;
		for(int k = 0; k <= (i-1); k++) {
			if(i==k) {
				aux = 1;
			}
			else {
				aux = A(i,k);
			}
			s = s + aux * y[k];
			y[i] = c[i] - s;		
		}		
	
	}
	
	x = retroSubstitution(A,y);
	
	res = b - A_orig*x;
	
	for(int i = 0; i <=n; i++) {
		
		aux2 = (int)p[i];
		c[i] = res[aux2];
	
	}
	
	y[0] = c[0];
	
	for(int i = 1; i <= n; i++) {
		s = 0;
		for(int k = 0; k <= (i-1); k++) {
			if(i==k) {
				aux = 1;
			}
			else {
				aux = A(i,k);
			}
			s = s + aux * y[k];
			y[i] = c[i] - s;		
		}		
	
	}
	
	e_0 = retroSubstitution(A,y);
	
	x = x + e_0;
	
	for(int i = 0; i < maxIterarions; i++) {
	
		res = b-A_orig*x;
		
		for(int i = 0; i <=n; i++) {
			
			aux2 = (int)p[i];
			c[i] = res[aux2];
		
		}
		
		y[0] = c[0];
		
		for(int i = 1; i <= n; i++) {
			s = 0;
			for(int k = 0; k <= (i-1); k++) {
				if(i==k) {
					aux = 1;
				}
				else {
					aux = A(i,k);
				}
				s = s + aux * y[k];
				y[i] = c[i] - s;		
			}		
		
		}
		
		e_1 = retroSubstitution(A,y);
		x = x + e_1;
		
		if((e_1-e_0).norm2()/e_1.norm2() < tolerance) {
		
			return x;
		
		}
	
		for (int j = 0; j <= n; ++j) {
			e_0[j] = e_1[j];
		}	
	
	}
	
	return x;

}

Vector DirectSolver::doolittleFactorization() {

	int n = A.getNumCols()-1;
	double aux;
	
	for(int i = 1; i <= n; i++) {	
		A(i,0) = A(i,0)/A(0,0);	
	}
	
	for(int i = 1; i <= n; i++) {	
		for(int j = i; j <= n; j++) {
			aux = A(i,j);
			for(int k = 0; k <= i-1; k++) {			
				aux = aux - A(i,k)*A(k,j);			
			}
			A(i,j) = aux;
		}
		if(i != n) {
			for(int j = i+1; j <= n; j++) {
				aux = A(j,i);
				for(int k = 0; k <= i-1; k++) {			
					aux = aux - A(j,k)*A(k,i);			
				}
				A(j,i) = aux/A(i,i);
			}		
		}	
	}
	
	for(int i = 1; i <= n; i++) {
		for(int k = 0; k <= i-1; k++) {	
			if(i == k)
				aux = 1;
			else
				aux = A(i,k);	
			b[i] = b[i]-aux*b[k];		
		}
	}
	
	b = retroSubstitution(A,b);
	
	return b;

}

Vector DirectSolver::choleskyFactorization() {
	
	int n = A.getNumCols()-1;
	double aux;
	
	A(0,0) = sqrt(A(0,0));
	
	for(int j = 1; j <= n; j++) {
	
		A(0,j) = A(0,j)/A(0,0);	
		
	}
	
	for(int i = 1; i <= n; i++) {
	
		for(int j = i; j <= n; j++) {
		
			aux = A(i,j);
			for(int k = 0; k <= i-1; k++) {			
				aux = aux - A(k,i)*A(k,j);			
			}
			
			if(j == i) {
			
				A(i,i) = sqrt(aux);
			
			}
			else {
			
				A(i,j) = aux/A(i,i);
				
			}
		
		}
	
	}
	
	b[0] = b[0]/A(0,0);
	
	for(int i = 1; i <= n; i++) {
		aux = b[i];
		for(int k = 0; k <= i-1; k++) {	
			aux = aux-A(k,i)*b[k];		
		}
		b[i] = aux/A(i,i);
	}
	
	b[n] = b[n]/A(n,n);
	
	for(int i = n-1; i >= 0; i--) {
	
		aux = b[i];
		for(int k = i+1; k <= n; k++) {	
			aux = aux-A(i,k)*b[k];		
		}
		b[i] = aux/A(i,i);
	
	}
	
	return b;
	
}

void DirectSolver::printSystem() {

	if(A.isNull() || b.isNull()) {	
		cout << "Error printSystem()!! Verify o matrix A and vector b!!" << endl;	
		return;
	}

	double num;
	char variable = 'a';
	
	for(int i = 0; i < A.getNumLines(); i++) {
		for(int j = 0; j < A.getNumCols(); j++) {		
			num = A(i,j);
			
			if(j == 0) {			
				cout << num << variable;	
				variable++;		
			}
			else {
				if(num > 0) {			
					cout << "+" << num << variable;
					variable++;			
				}
				else if (num < 0){
					cout << num << variable;			
					variable++;
				}
				else
					variable++;
			}	
		}		
		variable = 'a';	
		cout << " = " << b[i];		
		cout << endl;		
	}

}

