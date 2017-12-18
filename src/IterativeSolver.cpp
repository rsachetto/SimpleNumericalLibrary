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

#include "IterativeSolver.h"

IterativeSolver::IterativeSolver(Matrix &A, Vector &b){
	
	this->A = A;
	this->b = b;	
}

IterativeSolver::~IterativeSolver() {
}

Vector IterativeSolver::jacob(Vector x_0, int maxIteration, double tolerance) {

	int k = 1;
	int n = A.getNumCols()-1;
	double aux;
	Vector x_1(A.getNumCols());
	
	while(k <= maxIteration) {
	
		for(int i = 0; i <= n; i++) {
					
			aux = b[i];
			
			for(int j=0; j <= n; j++) {
			
				if(j != i) {				
					aux = aux - A(i,j)*x_0[j];				
				}				
			}
						
			x_1[i] = aux/A(i,i);
		}
				
		if(stopCriteriaJacob(x_0,x_1) <= tolerance) {
		
			return x_1;
		
		}
		else {
		
			for(int j = 0; j <= n; j++) {			
				x_0[j] = x_1[j];
			
			}		
		}
		k++;
	}
	
	return x_1;
	
}

Vector IterativeSolver::jacob(int maxIteration, double tolerance) {

	int k = 1;
	int n = A.getNumCols()-1;
	double aux;
	Vector x_1(A.getNumCols());
	Vector x_0(A.getNumCols());
	
	while(k <= maxIteration) {
	
		for(int i = 0; i <= n; i++) {
					
			aux = b[i];
			
			for(int j=0; j <= n; j++) {
			
				if(j != i) {				
					aux = aux - A(i,j)*x_0[j];				
				}				
			}
						
			x_1[i] = aux/A(i,i);
		}
				
		if(stopCriteriaJacob(x_0,x_1) <= tolerance) {
		
			return x_1;
		
		}
		else {
		
			for(int j = 0; j <= n; j++) {			
				x_0[j] = x_1[j];
			
			}		
		}
		k++;
	}
	
	return x_1;
	
}

Vector IterativeSolver::gaussSeidel(Vector x_0, int maxIteration, double tolerance) {

	int k = 1;
	int n = A.getNumCols()-1;
	double aux;
	Vector x_1(A.getNumCols());
	
	while(k <= maxIteration) {
	
		for(int i = 0; i <= n; i++) {
					
			aux = b[i];
			
			for(int j=0; j <= n; j++) {
			
				if(j > i) {				
					aux = aux - A(i,j)*x_0[j];				
				}
				else if(j < i) {					
					aux = aux - A(i,j)*x_1[j];
				}				
			}
						
			x_1[i] = aux/A(i,i);
		}
				
		if(stopCriteriaJacob(x_0,x_1) <= tolerance) {
		
			return x_1;
		
		}
		else {
		
			for(int j = 0; j <= n; j++) {			
				x_0[j] = x_1[j];			
			}		
		}
		k++;
	}
	
	return x_1;
	
}

Vector IterativeSolver::gaussSeidel(int maxIteration, double tolerance) {

	int k = 1;
	int n = A.getNumCols()-1;
	double aux;
	Vector x_1(A.getNumCols());
	Vector x_0(A.getNumCols());
	
	while(k <= maxIteration) {
	
		for(int i = 0; i <= n; i++) {
					
			aux = b[i];
			
			for(int j=0; j <= n; j++) {
			
				if(j > i) {				
					aux = aux - A(i,j)*x_0[j];				
				}
				else if(j < i) {					
					aux = aux - A(i,j)*x_1[j];
				}				
			}
						
			x_1[i] = aux/A(i,i);
		}
				
		if(stopCriteriaJacob(x_0,x_1) <= tolerance) {
		
			return x_1;
		
		}
		else {
		
			for(int j = 0; j <= n; j++) {			
				x_0[j] = x_1[j];			
			}		
		}
		k++;
	}
	
	return x_1;
	
}

Vector IterativeSolver::SOR(Vector x_0, int maxIteration, double tolerance, double omega) {

	int k = 1;
	int n = A.getNumCols()-1;
	double aux;
	Vector x_1(A.getNumCols());
	
	while(k <= maxIteration) {
	
		for(int i = 0; i <= n; i++) {
					
			aux = b[i];
			
			for(int j=0; j <= n; j++) {
			
				if(j > i) {				
					aux = aux - A(i,j)*x_0[j];				
				}
				else if(j < i) {					
					aux = aux - A(i,j)*x_1[j];
				}				
			}
						
			x_1[i] =(1-omega)*x_0[i] + omega*(aux/A(i,i));
		}
				
		if(stopCriteriaJacob(x_0,x_1) <= tolerance) {
		
			return x_1;
		
		}
		else {
		
			for(int j = 0; j <= n; j++) {			
				x_0[j] = x_1[j];			
			}		
		}
		k++;
	}
	
	return x_1;
	
}

Vector IterativeSolver::SOR(int maxIteration, double tolerance, double omega) {

	int k = 1;
	int n = A.getNumCols()-1;
	double aux;
	Vector x_1(A.getNumCols());
	Vector x_0(A.getNumCols());
	
	while(k <= maxIteration) {
	
		for(int i = 0; i <= n; i++) {
					
			aux = b[i];
			
			for(int j=0; j <= n; j++) {
			
				if(j > i) {				
					aux = aux - A(i,j)*x_0[j];				
				}
				else if(j < i) {					
					aux = aux - A(i,j)*x_1[j];
				}				
			}
						
			x_1[i] =(1-omega)*x_0[i] + omega*(aux/A(i,i));
		}
				
		if(stopCriteriaJacob(x_0,x_1) <= tolerance) {
		
			return x_1;
		
		}
		else {
		
			for(int j = 0; j <= n; j++) {			
				x_0[j] = x_1[j];			
			}		
		}
		k++;
	}
	
	return x_1;
	
}

Vector IterativeSolver::gradientMethod(Vector v_0, int maxIteration, double tolerance) {

int n = b.getSize();
	Vector r_0(n);
	Vector r_1(n);
	Vector v_1(n);
	
	int k = 1;
	double t_min;
	
	r_0 = A*v_0 - b;
	
	while(k <= maxIteration) {
	
		t_min = Vector::innerProduct(r_0,r_0)/Vector::innerProduct(A*r_0,r_0);
		
		v_1 = v_0 - (r_0*t_min);
		
		r_1 = r_0 - (A*r_0)*t_min;
		
		if(stopCriteriaGradient(v_0,v_1,r_1,tolerance)) {
			return v_1;
		}
		for (int i = 0;  i < n; i++) {
			r_0[i] = r_1[i];
			v_0[i] = v_1[i];
		}
		
		k++;	
	}

	return v_1;

}	

Vector IterativeSolver::gradientMethod(int maxIteration, double tolerance) {


	int n = b.getSize();
	Vector r_0(n);
	Vector r_1(n);
	Vector v_0(n);
	Vector v_1(n);
	
	int k = 1;
	double t_min;
	
	r_0 = A*v_0 - b;
	
	while(k <= maxIteration) {
	
		t_min = Vector::innerProduct(r_0,r_0)/Vector::innerProduct(A*r_0,r_0);
		
		v_1 = v_0 - (r_0*t_min);
		
		r_1 = r_0 - (A*r_0)*t_min;
		
		if(stopCriteriaGradient(v_0,v_1,r_1,tolerance)) {
			return v_1;
		}
		for (int i = 0;  i < n; i++) {
			r_0[i] = r_1[i];
			v_0[i] = v_1[i];
		}
		
		k++;	
	}

	return v_1;

}

Vector IterativeSolver::conjugatedGradientMethod(int maxIteration, double tolerance) {

	int n = b.getSize();
	Vector r_0(b.getSize());
	Vector r_1(b.getSize());
	Vector v_0(b.getSize());
	Vector v_1(b.getSize());
	Vector p_0(b.getSize());
	Vector p_1(b.getSize());
	
	int k = 2;
	double q, alpha;
	
	r_0 = A*v_0 - b;
	
	p_0 = -r_0;
	
	q = Vector::innerProduct(r_0,r_0)/Vector::innerProduct(A*r_0,r_0);
	
	v_1 = v_0 + p_0*q;
	r_1 = r_0 + (A*p_0)*q;
	
	if(stopCriteriaGradient(v_0,v_1,r_1,tolerance)) {
			return v_1;
	}
	
	while(k <= maxIteration) {
	
	
		alpha = Vector::innerProduct(r_1,r_1)/Vector::innerProduct(r_0,r_0);
		
		p_1 = -r_1 + p_0*alpha;
		
		
		q = Vector::innerProduct(r_1,r_1)/Vector::innerProduct(A*p_1,p_1);
	
		
		v_1 = v_0 + (p_1*q);
		r_1 = r_0 + (A*p_1)*q;
		
		
		if(stopCriteriaGradient(v_0,v_1,r_1,tolerance)) {
			return v_1;
		}
		for (int i = 0;  i < n; i++) {
			r_0[i] = r_1[i];
			v_0[i] = v_1[i];
			p_0[i] = p_1[i];
		}
		k++;	
	}

	return v_1;

}

Vector IterativeSolver::conjugatedGradientMethod(Vector v_0, int maxIteration, double tolerance) {

	int n = b.getSize();
	Vector r_0(b.getSize());
	Vector r_1(b.getSize());
	Vector v_1(b.getSize());
	Vector p_0(b.getSize());
	Vector p_1(b.getSize());
	
	int k = 2;
	double q, alpha;
	
	r_0 = A*v_0 - b;
	
	p_0 = -r_0;
	
	q = Vector::innerProduct(r_0,r_0)/Vector::innerProduct(A*r_0,r_0);
	
	v_1 = v_0 + p_0*q;
	r_1 = r_0 + (A*p_0)*q;
	
	if(stopCriteriaGradient(v_0,v_1,r_1,tolerance)) {
			return v_1;
	}
	
	while(k <= maxIteration) {
	
	
		alpha = Vector::innerProduct(r_1,r_1)/Vector::innerProduct(r_0,r_0);
		
		p_1 = -r_1 + p_0*alpha;
		
		
		q = Vector::innerProduct(r_1,r_1)/Vector::innerProduct(A*p_1,p_1);
	
		
		v_1 = v_0 + (p_1*q);
		r_1 = r_0 + (A*p_1)*q;
		
		
		if(stopCriteriaGradient(v_0,v_1,r_1,tolerance)) {
			return v_1;
		}
		for (int i = 0;  i < n; i++) {
			r_0[i] = r_1[i];
			v_0[i] = v_1[i];
			p_0[i] = p_1[i];
		}
		k++;	
	}

	return v_1;

}

double IterativeSolver::stopCriteriaJacob(Vector x_0, Vector x_1) {
	
	double max1, max2, aux;
	int n = x_0.getSize()-1;
	
	max1 = fabs(x_1[0] - x_0[0]);
	max2 = fabs(x_1[0]);
	
	for(int i = 1; i <= n; i++) {
	
		aux = fabs(x_1[i] - x_0[i]);
		
		if(aux > max1) {		
			max1 = aux;		
		}
		
		aux = fabs(x_1[i]);
		if(aux > max2) {		
			max2 = aux;		
		}
	}
		
	return max1/max2;
	
}

bool IterativeSolver::stopCriteriaGradient(Vector v_0, Vector v_1, Vector r_1, double tolerance) {

	return ((r_1.norm2() < tolerance) || (((v_1-v_0).norm2()/v_1.norm2()) < tolerance));

}

void IterativeSolver::printSystem() {

	if(A.isNull() || b.isNull()) {	
		cout << "Erro imprimeSistema()!! Verique o matriz A e o vetor b do sistema!!" << endl;	
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

Vector IterativeSolver::conjugatedGradientMethod2(Vector x0, int maxIteration, double tolerance) {
	
	int dim = x0.getSize();
	Vector x(dim),r(dim),v(dim),z(dim);
	double c,t,d;
	x = x0;
	r = b - A*x;
	v = r;
	c = Vector::innerProduct(r,r);
	
	for(int i = 0;i < dim; i++) {
		if(sqrt(Vector::innerProduct(v,v)) < tolerance) {
			cerr << "Error in ConjugateGradient: execution ";
			cerr << "of function terminated" << endl;
			break;
		}
	
		z = A*v;
		t = c/Vector::innerProduct(v,z);
		x = x + v*t;
		r = r - z*t;
		d = Vector::innerProduct(r,r);
		
		if(sqrt(d) < tolerance)
			break;
			
		v = r + v*(d/c);
		c = d;
	}
	return x;
}
