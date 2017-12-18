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

#include "GaussQuadrature.h"

GaussQuadrature::GaussQuadrature(function1d f)
{
	this->f = f;
}

GaussQuadrature::~GaussQuadrature()
{
}

bool GaussQuadrature::getQuadraturePoints(int n) {
	
	if(n == 2) {
		points = (double *)malloc(n*sizeof(double));
		weights = (double *)malloc(n*sizeof(double));		
		points[0] = -0.57735026918962576451;
		points[1] = 0.57735026918962576451;		
		weights[0] = 1.0;
		weights[1] = 1.0;
		return true;		
	}	
	
	else if(n == 3) {
		points = (double *)malloc(n*sizeof(double));
		weights = (double *)malloc(n*sizeof(double));	
			
		points[0] = 0.0;
		points[1] = -0.77459667;
		points[2] = 0.77459667;	
		weights[0] = 0.88888889;
		weights[1] = 0.55555555;
		weights[2] = 0.55555555;
		return true;		
	}
	
	else if(n == 4) {
			
		points = (double *)malloc(n*sizeof(double));
		weights = (double *)malloc(n*sizeof(double));
				
		points[0] = -0.33998104;
		points[1] = 0.33998104; 
		points[2] = -0.86113631; 
		points[3] = 0.86113631;
			
		weights[0] = 0.65214515;
		weights[1] = 0.65214515;
		weights[2] = 0.34785485;
		weights[3] = 0.34785485;
		return true;		
	}
	
	else if(n == 5) {
			
		points = (double *)malloc(n*sizeof(double));
		weights = (double *)malloc(n*sizeof(double));
				
		points[0] = 0.0;
		points[1] = 0.53846931; 
		points[2] = -0.53846931; 
		points[3] = 0.90617985;
		points[4] = -0.90617985;
		
		weights[0] = 0.56888889;
		weights[1] = 0.47862867;
		weights[2] = 0.47862867;
		weights[3] = 0.23692689;
		weights[4] = 0.23692689;
		return true;		
	}
	
	else if(n == 6) {
			
		points = (double *)malloc(n*sizeof(double));
		weights = (double *)malloc(n*sizeof(double));
				
		points[0] = -0.23861918;
		points[1] = 0.23861918; 
		points[2] = -0.66120939; 
		points[3] = 0.66120939;
		points[4] = -0.93246951;
		points[5] = 0.93246951;
		
		weights[0] = 0.46791393;
		weights[1] = 0.46791393;
		weights[2] = 0.36076157;
		weights[3] = 0.36076157;
		weights[4] = 0.17132449;
		weights[5] = 0.17132449;
		return true;		
	}
	
	else if(n == 7) {
			
		points = (double *)malloc(n*sizeof(double));
		weights = (double *)malloc(n*sizeof(double));
				
		points[0] = 0.0;
		points[1] = 0.40584515; 
		points[2] = -0.40584515; 
		points[3] = 0.74153119;
		points[4] = -0.74153119;
		points[5] = 0.94910791;
		points[6] = -0.94910791;
		
		weights[0] = 0.41795918;
		weights[1] = 0.38183005;
		weights[2] = 0.38183005;
		weights[3] = 0.27970539;
		weights[4] = 0.27970539;
		weights[5] = 0.12948497;
		weights[6] = 0.12948497;
		return true;		
	}
	
	else if(n == 8) {
			
		points = (double *)malloc(n*sizeof(double));
		weights = (double *)malloc(n*sizeof(double));
				
		points[0] = 0.18343464;
		points[1] = -0.18343464; 
		points[2] = 0.52553241; 
		points[3] = -0.52553241;
		points[4] = 0.79666648;
		points[5] = -0.79666648;
		points[6] = 0.96028986;
		points[7] = -0.96028986; 	
		
		weights[0] = 0.36268378;
		weights[1] = 0.36268378;
		weights[2] = 0.31370665;
		weights[3] = 0.31370665;
		weights[4] = 0.22238103;
		weights[5] = 0.22238103;
		weights[6] = 0.10122854;
		weights[7] = 0.10122854;
		return true;		
	}
	else
		return false;
	
}	

double GaussQuadrature::integrate(double a, double b, int n) {
	
	double dx = (b-a)/2.0;
	double integral = 0.0;
	double bpa2 = (b+a)/2;
	
	if(!getQuadraturePoints(n)) {
		cout << "Error getting Gauss' points" << endl << "Verify if n <= 8"<< endl;
		return 0.0;
	}
	
	for(int i = 0; i < n; i++) {
	
		integral += f(dx*points[i] + bpa2)*weights[i];
	
	}
	
	return integral*dx;
}
