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

#include "RootFinder.h"

RootFinder::RootFinder(function1d f) {
	this->f = f;
}

RootFinder::~RootFinder()
{
}

void RootFinder::setFunction(function1d f) {

	this->f = f;

}

double RootFinder::findRootByBisectionMethod(double xR, double xL,double tolerance, int maxIterations) {
	
	double x,x0 = 0.0, error = 1;
	x = (xR + xL)/2.0;
	int count = 0;
		
	while ( (error > tolerance) &&(count < maxIterations)) {
		
		if(f(x)*f(xR) > 0) {			
				xR = x;			
		}
		else
			xL = x;
			
		x0 = x;
		x = (xR + xL)/2.0;
			
		error = fabs(x0-x);
		count++;
	}	
	numIterations = count;
	root = x0;
	return x0;

}

double RootFinder::findRootByBisectionMethod(function1d f, double xR, double xL,double tolerance, int maxIterations) {
	
	double x,x0 = 0.0, error = 1;
	x = (xR + xL)/2.0;
	int count = 0;
		
	while ( (error > tolerance) &&(count < maxIterations)) {
		
		if(f(x)*f(xR) > 0) {			
				xR = x;			
		}
		else
			xL = x;
			
		x0 = x;
		x = (xR + xL)/2.0;
			
		error = fabs(x0-x);
		count++;
	}	
	
	
	return x0;

}

double RootFinder::findRootBySecantMethod(double xn_1, double xn,  double tolerance, int maxIterations) {

	double d;
    
	for (int n = 1; n <= maxIterations; n++) {
        d = (xn - xn_1) / (f(xn) - f(xn_1)) * f(xn);
        if (fabs(d) < tolerance) {
        	numIterations = n;
        	return xn;
        }
        xn_1 = xn;
        xn = xn - d;
    }
    numIterations = maxIterations;
    return xn;
}

double RootFinder::findRootBySecantMethod(function1d f, double xn_1, double xn,  double tolerance, int maxIterations) {

	double d;
    
	for (int n = 1; n <= maxIterations; n++) {
        d = (xn - xn_1) / (f(xn) - f(xn_1)) * f(xn);
        if (fabs(d) < tolerance)
            return xn;
        xn_1 = xn;
        xn = xn - d;
    }
    return xn;

}

double RootFinder::findRootByNewtonMethod(function1d df, double x0, double tolerance, int maxIterations) {
	
	double x, error = 1;
	int count = 0;
		
	while ( (error > tolerance) &&(count < maxIterations)) {
		
		x = x0 - f(x0)/df(x0);
		error = fabs(x0-x);
		x0 = x;
		count++;
	}	
	
	numIterations = maxIterations;
	return x0;
	

}

double RootFinder::findRootByNewtonMethod(function1d f, function1d df, double x0, double tolerance, int maxIterations) {
	
	double x, error = 1;
	int count = 0;
		
	while ( (error > tolerance) &&(count < maxIterations)) {
		
		x = x0 - f(x0)/df(x0);
		error = fabs(x0-x);
		x0 = x;
		count++;
	}	
	
	
	return x0;
	

}

int RootFinder::getNumIterations() {

	return numIterations;

}

double RootFinder::verifyRoot() {
	return f(root);
}


