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

#include "MatrixFileHandler.h"

MatrixFileHandler::~MatrixFileHandler()
{
	fclose(f);
}

bool MatrixFileHandler::readMatrixFile(string fileName, double **destMatrix, int numLines, int numCols){

	f = fopen(fileName.c_str(),"r");
	
	if(destMatrix == NULL) {
		cout << "error readMatrixFile, alocate space for destMatrix" << endl;
		return false;
	}
	
	if(f == NULL) {
		cout << "error readMatrixFile, could't open file " << fileName << endl;
		return false;
	}
		

	double element;
	for(int i = 0; i < numLines; ++i) {
		for(int j = 0; j < numCols; ++j) {
		
			fscanf(f,"%lf", &element);
			destMatrix[i][j] = element;
			
		}		
	}
	return true;
}

bool MatrixFileHandler::writeMatrixFile(string fileName, double **srcMatrix, int numLines, int numCols){

	f = fopen(fileName.c_str(),"w");
	
	if(srcMatrix == NULL) {
		cout << "error readMatrixFile, alocate space for srcMatrix" << endl;
		return false;
	}
	
	if(f == NULL) {
		cout << "error readMatrixFile, could't open file " << fileName << endl;
		return false;
	}		

	double element;
	for(int i = 0; i < numLines; ++i) {
		for(int j = 0; j < numCols; ++j) {			
			element = srcMatrix[i][j];		
			fprintf(f,"%lf ", element);			
		}
			fprintf(f,"\n");				
	}
	return true;
}
