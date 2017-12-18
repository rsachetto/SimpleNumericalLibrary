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

#include "VectorFileHandler.h"

VectorFileHandler::~VectorFileHandler()
{
	fclose(f);
}


bool VectorFileHandler::readVectorFile(string fileName, double *destVector, int size){

	f = fopen(fileName.c_str(),"r");
	
	if(destVector == NULL) {
		cout << "error readVectorFile, alocate space for destVector" << endl;
		return false;
	}
	
	if(f == NULL) {
		cout << "error readVectorFile, could't open file " << fileName << endl;
		return false;
	}		

	double element;
	
	for(int i = 0; i < size; ++i) {
		fscanf(f,"%lf", &element);
		destVector[i] = element;				
	}
	return true;
}

bool VectorFileHandler::writeVectorFile(string fileName, double *srcVector, int size){

	f = fopen(fileName.c_str(),"w");
	
	if(srcVector == NULL) {
		cout << "error readVectorFile, alocate space for destVector" << endl;
		return false;
	}
	
	if(f == NULL) {
		cout << "error readVectorFile, could't open file " << fileName << endl;
		return false;
	}		

	double element;
	for(int i = 0; i < size; ++i) {		
			element = srcVector[i];		
			fprintf(f,"%lf ", element);			
		}
	return true;
}
