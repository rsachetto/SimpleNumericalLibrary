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

#ifndef _IDENTITYMATRIX_H_
#define _IDENTITYMATRIX_H_

#include "Matrix.h"

//! An nxn Identity Matrix.
/*!
A class that represents a nxn Identity Matrix of doubles
*/

class IdentityMatrix:public Matrix
{
public:
	//! Constructor.
	/*!
		 Construct a Identity Matrix numLinesXnumCols
		\param numLines lines number of the matrix
		\param numCols columns number of the matrix
	*/
	IdentityMatrix(int numLines = 1,int numCols = 1);
	
};

#endif //_IDENTITYMATRIX_H_
