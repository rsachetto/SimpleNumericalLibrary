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

#include "Matrix.h"
#include "IdentityMatrix.h"
#include "util.h"

Matrix::Matrix(int numLines, int numCols)
{
	alocateMatrix(numLines, numCols);
	this->numLines = numLines;
	this->numCols = numCols;
}

Matrix::Matrix(const Matrix &m)
{
	if(m.matrix != NULL) {
		numLines = m.numLines;
		numCols = m.numCols;	
		alocateMatrix(m.numLines,m.numCols);
		copyMatrix(matrix, m.matrix);
	}
	else
		throw MatrixException("the source matrix is NULL");
}

Matrix::Matrix(double **mat, int numLines , int numCols ) {

	alocateMatrix(numLines, numCols);
	this->numLines = numLines;
	this->numCols = numCols;

	for (int i = 0; i < numLines; ++i) {
		for (int j = 0; j < numCols; ++j) {
			this->setIJ(i,j,mat[i][j]);			
		}

	}

}

Matrix::~Matrix() {

	freeMatrix();
}

int Matrix::getNumLines() const{
	return numLines;
}

int Matrix::getNumCols() const {
	return numCols;
}

double Matrix::getIJ(int numLine, int numColumn) const {

	if((numLine < this->numLines) && (numColumn < this->numCols))
		return matrix[numLine][numColumn];
	else {
		cerr << "Error getIJ()!! Verify the index"; 
		exit(0);
	}
}

void Matrix::setIJ(int numLine, int numColumn, double value) {	

	if((numLine < this->numLines) && (numColumn < this->numCols))
		matrix[numLine][numColumn] = value;	
	else {
		cerr << "Error setIJ()!! Verify the index" << endl;
		cerr << "I parametro = " << numLine << ", J parametro = " << numColumn << endl;
		cerr << "I matrix = " << this->numLines << ", J matrix = " << this->numCols << endl;
		exit(0);
	}
}

void Matrix::alocateMatrix(int numLines, int numCols) {

	matrix = (double**)calloc(numLines,sizeof(double*));

	if(!matrix){
		cerr << "Error on memory allocation";
		exit(0);
	}

	for(int i = 0; i < numLines; i++){
		matrix[i] = (double*)calloc(numCols,sizeof(double));
		if(!matrix[i]){
			cerr << "Error on memory allocation";
			exit(0);
		}
	}
}

void Matrix::printMatrix() const {

	cout << setw(10);

	for(int i = 0; i < numLines; i++) {
		for(int j = 0; j < numCols; j++) {
			cout << matrix[i][j] << setw(10);	
		}
		cout << endl;
	}

	cout << endl;
}

void Matrix::printMatrix(int p) const {

	cout << setw(2*p);
	cout << setprecision(p) << fixed;
	for(int i = 0; i < numLines; i++) {
		for(int j = 0; j < numCols; j++) {
			cout << matrix[i][j] << setw(2*p);	
		}
		cout << endl;
	}

	cout << endl;
}

void Matrix::copyMatrix(double **dest, double **src) {

	for(int i = 0; i < numLines; i++) {	
		for(int j = 0; j < numCols; j++) {
			dest[i][j] = src[i][j];
		}		
	}
}

const Matrix& Matrix::operator = (const Matrix &m) {

	if(this == &m)
		return *this;

	if(matrix)		
		freeMatrix();

	this->numLines = m.numLines;
	this->numCols = m.numCols;	
	alocateMatrix(m.numLines,m.numCols);
	copyMatrix(matrix, m.matrix);

	return *this;
}

Matrix Matrix::operator + (const Matrix &m) const {
	if((this -> numLines == m.numLines) && (this -> numCols == m.numCols) ) {
		Matrix temp(m.numLines, m.numCols);
		double aux = 0;
		for (int i = 0; i < temp.numLines; ++i) {
			for (int j = 0; j < temp.numCols; ++j) {
				aux = matrix[i][j] + m.matrix[i][j];
				temp.setIJ(i,j,aux); 
			}
		}

		return temp;
	}
	else {
		cerr << "Error na operacao + !! Verifique o numero de linhas e colunas da matrizes"; 
		exit(0);
	}
}

Matrix Matrix::operator - (const Matrix &m) const {

	if((this -> numLines == m.numLines) && (this -> numCols == m.numCols) ) {
		Matrix temp(m.numLines, m.numCols);
		double aux = 0;
		for (int i = 0; i < temp.numLines; ++i) {
			for (int j = 0; j < temp.numCols; ++j) {
				aux = matrix[i][j] - m.matrix[i][j];
				temp.setIJ(i,j,aux); 
			}
		}

		return temp;
	}
	else {
		cerr << "Error na operacao - !! Verifique o numero de linhas e colunas da matrizes"; 
		exit(0);
	}

}

Matrix Matrix::operator * (double a) {

	Matrix temp(this->numLines, this->numCols);
	double aux;
	for (int i = 0; i < temp.numLines; ++i) {
		for (int j = 0; j < temp.numCols; ++j) {
			aux = matrix[i][j]*a;
			temp.setIJ(i,j,aux); 
		}
	}		
	return temp;	
}

Matrix Matrix::operator / (double a) {

	Matrix temp(this->numLines, this->numCols);
	double aux;
	for (int i = 0; i < temp.numLines; ++i) {
		for (int j = 0; j < temp.numCols; ++j) {
			aux = matrix[i][j]/a;
			temp.setIJ(i,j,aux); 
		}
	}		
	return temp;	
}

Matrix Matrix::operator * (const Matrix &m) const {

	if(this -> numCols == m.numLines) {
		Matrix temp(this->numLines, m.numCols);
		double aux = 0;

		for (int i = 0; i < temp.numLines; i++) {
			for (int j = 0; j < temp.numCols; j++) {
				for(int z = 0; z < this->numCols; z++) {							
					aux += (matrix[i][z]*m.matrix[z][j]); 			
				}			
				temp.setIJ(i,j,aux);		
				aux = 0;
			}
		}

		return temp;	
	}
	else {
		cerr << "Error na operacao * !! Verifique se o numero de linhas da primeira matrix = numero de colunas da segunda"; 
		exit(0);
	}
}

Vector Matrix::operator * (Vector &rhs) const {

	if(this->numCols == rhs.getSize()) {
		Vector temp(this->numLines);
		double aux = 0;

		for (int i = 0; i < this->numLines; i++) {
			for (int j = 0; j < this->numCols; j++) {
				aux += (matrix[i][j]*rhs[j]); 				
				temp[i] = aux;				
			}
			aux = 0;
		}

		return temp;	
	}
	else {
		cerr << "Error na operacao * !! Verifique se o numero de linhas da primeira matrix = numero de colunas da segunda"; 
		exit(0);
	}

}

double &Matrix::operator () (int i, int j) {

	if((i < this->numLines) && (j < this->numCols))
		return matrix[i][j];
	else {
		cout << "Matrix Index out of bounds: "<< "i = " << i << ", j = " << j; 
		exit(0);
	}
}

double Matrix::operator () (int i, int j) const {

	if((i < this->numLines) && (j < this->numCols))
		return matrix[i][j];
	else {
		cout << "Matrix Index out of bounds: "<< "i = " << i << ", j = " << j; 
		exit(0);
	}
}

bool Matrix::isNull() {

	if(matrix)
		return false;
	else
		return true;

}

void Matrix::freeMatrix() {

	if(matrix) {
		for(int i = 0; i < numLines; i++) {
			if(matrix[i])
				free(matrix[i]);
		}
		if(matrix)
			free(matrix);
	}		
}

Matrix Matrix::getTransposed() const {

	Matrix temp(numCols, numLines);
	double aux;

	for(int i = 0; i < temp.numLines; i++) {
		for(int j = 0; j < temp.numCols; j++) {
			aux = matrix[j][i];
			temp.setIJ(i,j,aux);			
		}
	}

	return temp;
}

Matrix Matrix::getInverseMatrix() {

	Matrix A(*this);
	int n = A.getNumCols()-1;
	double aux;
	Vector *b;
	Matrix inverse(n+1,0);

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


	for(int i = 0; i <= n; i++) {

		b = new Vector(n+1);
		b->setI(i,1);

		for(int i = 1; i <= n; i++) {
			for(int k = 0; k <= i-1; k++) {	
				if(i == k)
					aux = 1;
				else
					aux = A(i,k);	
				b->setI(i,(b->getI(i)-aux*b->getI(k)));
				//b[i] = &b[i]-aux*(b[k]);		
			}
		}	

		*b = retroSubstitution(A,*b);
		inverse.addColumnRight(*b);

	}	

	return inverse;
}

double Matrix::getDeterminant() const {


	Matrix aux(this->getDiagonalMatrix());

	if(numLines != numCols) {
		cerr << "Error getDeterminant()!! ";
		exit(0);
	}

	double det = 1;		

	for(int j = 0; j < numCols; j++) { 
		det *= aux.getIJ(j,j);;						
	}	
	return det;	
}

Matrix Matrix::reducedMatrix (Matrix &m, int offset) {


	Matrix aux(m.numLines-1, m.numCols-1);
	int iAux = 0, jAux = 0;


	for(int i = 1; i < m.numLines; i++) {
		for(int j = 0; j < m.numCols; j++) {

			if(j != offset) {
				aux.setIJ(iAux,jAux,m.matrix[i][j]);	
				jAux += 1;					
			}

			if((jAux == aux.getNumCols()) && (iAux < aux.getNumLines())){
				iAux+=1;
				jAux = 0;								
			}							
		}	

	}

	return aux;

}

Vector Matrix::getPrincipalDiagonal() const {

	if(numLines != numCols) {
		cerr << "Error getPrincipalDiagonal()!! " << endl;
		exit(0);
	}

	Vector result(numLines);
	double aux;

	for(int i = 0; i < numLines; i++) {

		aux = matrix[i][i];
		result.setI(i, aux);

	}

	return result;

}

Vector Matrix::getSecundaryDiagonal() const {

	if(numLines != numCols) {
		cerr << "Error getSecundaryDiagonal()!! " << endl;
		exit(0);
	}

	Vector result(numLines);
	double aux;

	for(int i = numLines-1; i >= 0; i--) {

		aux = matrix[(numLines-1)-i][i];
		result.setI((numLines-1)-i, aux);

	}

	return result;

}

Vector Matrix::getLine(int linha) const {

	Vector result(this->numCols);

	for(int i = 0; i < numCols; i++) {

		result.setI(i,matrix[linha][i]);  

	}

	return result;

}

Vector Matrix::getColumn(int coluna) const {

	Vector result(this->numLines);

	for(int i = 0; i < numLines; i++) {

		result.setI(i,matrix[i][coluna]);  

	}

	return result;

}

void Matrix::readMatrixFromFile(const string fileName) {

	MatrixFileHandler fileHandler;
	fileHandler.readMatrixFile(fileName, matrix, numLines, numCols);	

}

void Matrix::writeMatrixToFile(const string fileName) const {

	MatrixFileHandler fileHandler;
	fileHandler.writeMatrixFile(fileName, matrix, numLines, numCols);	

}

void Matrix::multiplyLineEscalar(int line, double escalar){

	for (int j = 0; j < numCols; ++j) {
		matrix[line][j] = matrix[line][j]*escalar;
	}

}

void Matrix::sumLines(int destLine, int line1, int line2) {
	for (int j = 0; j < numCols; ++j) {
		matrix[destLine][j] = matrix[line1][j]+matrix[line2][j];
	}
}

void Matrix::changeLines(int line1, int line2) {

	double aux;

	for (int j = 0; j < numCols; ++j) {		
		aux = matrix[line1][j];
		matrix[line1][j] = matrix[line2][j];
		matrix[line2][j] = aux;
	}
}

void Matrix::sumLineVector(int linhaDestino,int line, Vector v) {

	for (int j = 0; j < numCols; ++j) {
		matrix[linhaDestino][j] = matrix[line][j] + v[j];
	}

}

void Matrix::addColumnRight(Vector &b) {

	if(numLines != b.getSize()) {
		cerr << "Error addColumn(Vector &b)!! " << endl;
		exit(0);
	}

	Matrix aux(numLines, numCols+1);
	double val;

	for (int i = 0; i < numLines; ++i) {
		for (int j = 0; j < numCols; ++j) {
			val = this->getIJ(i,j);
			aux.setIJ(i,j, val);
		}	
	}

	for (int i = 0; i < numLines; ++i) {
		val = b.getI(i);
		aux.setIJ(i,numCols, val);
	}	

	*this = aux;

}

void Matrix::addColumnLeft(Vector &b) {

	if(numLines != b.getSize()) {
		cerr << "Error addColumn(Vector &b)!! " << endl;
		exit(0);
	}

	Matrix aux(numLines, numCols+1);
	double val;

	for (int i = 0; i < numLines; ++i) {
		for (int j = 0; j < numCols; ++j) {
			val = this->getIJ(i,j);
			aux.setIJ(i,j+1, val);
		}	
	}

	for (int i = 0; i < numLines; ++i) {
		val = b.getI(i);
		aux.setIJ(i,0, val);
	}	

	*this = aux;

}

Matrix Matrix::getDiagonalMatrix() const {

	Matrix diagonal(*this);
	Vector aux(numCols);
	int numLine = 0;
	double m, mAux;

	while (numLine < numLines) {

		for(int k = numLine; k < numLines-1; k++) {

			mAux = diagonal.getIJ(numLine,numLine);

			if(mAux != 0) {
				m = diagonal.getIJ(k+1,numLine)/mAux;				
				aux = diagonal.getLine(numLine)*(-m);
				diagonal.sumLineVector(k+1, k+1,aux);		
			}
		}

		for(int k = numLine; k > 0; k--) {

			mAux = diagonal.getIJ(numLine,numLine);

			if(mAux != 0) {
				m = diagonal.getIJ(k-1,numLine)/mAux;
				aux = diagonal.getLine(numLine)*(-m);
				diagonal.sumLineVector(k-1, k-1,aux);
			}
		}				
		numLine++;			
	}
	return diagonal;			
}

Vector Matrix::retroSubstitution(Matrix A, Vector b) {

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

double Matrix::getNorm() const{

	double sum = 0;

	for(int i = 0; i < this->getNumLines(); i++) {
		for(int j = 0; j < this->getNumCols(); j++) {
			sum += pow(this->getIJ(i,j),2);
		}
	}

	return sqrt(sum);

}

double Matrix::getCond() {

	return  (this->getNorm())*(this->getInverseMatrix().getNorm());

}

