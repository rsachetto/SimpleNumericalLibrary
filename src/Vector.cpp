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

#include "Vector.h"
#include "VectorFileHandler.h"

Vector::Vector(int size) {	
	this->size = size;
	alocateVector();
}

Vector::Vector(const Vector &copy) {	
	
	if(copy.vector != NULL) {
		this->size = copy.size;
		alocateVector();
		copyVector(vector, copy.vector);
	}
	else
		throw VectorException("the source vector is NULL");
}

Vector::Vector(double *vec, int size) {	
	
	if(vec != NULL) {
		this->size = size;
		alocateVector();
			
		for(int i = 0; i < size; i++)
			this->setI(i,vec[i]);
	}
	else 
		throw VectorException("the source vector is NULL");
}

Vector::~Vector(){
	if(vector)
		freeVector();

}	

int Vector::getSize() const {

	return this->size;

}

double Vector::getI( int index ) const {
	if(index < this->size)
		return vector[index];
	else {
		cout << "Vector index out of bounds"<< index << endl;
		exit(0);
	}
		
}

void Vector::setI( int index, double value ){
	if(index < this->size)
		vector[index] = value;
	else {
		cout << "Vector index out of bounds" << index << endl;
		exit(0);
	}
}

const Vector& Vector::operator = (const Vector & rhs){
	
	if(this == &rhs)
		return *this;
	
	if(vector)
		freeVector();
	
	this->size = rhs.size;	
	alocateVector();
	copyVector(vector, rhs.vector);
	
	return *this;
}

Vector Vector::operator + (const Vector & rhs) const {
	
	Vector temp(rhs.size);
	double aux = 0;
	for (int i = 0; i < temp.size; ++i) {		
			aux = vector[i] + rhs.vector[i];
			temp.setI(i,aux); 		
	}		
	
	return temp;	
}

Vector Vector::operator - (const Vector & rhs) const{
	
	Vector temp(rhs.size);
	double aux = 0;
	for (int i = 0; i < temp.size; ++i) {		
			aux = vector[i] - rhs.vector[i];
			temp.setI(i,aux); 		
	}		
	
	return temp;	
}

Vector Vector::operator - () const {

	Vector temp(this->size);
	double aux;
	for (int i = 0; i < temp.size; ++i) {		
			aux = -vector[i];
			temp.setI(i,aux); 		
	}	

	return temp;

}

Vector Vector::operator * (const Vector & rhs) const{
	
	Vector temp(rhs.size);
	double aux = 0;
	for (int i = 0; i < temp.size; ++i) {		
			aux = vector[i] * rhs.vector[i];
			temp.setI(i,aux); 		
	}		
	
	return temp;	
}

Vector Vector::operator * (double escalar){
	
	Vector temp(this->size);
	double aux = 0;
	for (int i = 0; i < temp.size; ++i) {		
			aux = vector[i]*escalar;
			temp.setI(i,aux); 		
	}		
	
	return temp;	
}

Vector Vector::operator / (double escalar){
	
	Vector temp(this->size);
	double aux = 0;
	for (int i = 0; i < temp.size; ++i) {		
			aux = vector[i]/escalar;
			temp.setI(i,aux); 		
	}		
	
	return temp;	
}

double &Vector::operator [] (int i)  {

	if(i < this->size)
		return vector[i];
	else {
		cout << "Vector index out of bounds: " << i << endl;
		exit(0);
	}

}


double Vector::operator [] (int i) const {

	if(i < this->size)
		return vector[i];
	else {
		cout << "Vector index out of bounds: " << i << endl;
		exit(0);
	}

}

double Vector::innerProduct(Vector a, Vector b) {

	double aux;
	aux = 0;
	
	for(int i = 0; i < a.getSize(); i++)
		aux += a[i]*b[i];
		
	return aux; 
}

bool Vector::isNull(){

	if(vector)
		return false;
	else
		return true;

}

void Vector::printVector() const {
	for (int i = 0; i < size; ++i) {		
			if(i != size-1)
				cout << vector[i] << ", ";		
			else
				cout << vector[i];
	}			
}

double Vector::norm2() const {

	double aux;
	aux = 0;
	
	for(int i = 0; i < this->getSize(); i++)
		aux += pow(this->getI(i),2);
		
	return sqrt(aux);

}

void Vector::alocateVector() {		
	vector = (double*)calloc(size, sizeof(double));
}

void Vector::freeVector() {
	free(vector);
}

void Vector::copyVector(double *dest, double *src) {
	
	for(int i = 0; i < size; i++) {
	
		dest[i] = src[i];
	
	}

}

double Vector::multipilyComponents() {
	
	double result = 1;
	
	for(int i = 0; i < size; i++) {
	
		result *= vector[i];
	
	}
	
	return result;	
}

void Vector::readVectorFromFile(const string fileName){
	
	VectorFileHandler fileHandler;
	fileHandler.readVectorFile(fileName, vector, size);

}
void Vector::writeVectorToFile(const string fileName) const {
	

	VectorFileHandler fileHandler;
	fileHandler.writeVectorFile(fileName, vector, size);

}

ostream &operator<<( ostream &output, const Vector &a ) {
	int i;

    for ( i = 0; i < a.size; i++ ) {
    	output << setw( 10 ) << a.vector[ i ];

        if ( ( i + 1 ) % 10 == 0 ) // 4 numbers per row of output
        	output << endl;
     }

     if ( i % 10 != 0 ) // end last line of output
     	output << endl;

     return output;
}

