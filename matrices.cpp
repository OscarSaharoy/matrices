#include "matrices.h"
#include <ctime>
#include <cstdlib>
#include <stdexcept>
#include <cmath>

Matrix::Matrix() : 
	r(0), 
	c(0), 
	l(0) {

	// construct a blank matrix

	// assign array for elements in matrix
	mp = new double[0]; // memory pointer
};

Matrix::Matrix(const int& r, const int& c) : 
	r(r), 
	c(c), 
	l(r*c) {

	// check that number of rows and columns is positive
	if(r < 0 || c < 0) {

		throw std::invalid_argument("Both args must be positive (number of rows and number of columns in matrix)");
	}

	// assign array for elements in matrix
	mp = new double[l];
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> init_list) {

	// contruct matrix from nested initialiser lists

	// set dimensions of matrix
	r = init_list.end() - init_list.begin();
	c = init_list.begin() -> size();
	l = r*c;

	// intialise iterators
	auto iter1 = init_list.begin();
	auto iter2 = (*iter1).begin();
	auto end2  = (*iter1).end();

	// allocate array
	mp = new double[l];

	// populate array
	for(int t = 0; t < l; ++t) {

		mp[t] = *iter2;

		++iter2;

		// set iterators values to allow 1-dimensional loop to handle 2d initialiser lists
		if(iter2 == end2) {

			++iter1;
			iter2 = (*iter1).begin();
			end2  = (*iter1).end();
		}
	}
}

Matrix::Matrix(const Matrix& M) : 
	r(M.r), 
	c(M.c), 
	l(M.l) {

	// copy constructor

	// allocate array
	mp = new double[M.l];

	// use memcpy to copy elements of orignal array into this one
	memcpy(mp, M.mp, sizeof(double) * l);
}

Matrix::~Matrix() {

	// destructor

	// deallocate array
	delete[] mp;
}


Matrix Matrix::dims() const {

	// return matrix containing dimensions of matrix as [[   rows], 
	//													 [columns]]

	return Matrix( {{(double)r}, {(double)c}} );
}

Matrix Matrix::T() const {

	// returns transposed matrix

	// create new array to fill with transposed elements
	Matrix temp(c, r);

	// set indices
	int i=0;
	int j=0;

	for(int t=0; t < l; ++t) {

		// set element in temp matrix
		temp.mp[t] = mp[j*c + i];

		// update indices
		++j;

		if(j==temp.c) {

			++i;
			j = 0;
		}
	}

	return temp;
}

double Matrix::det() const {

	// find the determinant of matrix - must be square

	// temporary variable to store output
	double temp = 0;

	// if its a 2x2 matrix, return the determinant directly
	if(r == 2 && c == 2) {

		return mp[0] * mp[3] - mp[1] * mp[2];
	}

	// loop over first row, recursively calling det on the minors
	for(int t=0; t<c; ++t) {

		temp += (-t%2 * 2 + 1) * mp[t] * minor(0, t).det();
	}

	return temp;
}

Matrix Matrix::inv() const {

	// invert matrix - square matrices only
	Matrix temp(r, c);
	Matrix transpose = T();
	double determinant = det();

	// if its a 2x2 matrix return the inverse
	if(r == 2 && c == 2) {

		temp.mp[0] =   1/determinant * mp[3];
		temp.mp[1] = - 1/determinant * mp[1];
		temp.mp[2] = - 1/determinant * mp[2];
		temp.mp[3] =   1/determinant * mp[0];

		return temp;
	}

	for(int t=0; t<l; ++t) {

		temp.mp[t] = (-t%2 * 2 + 1)/determinant * transpose.minor(t).det();
	}

	return temp;
}

Matrix Matrix::normalised() const {

	Matrix temp(r, c);

	double mag = magnitude();

	for(int t=0; t<l; ++t) {

		temp.mp[t] = mp[t] / mag;
	}

	return temp;
}

void Matrix::for_each(double func(double&)) {

	// applies func to each element in matrix

	for(int t=0; t<l; ++t) {

		// set each element in the oputput matrix as that of the original matrix with func applied
		mp[t] = func(mp[t]);
	}
}

Matrix Matrix::squared() const {

	// square each element in the matrix

	// create result matrix
	Matrix temp(r, c);

	for(int t=0; t<l; ++t) {

		// set each element in the output matrix as the square of that in the original matrix
		temp.mp[t] = mp[t] * mp[t];
	}

	return temp;
}

Matrix Matrix::sqrt() const {

	// square root each element in the matrix

	// create result matrix
	Matrix temp(r, c);

	for(int t=0; t<l; ++t) {

		// set each element in the output matrix as the square root of that in the original matrix
		temp.mp[t] = std::pow(mp[t], 0.5);
	}

	return temp;
}

double Matrix::sum() const {

	// return the sum of all the elements in the matrix

	// create temp double to store the sum
	double temp = 0;

	for(int t=0; t<l; ++t) {

		// add each element to temp before returning it
		temp += mp[t];
	}

	return temp;
}

double Matrix::magnitude() const {

	double temp = 0;

	for(int t=0; t<l; ++t) {

		temp += mp[t]*mp[t];
	}

	return std::sqrt(temp);
}

void Matrix::randomise(const double& mean, const double& range) {

	for(int i = 0; i < l; ++i) {

		// set each element to a random float with given mean and range using cstdlib random function 
		mp[i] = (double)std::rand() / RAND_MAX * range + mean - range/2;
	}
}

void Matrix::identity() {

	for(int i = 0; i < l; ++i) {

		// set each element to 1 or 0 depending on position
		mp[i] = i%c == i/r ? 1 : 0;
	}
}

Matrix& Matrix::operator=(const Matrix& other) {

	// assignment operator

	// set rows, columns and length to that of the other matrix
	r = other.r;
	c = other.c;
	l = r*c;

	// delete and reallocate array for elements
	delete[] mp;
	mp = new double[l];

	// copy elements from other matrix into this one's array
	memcpy(mp, other.mp, sizeof(double) * l);

	// im not sure how this works or what it does
	return *this;
}


Matrix Matrix::operator+(const Matrix& other) const {

	// add 2 matrices element wise

	// create result matrix
	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		// element in result matrix is sum of that in 2 input matrices
		temp.mp[i] = mp[i] + other.mp[i];
	}

	return temp;
}

Matrix Matrix::operator-(const Matrix& other) const {

	// subtract 2 matrices element wise

	// create result matrix
	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		// element in result matrix is element in first matrix minus element in other matrix
		temp.mp[i] = mp[i] - other.mp[i];
	}

	return temp;
}


Matrix Matrix::operator*(const Matrix& other) const {

	// check dimensons are valid for matrix multiply
	if(c != other.r) {

		throw std::invalid_argument("Dimensions invalid for matrix multiply.");
	}

	// initialise result matrix
	Matrix temp(r, other.c);

	// indices of entry in result matrix
	int i = 0;
	int j = 0;

	for(int t=0; t < temp.l; ++t) {

		// sum for entry in result matrix
		double sum = 0;

		// calculate entry
		for(int k=0; k < other.r; k++) {

			sum += mp[i*c + k] * other.mp[k*other.c + j];
		}

		// assign entry to value of sum
		temp.mp[t] = sum;

		// update indices
		++j;

		if(j==temp.c) {

			++i;
			j = 0;
		}
	}

	return temp;
}

Matrix Matrix::operator+(const double& other) const {

	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		temp.mp[i] = mp[i] + other;
	}

	return temp;
}

Matrix Matrix::operator-(const double& other) const {

	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		temp.mp[i] = mp[i] - other;
	}

	return temp;
}

Matrix Matrix::operator*(const double& other) const {

	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		temp.mp[i] = mp[i] * other;
	}

	return temp;
}

Matrix Matrix::operator/(const double& other) const {

	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		temp.mp[i] = mp[i] / other;
	}

	return temp;
}

Matrix Matrix::operator&&(const Matrix& other) const {

	Matrix temp(r, c);

	for(int t=0; t < l; ++t) {

		temp.mp[t] = mp[t] && other.mp[t];
	}

	return temp;
}

Matrix Matrix::operator||(const Matrix& other) const {

	Matrix temp(r, c);

	for(int t=0; t < l; ++t) {

		temp.mp[t] = mp[t] || other.mp[t];
	}

	return temp;
}

Matrix Matrix::operator!() const {
	
	Matrix temp(r, c);

	for(int t=0; t < l; ++t) {

		temp.mp[t] = !mp[t];
	}

	return temp;
}

void Matrix::operator+=(const Matrix& other) {

	for(int t=0; t < l; ++t) {

		mp[t] += other.mp[t];
	}
}

void Matrix::operator-=(const Matrix& other) {

	for(int t=0; t < l; ++t) {

		mp[t] -= other.mp[t];
	}
}

void Matrix::operator*=(const Matrix& other) {

	if(r != other.c) {

		throw std::invalid_argument("Dimensions invalid for matrix multiply.");
	}

	// initialise result array
	Matrix temp(other.r, c);

	// indices of entry in result matrix
	int i = 0;
	int j = 0;

	for(int t=0; t < temp.l; ++t) {

		// sum for entry in result matrix
		double sum = 0;

		// calculate entry
		for(int k=0; k < r; k++) {

			sum += other.mp[i*other.c + k] * mp[k*c + j];
		}

		// assign entry to value of sum
		temp.mp[t] = sum;

		// update indices
		++j;

		if(j==temp.c) {

			++i;
			j = 0;
		}
	}

	// overwrite this matrix with result of multiplication
	*this = temp;
}

void Matrix::operator+=(const double& other) {

	for(int i=0; i<l; ++i) {

		mp[i] += other;
	}
}

void Matrix::operator-=(const double& other) {

	for(int i=0; i<l; ++i) {

		mp[i] -= other;
	}
}

void Matrix::operator*=(const double& other) {

	for(int i=0; i<l; ++i) {

		mp[i] *= other;
	}
}

void Matrix::operator/=(const double& other) {

	for(int i=0; i<l; ++i) {

		mp[i] /= other;
	}
}

double& Matrix::operator()(const int& i, const int& j) {

	// return reference to element at index i, j
	return mp[i*c + j];
}

double& Matrix::operator[](const int& t) {

	// return reference to element at index t
	return mp[t];
}

Matrix Matrix::slice(const int& r1, const int& c1, const int& r2, const int& c2) const {

	// check indices are valid to make a slice
	if(r2 < r1 || c2 < c1) {

		throw std::invalid_argument("index of bottom left corner (2nd 2 args) must be greater than index of first corner (first 2 args)");
	}

	// check indices are inside matrix
	if(r1 > r || c1 > c) {

		throw std::invalid_argument("indices must all be inside matrix");
	}

	// create matrix with correct number of rows and columns
	Matrix temp(r2-r1, c2-c1);

	// indices used to transfer elements into temp matrix
	int i = 0;
	int j = 0;

	for(int t=0; t < temp.l; ++t) {

		// set element in temp
		temp.mp[t] = mp[(i + r1) * c + (j + c1)];

		// update indices
		++j;

		if(j==temp.c) {

			++i;
			j = 0;
		}
	}

	return temp;
}

Matrix Matrix::minor(const int& r1, const int& c1) const {

	// temporary matrix to store result
	Matrix temp(r-1, c-1);

	// index of current element in the minor matrix
	int im = 0;

	// iterate over elements of temp matrix
	for(int t=0; t < temp.l; ++t, ++im) {

		// skip if im is inside the crossed out row/column
		while(im%c == c1 || (im/c)%r == r1) {

			++im;
		}

		temp.mp[t] = mp[im];
	}

	return temp;
}

Matrix Matrix::minor(const int& t1) const {

	// just find row and column of t1 and call other function
	int r1 = (t1/c)%r;
	int c1 = t1%c;

	return minor(r1, c1);
}


std::ostream& operator<<(std::ostream& os, const Matrix& M) {

	// print matrix

	// double open brackets at start
	os << "[[";

	for(int i=0; i < M.l-1; ++i) {

		// print element up to but not including last element in matrix
		os << M.mp[i];

		// print comma between elements or close bracket, comma, newline and open bracket after last element in row
		os << ( (i+1) % M.c != 0 ? ", " : "],\n [");
	}

	// different for last element so print it here
	os << M.mp[M.l-1];

	// double close bracket at end
	os << "]]";

	return os;
}

Matrix operator*(const double& d, const Matrix& M) {

	// multiply each element in matrix by double d
	Matrix temp(M.r, M.c);

	for(int i = 0; i < M.l; ++i) {

		temp.mp[i] = M.mp[i] * d;
	}

	return temp;
}

Matrix array_mult(const Matrix& M1, const Matrix& M2) {

	// multiply 2 matrices element wise

	// check matrices dimesions are the same
	if(M1.c != M2.c || M1.r != M2.r) {

		throw std::invalid_argument("Dimensions invalid for array multiply.");
	}

	// temp matrix to store output
	Matrix temp(M1.r, M1.c);

	for(int t=0; t < M1.l; ++t) {

		// output element is product of 2 elements in input matrices
		temp.mp[t] = M1.mp[t] * M2.mp[t];
	}

	return temp;
}

Matrix array_div(const Matrix& M1, const Matrix& M2) {

	// divide 2 matrices element wise

	// check matrices dimesions are the same
	if(M1.c != M2.c || M1.r != M2.r) {

		throw std::invalid_argument("Dimensions invalid for array multiply.");
	}

	// temp matrix to store output
	Matrix temp(M1.r, M1.c);

	for(int t=0; t < M1.l; ++t) {

		// output element is element in M1 divided by element in M2
		temp.mp[t] = M1.mp[t] / M2.mp[t];
	}

	return temp;
}


Matrix join(const Matrix& M1, const Matrix& M2, const int& axis) {

	Matrix temp;

	if(axis == 0) {

		// initialise output matrix
		int r = M1.r + M2.r;
		int c = M1.c;
		temp = Matrix(r, c);

		// loop over elements of temp, filling with the correct element of M1 or M2
		for(int t=0; t<temp.l; ++t) {

			temp.mp[t] = t < M1.l ? M1.mp[t] : M2.mp[t - M1.l];
		}
	}

	else { // axis == 1

		// initialise output matrix
		int r = M1.r;
		int c = M1.c + M2.c;
		temp = Matrix(r, c);

		// loop over elements of temp, filling with the correct element of M1 or M2
		for(int t=0; t<temp.l; ++t) {

			temp.mp[t] = t%c < M1.c ? M1.mp[t - t/c * M2.c] : M2.mp[t - (t/c + 1) * M1.c];
		}
	}

	return temp;
}

double dot(const Matrix& V1, const Matrix& V2) {

	// dot product of 2 vectors but also works on rectangular matrices

	double temp = 0;

	// for each element of V1, increment temp by its product with the corresponding element in V2
	for(int t=0; t < V1.l; ++t) {

		temp += V1.mp[t] * V2.mp[t];
	}

	return temp;
}

Matrix eigenvectors(const Matrix& M) {

	// approximate eigenvectors of M if M is square
	// expensive function

	// stores output matrix
	Matrix temp;
	Matrix inverse = M.inv();

	// randomise min and max eigenvectors initially
	Matrix umax(M.r, 1);
	Matrix umin(M.r, 1);
	umax.randomise(0.0f, 1.0f);
	umin.randomise(0.0f, 1.0f);

	// iteration to bring random vector toward max eigenvector
	for(int t=0; t<50; t++) {
		
		umax = (M * umax).normalised();
	}

	// get eigenvalue
	double lmax = (M * umax).mp[0] / umax.mp[0];

	// iteration to bring random vector toward min eigenvector
	for(int t=0; t<50; t++) {
		
		umin = (inverse * umin).normalised();
	}

	umax = umax / umax.mp[0];
	umin = umin / umin.mp[0];

	// get eigenvalue
	double lmin = (M * umin).mp[0] / umin.mp[0];

	temp = umin;

	for(int t=0; t<M.c-2; ++t) {

		// guess eigenvalue - interpolate between lmax and lmin
		double lguess = lmin + (lmax - lmin) * (float)(t+1) / (float)(M.c-1);
		Matrix ident = Matrix(M.r, M.c);
		ident.identity();

		// modified matrix to use to approach eigenvalue
		Matrix Mmod = (lguess * ident - M).inv();

		// eigenvector - start random
		Matrix u(M.r, 1);
		u.randomise(0, 1);

		// iterate to get eigenvector
		for(int t=0; t<50; t++) {
		
			u = (Mmod * u).normalised();
		}

		// fix u to have first element 1 and join to temp
		u = u / u.mp[0];
		temp = join(temp, u, 1);
	}

	// join umax to temp and return
	return join(temp, umax, 1);
}

Matrix eigenvalues(const Matrix& M, const Matrix& eigM) {

	// get eigenvalues from matrix and eigenvectors matrix

	// stores output matrix
	Matrix temp(eigM.r, eigM.c);
	temp.identity();

	// loop over diagonal elements of temp, filling them with the eigenvalues
	for(int t=0; t<temp.c; ++t) {
	
		Matrix u = eigM.slice(0, t, eigM.r, t+1);
		temp(t, t) = (M * u).mp[0] / u.mp[0];
	}

	return temp;
}

Matrix cross(const Matrix& V1, const Matrix& V2) {

	// cross product - only defined for 3x1 matrices!!
	Matrix temp(3, 1);

	temp.mp[0] = V1.mp[1] * V1.mp[2] - V1.mp[2] * V2.mp[1];
	temp.mp[1] = V1.mp[2] * V1.mp[0] - V1.mp[0] * V2.mp[2];
	temp.mp[2] = V1.mp[0] * V1.mp[1] - V1.mp[1] * V2.mp[0];

	return temp;
}

double det(const Matrix& M) {

	// find the determinant of a matrix - must be square

	// temporary variable to store output
	double temp = 0;

	// if its a 2x2 matrix, return the determinant directly
	if(M.r == 2 && M.c == 2) {

		return M.mp[0] * M.mp[3] - M.mp[1] * M.mp[2];
	}

	// loop over first row, recursively calling det on the minors
	for(int t=0; t<M.c; ++t) {

		temp += (-t%2 * 2 + 1) * M.mp[t] * det(M.minor(0, t));
	}

	return temp;
}