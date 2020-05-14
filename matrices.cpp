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


std::vector<int> Matrix::dims() {

	// return vector containing dimensions of matrix as (rows, columns)
	return std::vector<int>({r, c});
}

Matrix Matrix::T() {

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

double Matrix::det() {

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

Matrix Matrix::inv() {

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

Matrix Matrix::normalised() {

	Matrix temp(r, c);

	double mag = magnitude();

	for(int t=0; t<l; ++t) {

		temp.mp[t] = mp[t] / mag;
	}

	return temp;
}

Matrix Matrix::for_each(double func(double&)) {

	// applies func to each element in matrix

	// create result matrix
	Matrix temp(r, c);

	for(int t=0; t<l; ++t) {

		// set each element in the oputput matrix as that of the original matrix with func applied
		temp.mp[t] = func(mp[t]);
	}

	return temp;
}

Matrix Matrix::squared() {

	// square each element in the matrix

	// create result matrix
	Matrix temp(r, c);

	for(int t=0; t<l; ++t) {

		// set each element in the output matrix as the square of that in the original matrix
		temp.mp[t] = mp[t] * mp[t];
	}

	return temp;
}

Matrix Matrix::sqrt() {

	// square root each element in the matrix

	// create result matrix
	Matrix temp(r, c);

	for(int t=0; t<l; ++t) {

		// set each element in the output matrix as the square root of that in the original matrix
		temp.mp[t] = std::pow(mp[t], 0.5);
	}

	return temp;
}

double Matrix::sum() {

	// return the sum of all the elements in the matrix

	// create temp double to store the sum
	double temp = 0;

	for(int t=0; t<l; ++t) {

		// add each element to temp before returning it
		temp += mp[t];
	}

	return temp;
}

double Matrix::magnitude() {

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


Matrix Matrix::operator+(const Matrix& other) {

	// add 2 matrices element wise

	// create result matrix
	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		// element in result matrix is sum of that in 2 input matrices
		temp.mp[i] = mp[i] + other.mp[i];
	}

	return temp;
}

Matrix Matrix::operator-(const Matrix& other) {

	// subtract 2 matrices element wise

	// create result matrix
	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		// element in result matrix is element in first matrix minus element in other matrix
		temp.mp[i] = mp[i] - other.mp[i];
	}

	return temp;
}


Matrix Matrix::operator*(const Matrix& other) {

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

Matrix Matrix::operator+(const double& other) {

	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		temp.mp[i] = mp[i] + other;
	}

	return temp;
}

Matrix Matrix::operator-(const double& other) {

	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		temp.mp[i] = mp[i] - other;
	}

	return temp;
}

Matrix Matrix::operator*(const double& other) {

	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		temp.mp[i] = mp[i] * other;
	}

	return temp;
}

Matrix Matrix::operator/(const double& other) {

	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		temp.mp[i] = mp[i] / other;
	}

	return temp;
}

Matrix Matrix::operator&&(const Matrix& other) {

	Matrix temp(r, c);

	for(int t=0; t < l; ++t) {

		temp.mp[t] = mp[t] && other.mp[t];
	}

	return temp;
}

Matrix Matrix::operator||(const Matrix& other) {

	Matrix temp(r, c);

	for(int t=0; t < l; ++t) {

		temp.mp[t] = mp[t] || other.mp[t];
	}

	return temp;
}

Matrix Matrix::operator!() {
	
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

	if(c != other.r) {

		throw std::invalid_argument("Dimensions invalid for matrix multiply.");
	}

	// initialise result array
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

	int r, c;

	if(axis == 0) {
		if(M1.c != M2.c) {

			throw std::invalid_argument("Matrices should have identical number of columns");
		}

		c = M1.c;
		r = M1.r + M2.r;
	}

	else if(axis == 1) {
		if(M1.r != M2.r) {

			throw std::invalid_argument("Matrices should have identical number of rows");
		}

		c = M1.c + M2.c;
		r = M1.r;
	}

	else {
		throw std::invalid_argument("axis must be equal to 0 or 1");
	}

	Matrix temp(r, c);

	if(axis == 1) {

		for(int t=0; t<temp.l; ++t) {

			temp.mp[t] = t%c < M1.c ? M1.mp[t - t/c * M2.c] : M2.mp[t - (t/c + 1) * M1.c];
		}

		return temp;
	}

	int cr = -1; // current row

	for(int t=0; t < M1.l; ++t) {

		if(t % M1.c == 0) {

			++cr;
		}

		temp.mp[t + cr*M1.c] = M1.mp[t];
	}

	cr = -1; // current row

	for(int t=0; t < M2.l; ++t) {

		if(t % M2.c == 0) {

			++cr;
		}

		temp.mp[t + cr*(M2.c - 1) + M1.c] = M2.mp[t];

	}

	return temp;
}

Matrix join(const std::vector<Matrix*>& Ms, const int& axis) {

	int r = 0;
	int c = 0;

	if(axis == 0) {

		for(auto X : Ms) {

			if(X->l != Ms[0]->l) {
				throw std::invalid_argument("Matrices should have identical number of columns");
			}

			r += X->r;
		}

		c = Ms[0]->c;
	}

	else if(axis == 1) {

		for(auto X : Ms) {
			
			if(X->l != Ms[0]->l) {
				throw std::invalid_argument("Matrices should have identical number of columns");
			}

			c += X->c;
		}

		r = Ms[0]->r;
	}

	else {
		throw std::invalid_argument("axis must be equal to 0 or 1");
	}

	Matrix temp(r, c);

	int c_sum = 0; // sum of columns for offset in array

	for(auto X : Ms) {

		int cr = -1; // current row

		for(int t=0; t < X->l; ++t) {

			if(t % X->c == 0) {

				++cr;
			}

			temp.mp[t + cr*(temp.c - 1) + c_sum] = X->mp[t];
		}

		c_sum += X->c;
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

Matrix eigenvalues(Matrix& M) {

	// approximate eigenvalues of M if M is square

	// stores output matrix
	Matrix temp(M.r, M.c);
	Matrix eigenMatrix = eigenvectors(M);
	Matrix u1 = eigenMatrix.slice(0,0,2,1);
	Matrix u2 = eigenMatrix.slice(0,1,2,2);

	temp.mp[0] = (M * u1).mp[0] / u1.mp[0];
	temp.mp[1] = 0;
	temp.mp[2] = 0;
	temp.mp[3] = (M * u2).mp[0] / u2.mp[0];

	return temp;
}

Matrix eigenvectors(Matrix& M) {

	// approximate eigenvectors of M if M is square

	// stores output matrix
	Matrix temp;
	Matrix inverse = M.inv();

	Matrix u1(M.r, 1);
	Matrix u2(M.r, 1);
	Matrix u3(M.r, 1);
	u1.randomise(0.0f, 1.0f);
	u2.randomise(0.0f, 1.0f);
	u3.randomise(0.0f, 1.0f);

	for(int t=0; t<20; t++) {
		
		u1 = (M * u1).normalised();
	}

	for(int t=0; t<20; t++) {
		
		u2 = (inverse * u2).normalised();
	}

	u1 = u1 / u1.mp[0];
	u2 = u2 / u2.mp[0];

	return join(u1, u2, 1);
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