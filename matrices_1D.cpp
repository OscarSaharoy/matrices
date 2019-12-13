#include "matrices_1D.h"
#include <ctime>
#include <cstdlib>
#include <stdexcept>

Matrix::Matrix() : 
	r(0), 
	c(0), 
	l(0) {

	mp = new double[0];
};

Matrix::Matrix(const int& r, const int& c) : 
	r(r), 
	c(c), 
	l(r*c) {

	if(r < 0 || c < 0) {

		throw std::invalid_argument("Both args must be positive (number of rows and number of columns in matrix)");
	}

	mp = new double[l];
}

Matrix::Matrix(std::initializer_list<std::vector<double>> init_list) {

	r = init_list.end() - init_list.begin();
	c = init_list.begin() -> size();
	l = r*c;

	// intialise iterator and get dimensions of matrix
	auto iter = init_list.begin();

	// allocate array
	mp = new double[l];

	// populate array
	for(int i = 0; i < r; ++i, ++iter) {

		for(int j = 0; j < c; ++j) {

			mp[i*c + j] = (*iter)[j];
		}
	}
}

Matrix::Matrix(const Matrix& M) : 
	r(M.r), 
	c(M.c), 
	l(M.l) {

	mp = new double[M.l];

	memcpy(mp, M.mp, sizeof(double) * l);
}

Matrix::~Matrix() {

	delete[] mp;
}

std::vector<int> Matrix::dims() {

	return std::vector<int>(r, c);
}

Matrix Matrix::T() {

	Matrix temp(c, r);

	int i=0;
	int j=0;

	for(int t=0; t < l; ++t) {

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

double Matrix::sum() {

	double temp = 0;

	for(int t=0; t<l; ++t) {

		temp += mp[t];
	}

	return temp;
}

void Matrix::randomise(const double& a, const double& b) {

	for(int i = 0; i < l; ++i) {

		mp[i] = (double)std::rand() / RAND_MAX + a * (b-a);
	}
}

Matrix& Matrix::operator=(const Matrix& other) {

	r = other.r;
	c = other.c;
	l = r*c;

	delete[] mp;

	mp = new double[l];

	memcpy(mp, other.mp, sizeof(double) * l);

	return *this;
}


Matrix Matrix::operator+(const Matrix& other) {

	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		temp.mp[i] = mp[i] + other.mp[i];
	}

	return temp;
}

Matrix Matrix::operator-(const Matrix& other) {

	Matrix temp(r, c);

	for(int i=0; i<l; ++i) {

		temp.mp[i] = mp[i] - other.mp[i];
	}

	return temp;
}


Matrix Matrix::operator*(const Matrix& other) {

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


Matrix Matrix::for_each(double func(double&)) {

	Matrix temp(r, c);

	for(int t=0; t<l; ++t) {

		temp.mp[t] = func(mp[t]);
	}

	return temp;
}

double& Matrix::operator()(const int& i, const int& j) {

	double& temp = mp[i*c + j];

	return temp;
}

Matrix Matrix::slice(const int& r1, const int& c1, const int& r2, const int& c2) {

	if(r2 < r1 || c2 < c1) {

		throw std::invalid_argument("index of bottom left corner (2nd 2 args) must be greater than index of first corner (first 2 args)");
	}

	if(r1 > r || c1 > c) {

		throw std::invalid_argument("indices must all be inside matrix");
	}

	Matrix temp(r2-r1, c2-c1);

	int i = 0;
	int j = 0;

	for(int t=0; t < temp.l; ++t) {

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


std::ostream& operator<<(std::ostream& os, const Matrix& M) {

	os << "[[";

	for(int i=0; i < M.l-1; ++i) {

		os << M.mp[i];

		os << ( (i+1) % M.c != 0 ? ", " : "],\n [");
	}

	os << M.mp[M.l-1];

	os << "]]";

	return os;
}

Matrix array_mult(const Matrix& M1, const Matrix& M2) {

	Matrix temp(M1.r, M1.c);

	for(int t=0; t < M1.l; ++t) {

		temp.mp[t] = M1.mp[t] * M2.mp[t];
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