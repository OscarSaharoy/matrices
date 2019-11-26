#include "matrices.h"

Matrix::Matrix() {};

Matrix::Matrix(std::vector<std::vector<double>> vec) {

	for(auto iter = vec.begin(); iter != vec.end(); ++iter) {

		this->push_back(*iter);
	}
}

Matrix::Matrix(std::initializer_list<std::vector<double>> init_list) {

	for(auto iter = init_list.begin(); iter != init_list.end(); ++iter) {

		this->push_back(*iter);
	}
}

void Matrix::display() {

	std::cout << "[";

	for(auto iter1 = begin(); iter1 != end(); ++iter1) {

		std::cout << "[";

		for(auto iter2 = iter1->begin(); iter2 != iter1->end(); ++iter2) {

			std::cout << *iter2 << ", ";
		}

		std::cout << "]";

		std::cout << std::endl;
	}

	std::cout << "]";
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {

	std::vector<int> dims_vector = m.dims();

	os << "[";

	for(int i=0; i < dims_vector[0]; ++i) {

		os << "[";

		for(int j=0; j < dims_vector[1]; ++j) {

			os << m[i][j];

			os << (j != dims_vector[1]-1 ? ", " : "]");

		}

		os << (i != dims_vector[0]-1 ? ",\n " : "]\n");
	}

	return os;
}

Matrix array_mult(const Matrix& m1, const Matrix& m2) {

	Matrix temp;

	for(int i=0; i < m1.size(); i++) {

		std::vector<double> emp;
		temp.push_back(emp);

		for(int j=0; j < m1[i].size(); j++) {

			temp[i].push_back( m1[i][j] * m2[i][j] );
		}
	}

	return temp;

}

Matrix Matrix::T() {

	Matrix temp;
	std::vector<int> dims_vector = dims();

	for(int i=0; i < dims_vector[1]; i++) {

		std::vector<double> row;

		for(int j=0; j < dims_vector[0]; j++) {

			row.push_back( (*this)[j][i] );
		}

		temp.push_back(row);
	}

	return temp;
}

double Matrix::det() {

	Matrix &t = *this;

	double d1 = t[0][0] * ( t[1][1]*t[2][2] - t[1][2]*t[2][1]);
	double d2 = t[0][1] * ( t[1][0]*t[2][2] - t[1][2]*t[2][0]);
	double d3 = t[0][2] * ( t[1][0]*t[2][1] - t[1][1]*t[2][0]);

	return d1 - d2 + d3;
}

std::vector<int> Matrix::dims() const {

	std::vector<int> temp;

	temp.push_back( size() );
	temp.push_back( (*this)[0].size() );

	return temp;
}

Matrix Matrix::operator+(const Matrix& other) {

	Matrix temp;

	for(int i=0; i < size(); i++) {

		std::vector<double> emp;
		temp.push_back(emp);

		for(int j=0; j < (*this)[i].size(); j++) {

			temp[i].push_back( (*this)[i][j] + other[i][j] );
		}
	}

	return temp;
}

Matrix Matrix::operator-(const Matrix& other) {

	Matrix temp;

	for(int i=0; i < size(); i++) {

		std::vector<double> emp;
		temp.push_back(emp);

		for(int j=0; j < (*this)[i].size(); j++) {

			temp[i].push_back( (*this)[i][j] - other[i][j] );
		}
	}

	return temp;
}

Matrix Matrix::operator*(const Matrix& other) {

	Matrix temp;

	for(int i=0; i < size(); i++) {

		std::vector<double> in_temp;

		for(int j=0; j < other[0].size(); j++) {

			double sum = 0;

			for(int k=0; k < other.size(); k++) {

				sum += (*this)[i][k] * other[k][j];
			}

			in_temp.push_back(sum);
		}

		temp.push_back(in_temp);
	}

	return temp;
}

Matrix Matrix::operator+(const double& other) {

	Matrix temp;

	for(int i=0; i < size(); i++) {

		std::vector<double> emp;
		temp.push_back(emp);

		for(int j=0; j < (*this)[i].size(); j++) {

			temp[i].push_back( (*this)[i][j] + other );
		}
	}

	return temp;
}

Matrix Matrix::operator-(const double& other) {

	Matrix temp;

	for(int i=0; i < size(); i++) {

		std::vector<double> emp;
		temp.push_back(emp);

		for(int j=0; j < (*this)[i].size(); j++) {

			temp[i].push_back( (*this)[i][j] - other );
		}
	}

	return temp;
}

Matrix Matrix::operator*(const double& other) {

	Matrix temp;

	for(int i=0; i < size(); i++) {

		std::vector<double> emp;
		temp.push_back(emp);

		for(int j=0; j < (*this)[i].size(); j++) {

			temp[i].push_back( (*this)[i][j] * other );
		}
	}

	return temp;
}

Matrix Matrix::operator/(const double& other) {

	Matrix temp;

	for(int i=0; i < size(); i++) {

		std::vector<double> emp;
		temp.push_back(emp);

		for(int j=0; j < (*this)[i].size(); j++) {

			temp[i].push_back( (*this)[i][j] / other );
		}
	}

	return temp;
}

Matrix Matrix::operator&&(const Matrix& other) {

	Matrix temp;

	for(int i=0; i < size(); i++) {

		std::vector<double> emp;
		temp.push_back(emp);

		for(int j=0; j < (*this)[i].size(); j++) {

			temp[i].push_back( (*this)[i][j] && other[i][j] );
		}
	}

	return temp;
}

Matrix Matrix::operator||(const Matrix& other) {

	Matrix temp;

	for(int i=0; i < size(); i++) {

		std::vector<double> emp;
		temp.push_back(emp);

		for(int j=0; j < (*this)[i].size(); j++) {

			temp[i].push_back((*this)[i][j] || other[i][j]);
		}
	}

	return temp;
}

Matrix Matrix::operator!() {

	Matrix temp;

	for(int i=0; i < size(); i++) {

		std::vector<double> emp;
		temp.push_back(emp);

		for(int j=0; j < (*this)[i].size(); j++) {

			temp[i].push_back( !(*this)[i][j] );
		}
	}

	return temp;
}

double Matrix::operator()(const int& i1, const int& i2) {

	return (*this)[i1][i2];
}