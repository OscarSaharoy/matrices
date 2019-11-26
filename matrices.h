
#include <vector>
#include <initializer_list>
#include <iostream>

class Matrix: public std::vector<std::vector<double>> {
public:

	Matrix();
	Matrix(std::vector<std::vector<double>>);
	Matrix(std::initializer_list<std::vector<double>> init_list);

	void display();
	friend std::ostream& operator<<(std::ostream& os, const Matrix& m);
	friend Matrix array_mult(const Matrix& m1, const Matrix& m2);

	Matrix T();
	double det();
	std::vector<int> dims() const;

	Matrix operator+(const Matrix& other);
	Matrix operator-(const Matrix& other);
	Matrix operator*(const Matrix& other);

	Matrix operator+(const double& other);
	Matrix operator-(const double& other);
	Matrix operator*(const double& other);
	Matrix operator/(const double& other);

	Matrix operator&&(const Matrix& other);
	Matrix operator||(const Matrix& other);
	Matrix operator!();

	double operator()(const int& i1, const int& i2);
};