#include <iostream>
#include <vector>

class Matrix {

public:

	Matrix();
	Matrix(const int& r, const int& c);
	Matrix(std::initializer_list<std::vector<double>> init_list);
	Matrix(const Matrix& M);

	~Matrix();

	// todo: det, inv

	std::vector<int> dims();
	Matrix T();
	double sum();
	void randomise(const double& mean, const double& range);

	Matrix& operator=(const Matrix& other);

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

	Matrix for_each(double func(double&));

	double& operator()(const int& i, const int& j);
	Matrix slice(const int& r1, const int& c1, const int& r2, const int& c2);

	friend std::ostream& operator<<(std::ostream& os, const Matrix& M);
	friend Matrix array_mult(const Matrix& M1, const Matrix& M2);
	friend Matrix join(const Matrix& M1, const Matrix& M2, const int& axis);

	double* mp; // pointer to matrix array
	int r;      // number of rows
	int c;      // number of columns
	int l;      // length of array r*c
};