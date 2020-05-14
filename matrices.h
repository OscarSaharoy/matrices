#include <iostream>
#include <vector>

class Matrix {

public:

	Matrix(); // constructs a blank matrix
	Matrix(const int& r, const int& c); // constructs a matrix with r rows and c columns
	Matrix(std::initializer_list<std::initializer_list<double>> init_list); // constructs a matrix from nested initialiser lists
	Matrix(const Matrix& M); // copy constructor

	~Matrix(); // destructor

	// todo: make eigenvalue functions better

	std::vector<int> dims(); // returns vector of dimensions of matrix (rows, columns)
	Matrix T(); // returns transpose of matrix
	double det(); // returns the determinant of square matrix
	Matrix inv(); // returns inverse of matrix
	Matrix normalised(); // returns vector of length 1
	Matrix for_each(double func(double&)); // applies a function to each element in matrix
	Matrix squared(); // squares each element in matrix
	Matrix sqrt(); // square root of each element in matrix
	double sum(); // sums matrix
	double magnitude(); // returns magnitude of vector
	void randomise(const double& mean, const double& range); /// randomises values in matrix with uniform distribution

	Matrix& operator=(const Matrix& other); // assingment operator

	Matrix operator+(const Matrix& other); // adds 2 matrices element wise
	Matrix operator-(const Matrix& other); // subtracts 2 matrices element wise
	Matrix operator*(const Matrix& other); // matrix multiplication

	Matrix operator+(const double& other); // adds a double to each element in matrix
	Matrix operator-(const double& other); // subtracts a double from each element in matrix
	Matrix operator*(const double& other); // multiplies each element by a double
	Matrix operator/(const double& other); // divides each element by a double

	Matrix operator&&(const Matrix& other); // logical and - returns matrix of 1 and 0
	Matrix operator||(const Matrix& other); // logical or - returns matrix of 1 and 0
	Matrix operator!(); // logical not - returns matrix of 1 and 0

	void operator+=(const Matrix& other);
	void operator-=(const Matrix& other);
	void operator*=(const Matrix& other);

	void operator+=(const double& other);
	void operator-=(const double& other);
	void operator*=(const double& other);
	void operator/=(const double& other);

	double& operator()(const int& i, const int& j); // get element in matrix by index eg M(0,4)
	Matrix slice(const int& r1, const int& c1, const int& r2, const int& c2) const; // returns a new matrix that is a rectangle cut out of the old one
	Matrix minor(const int& r1, const int& c1) const; // returns the minor of the matrix from r1,c1
	Matrix minor(const int& t1) const; // returns the minor of the matrix from tM

	friend std::ostream& operator<<(std::ostream& os, const Matrix& M); // prints matrix to console
	friend Matrix operator*(const double& d, const Matrix& M); // multiplies each element by a double
	friend Matrix array_mult(const Matrix& M1, const Matrix& M2); // multiplies each element of a matrix by the corresponding element of another matrix
	friend Matrix array_div(const Matrix& M1, const Matrix& M2); // divides each element of a matrix by the corresponding element of another matrix
	friend Matrix join(const Matrix& M1, const Matrix& M2, const int& axis); // joins 2 matrices along a given axis
	friend double det(const Matrix& M1); // returns the determinant of a square matrix
	friend Matrix eigenvalues(const Matrix& M); // returns diagonal matrix of eigenvalues of M
	friend Matrix eigenvectors(Matrix& M); // returns matrix of eigenvectors of M

	friend double dot(const Matrix& V1, const Matrix& V2); // vector dot product
	friend Matrix cross(const Matrix& V1, const Matrix& V2); // vector cross product

	double* mp; // pointer to matrix array
	int r;      // number of rows
	int c;      // number of columns
	int l;      // length of array r*c
};