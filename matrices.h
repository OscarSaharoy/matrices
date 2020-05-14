#include <iostream>
#include <vector>

class Matrix {

public:

	Matrix();                                                               // constructs a blank matrix
	Matrix(const int& r, const int& c);                                     // constructs a matrix with r rows and c columns
	Matrix(std::initializer_list<std::initializer_list<double>> init_list); // constructs a matrix from nested initialiser lists
	Matrix(const Matrix& M);                                                // copy constructor

	~Matrix(); // destructor

	Matrix dims() const;       // returns matrix of dimensions of matrix (rows, columns)

	Matrix T() const;          // returns transpose of matrix
	Matrix inv() const;        // returns inverse of matrix

	double det() const;        // returns the determinant of square matrix
	double sum() const;        // sums matrix
	double magnitude() const;  // returns magnitude of vector

	Matrix normalised() const; // returns vector of length 1
	Matrix squared() const;    // squares each element in matrix
	Matrix sqrt() const;       // square root of each element in matrix

	void for_each(double func(double&));                     // applies a function to each element in matrix
	void identity();                                         // sets diagonal to 1, other elements to 0
	void randomise(const double& mean, const double& range); // randomises values in matrix with uniform distribution

	Matrix& operator=(const Matrix& other);       // assingment operator

	Matrix operator+(const Matrix& other) const;  // adds 2 matrices element wise
	Matrix operator-(const Matrix& other) const;  // subtracts 2 matrices element wise
	Matrix operator*(const Matrix& other) const;  // matrix multiplication

	Matrix operator+(const double& other) const;  // adds a double to each element in matrix
	Matrix operator-(const double& other) const;  // subtracts a double from each element in matrix
	Matrix operator*(const double& other) const;  // multiplies each element by a double
	Matrix operator/(const double& other) const;  // divides each element by a double

	Matrix operator&&(const Matrix& other) const; // logical and - returns matrix of 1 and 0
	Matrix operator||(const Matrix& other) const; // logical or - returns matrix of 1 and 0
	Matrix operator!() const;                     // logical not - returns matrix of 1 and 0

	void operator+=(const Matrix& other);
	void operator-=(const Matrix& other);
	void operator*=(const Matrix& other); // this = other * this

	void operator+=(const double& other);
	void operator-=(const double& other);
	void operator*=(const double& other);
	void operator/=(const double& other);

	Matrix minor(const int& r1, const int& c1) const; // returns the minor of the matrix from r1,c1
	Matrix minor(const int& t1) const;                // returns the minor of the matrix from tM
	double& operator()(const int& i, const int& j);   // get element in matrix by index eg M(0,4)
	double& operator[](const int& t);                 // get element in matrix as if the matrix were a normal array
	Matrix slice(const int& r1, const int& c1, const int& r2, const int& c2) const; // returns a new matrix that is a rectangle cut out of the old one

	friend std::ostream& operator<<(std::ostream& os, const Matrix& M);      // prints matrix to console
	friend Matrix operator*(const double& d, const Matrix& M);               // multiplies each element by a double
	friend Matrix array_mult(const Matrix& M1, const Matrix& M2);            // multiplies each element of a matrix by the corresponding element of another matrix
	friend Matrix array_div(const Matrix& M1, const Matrix& M2);             // divides each element of a matrix by the corresponding element of another matrix
	friend Matrix join(const Matrix& M1, const Matrix& M2, const int& axis); // joins 2 matrices along a given axis
	friend double det(const Matrix& M1);                                     // returns the determinant of a square matrix

	friend Matrix eigenvectors(const Matrix& M);                    // returns matrix of eigenvectors of M
	friend Matrix eigenvalues(const Matrix& M, const Matrix& eigM); // returns diagonal matrix of eigenvalues from a matrix of eigenvectors

	friend double dot(const Matrix& V1, const Matrix& V2);   // vector dot product
	friend Matrix cross(const Matrix& V1, const Matrix& V2); // vector cross product

	double* mp; // pointer to matrix array
	int r;      // number of rows
	int c;      // number of columns
	int l;      // length of array r*c
};