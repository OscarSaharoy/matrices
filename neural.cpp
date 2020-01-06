#include "matrices.cpp"
#include <cmath>
#include <vector>

double act(double& x) {

	return (x > 0) ? x : 0.01*x;
}

double der(double& x) {

	return (x > 0) ? 1 : 0.01;
}

double erf(Matrix& X, Matrix& Y) {

	double sum = 0;

	for(int t=0; t < X.l; ++t) {

		sum += pow(X.mp[t] - Y.mp[t], 2);
	}

	return sum / (2 * X.l);
}

int main() {

	Matrix Ain({{0},
		        {0}});

	Matrix Aout({{0}});

	Matrix Bin({{0},
		        {1}});

	Matrix Bout({{1}});

	Matrix Cin({{1},
		        {0}});

	Matrix Cout({{1}});

	Matrix Din({{1},
		        {1}});

	Matrix Dout({{0}});

	std::vector<Matrix> training = { Ain,  Bin,  Cin,  Din};
	std::vector<Matrix> targets  = {Aout, Bout, Cout, Dout};

	Matrix W1(3, 2);
	W1.randomise(0, 0.1);

	Matrix W2(1, 3);
	W2.randomise(0, 0.1);

	double lr  = 0.1;
	int setlen = training.size();

	for(int i = 0; i < setlen*1000; ++i) {

		Matrix I1 = training[i%4];
		Matrix T1 = targets[i%4];

		Matrix H1a = W1 * I1;
		Matrix H1b = H1a.for_each(act);

		Matrix O1  = W2 * H1b;

		double error = std::pow( (O1 - T1).sum(), 2 );

		Matrix dW2 = (O1 - T1) * H1b.T();
		Matrix dW1 = array_mult( (W2.T() * (O1 - T1)), H1a.for_each(der) ) * I1.T();

		W2 = W2 - dW2 * lr;
		W1 = W1 - dW1 * lr;
	}

	std::cout << "Training Results\n" << std::endl;

	for(int i = 0; i < setlen; ++i) {

		Matrix I1 = training[i];
		Matrix T1 = targets[i];

		Matrix H1a = W1 * I1;
		Matrix H1b = H1a.for_each(act);

		Matrix O1  = W2 * H1b;

		std::cout << "Input: \n" << I1 << std::endl;
		std::cout << "Output : \n" << O1 << std::endl << std::endl;
	}
}
