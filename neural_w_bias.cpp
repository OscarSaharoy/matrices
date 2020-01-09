#include "matrices.cpp"
#include <cmath>
#include <vector>

// Leaky ReLU

double act(double& x) {

	return (x > 0) ? x : 0.01*x;
}

double der(double& x) {

	return (x > 0) ? 1 : 0.01;
}


// tanh

// double act(double& x) {

// 	return std::tanh(x);
// }

// double der(double& x) {

// 	return 1 / (std::cosh(x) * std::cosh(x));
// }

double clip(const double& x, const double& l) {

	return std::abs(x) < l ? x : l * x/std::abs(x);
}

// adds bias

Matrix add_bias(Matrix& M) {

	Matrix temp(M.r + 1, M.c);

	memcpy(temp.mp + 1, M.mp, M.l * sizeof(double));

	temp.mp[0] = 1;

	return temp;
}

Matrix act(Matrix& M) {

	Matrix temp = M.for_each(act);

	temp.mp[0] = 1;

	return temp;
}

Matrix der(Matrix& M) {

	Matrix temp = M.for_each(der);

	temp.mp[0] = 0;

	return temp;
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

	Matrix Din({{3},
		        {2}});

	Matrix Dout({{0}});

	std::vector<Matrix> training = { Ain,  Bin,  Cin,  Din};
	std::vector<Matrix> targets  = {Aout, Bout, Cout, Dout};

	Matrix W1(5, 3);
	W1.randomise(0, 0.05);

	Matrix W2(1, 5);
	W2.randomise(0, 0.05);

	double lr = 0.01;
	int epoch = training.size();

	for(int i = 0; i < epoch*10000; ++i) {

		Matrix I1    = add_bias(training[i%4]);
		Matrix T1    = targets[i%4];

		Matrix H1a   = W1 * I1;
		Matrix H1b   = act(H1a);

		Matrix O1    = W2 * H1b;

		Matrix dW2   = (O1 - T1) * H1b.T();
		Matrix dW1   = array_mult( W2.T() * (O1 - T1), der(H1b) ) * I1.T();

		W2 = W2 - lr * dW2;
		W1 = W1 - lr * dW1;

		double e     = 0.5 * std::pow( (O1(0,0) - T1(0,0)), 2 );

		if(i % 1000 == 0) {

			std::cout << e << std::endl;
		}
	}

	std::cout << "Training Results\n" << std::endl;

	for(int i = 0; i < epoch; ++i) {

		Matrix I1 = add_bias(training[i]);
		Matrix T1 = targets[i];

		Matrix H1a = W1 * I1;
		Matrix H1b = act(H1a);

		Matrix O1  = W2 * H1b;

		std::cout << "Input: \n"  << training[i] << std::endl;
		std::cout << "Output: \n" << O1 << "\n" << std::endl;
	}

	std::cout << "Weights:\n" << std::endl;
	std::cout << W1 << "\n\n" << W2 << "\n" << std::endl;
}
