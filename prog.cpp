#include "mpi.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <functional>

using namespace std;


int world_rank;
int world_size;

int lattice_n;
int lattice_m;

int ComputeLatticeId(int i, int j) {
	if (i > lattice_n) {
		throw "DFSDFSDF";
	}

	if (j > lattice_m) {
		throw "AAAAAAA";
	}

	return lattice_m * i + j;
}

void ComputeLatticeCoord(int id, int* i, int* j) {
	*j = id % lattice_m;
	*i = id / lattice_n;
}

void GetLatticeParams() {
	lattice_n = 1;
	lattice_m = 1;
}

void ComputeBlockSize(int N, int M, int* block_h, int* block_w) {
	*block_h = N;
	*block_w = M;
}

using V = vector<double>;
void PrintMatrix(const V& x, int N, int M) {
	if (x.size() != N * M) {
		throw "asdfasdf";
	}

	for (int i = 0; i < N; ++i) {
		for(int j = 0; j < M; ++j) {
			cout << x[j + i*M] << " ";
		}
		cout << endl;
	}
}

V apply(function<double(double, double)> foo, const V& x, const V& y) {
	if (x.size() != y.size()) {
		throw "asdfasd";
	}
	
	V result(x.size());
	for (size_t i = 0; i < x.size(); ++i) {
		result[i] = foo(x[i], y[i]);
	}
	
	return result;
}


double sqr(double x) {
	return x*x;
}

double func_k(double x, double y) {
	return x + 4;
}

double func_q(double x, double y) {

	auto tmp = x + y;
	return (tmp >= 0) * sqr(x + y);
}

double func_p(double x, double y) {
	return -2 * (x + y);
}

double func_u(double x, double y) {
	return exp(1 - sqr(x + y));
}

auto func_phi = func_u;

double func_F(double x, double y) {
	auto inner = (-2 * (sqr(func_p(x, y)) - 2) * func_k(x, y) - func_p(x, y)  + func_q(x, y));

	return func_u(x, y) * inner;
}


void MeshGrid(V& x, V& y, double x_min, double x_max, double y_min, double y_max, int N, int M, double* h1, double* h2) {
	*h1 = (x_max - x_min)/(N - 1);
	*h2 = (y_max - y_min)/(M - 1);

	x.resize(N * M);
	y.resize(N * M);

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < M; ++j) {
			x[j + M * i] = x_min + i * *h1;
			y[i + N * j] = y_min + j * *h2;
		}
	}
}



void MeshGrid2(V& x, V& y, double x_min, double y_min, int N, int M, double h1, double h2) {

	x.resize(N * M);
	y.resize(N * M);

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < M; ++j) {
			x[j + M * i] = x_min + i * h1;
		}
	}

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < M; ++j) {
			y[j + M * i] = y_min + j * h2;
		}
	}



}

V operator + (const V& a, const V& b) {
	if (a.size() != b.size()) {
		throw "SDFSDFS";
	}

	V result(a.size());
	for (size_t i = 0; i < a.size(); ++i) {
		result[i] = a[i] + b[i];
	}

	return result;
}

V operator * (const V& a, double num) {

	V result(a.size());
	for (size_t i = 0; i < a.size(); ++i) {
		result[i] = a[i] * num;
	}

	return result;
}

V operator / (const V& a, double num) {

	V result(a.size());
	for (size_t i = 0; i < a.size(); ++i) {
		result[i] = a[i] / num;
	}

	return result;
}


class LocalOperator {
public:
	const int N, M;
protected:
	int row_pos, col_pos;
	int block_h, block_w;
	double h1, h2;
	double my_x_min, my_y_min;
	V x_, y_;
public:
	LocalOperator(int my_i, int my_j, int block_h_, int block_w_, int N_, int M_, double x_min, double x_max, double y_min, double y_max):
		N(N_), M(M_)
	 {
		row_pos = my_i;
		col_pos = my_j;
		block_h = block_h_;
		block_w = block_w_;

		h1 = (x_max - x_min)/(N - 1);
		h2 = (y_max - y_min)/(M - 1);
		
		my_x_min = x_min + row_pos * h1;
		my_y_min = y_min + col_pos * h2;
		
		MeshGrid2(x_, y_, my_x_min, my_y_min, block_h, block_w, h1, h2);
		


		BuildLhs();
	 }

	void BuildLhs() {
		V w_ij_coefs(x_.size());
		V w_ip1j_coefs(x_.size());
		V w_im1j_coefs(x_.size());
		V w_ijp1_coefs(x_.size());
		V w_ijm1_coefs(x_.size());

		cout << h1 << " " << h2 << endl;
		cout << "=========\n";
		

		for (int i = 0; i < block_h; ++i) {
			for (int j = 0; j < block_w; ++j) {
				int pos = j + i * M;
				w_ij_coefs[pos] = 
(func_k(x_[pos] + 0.5 * h1, y_[pos]) +  func_k(x_[pos] - 0.5 * h1, y_[pos]))/sqr(h1) + 
(func_k(x_[pos], y_[pos] + 0.5 * h2) +  func_k(x_[pos], y_[pos] - 0.5 * h2))/sqr(h2) + func_q(x_[pos], y_[pos]);
		
				w_ip1j_coefs[pos] = -func_k(x_[pos] + 0.5 * h1, y_[pos])/sqr(h1);
				w_im1j_coefs[pos] = -func_k(x_[pos] - 0.5 * h1, y_[pos])/sqr(h1);
				w_ijp1_coefs[pos] = -func_k(x_[pos], y_[pos] + 0.5 * h2)/sqr(h2);
				w_ijm1_coefs[pos] = -func_k(x_[pos], y_[pos] - 0.5 * h2)/sqr(h2);
			}
		} 
		PrintMatrix(w_ij_coefs, block_h, block_w);


 
	}
};

void SendBorder(const vector<double>& border, int dst, int me) {
	MPI_Send(border.data(), border.size(), MPI_DOUBLE, dst, me, MPI_COMM_WORLD);
}

void RecvBorder(vector<double>& border, int src) {
	MPI_Status status;
	MPI_Recv(border.data(), border.size(), MPI_DOUBLE, src, src, MPI_COMM_WORLD, &status);
}

vector<double> SyncBorder(int my_i, int my_j, int other_i, int other_j, vector<double> my_border) {
	vector<double> other_border(my_border.size());

	
	SendBorder(my_border, ComputeLatticeId(other_i, other_j), ComputeLatticeId(my_i, my_j));
	RecvBorder(other_border, ComputeLatticeId(other_i, other_j));

	return other_border;
}
/*
class Block {
	vector<double> data;
	vector<double> border1;
	...
	int block_row, int block_col, int width, int height;

	operator +;
	operator -;
	operator *;
}


Block Step(const Block& w, int block_row, int block_col, int width, int height) {
	auto Aw = OperatorA(w, w.block_row, w.block_col, w.width, w.height);
	auto rhs = RHS(w.block_row, w.block_col, w.width, w.height);
	residual = Aw - rhs;
	
	auto Ar = Operator(residual, w.block_row, w.block_col, ...);
	
	double tau = Dot(r, Ar)/Dot(Ar, Ar);

	auto new_w = w - tau*r;
	
	new_w.setOtherUpperBorder((SyncBorder(my_i, my_j, my_i + 1, my_j, new_w.getMyUpperBorder()));
	


	return new_w;
}

*/

void ParseArgs(int argc, char **argv, int* N, int* M) {
	*N = 4;
	*M = 7;
}

const double X_MIN = -1, X_MAX = 2;
const double Y_MIN = -2, Y_MAX = 2;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int N, M;
	ParseArgs(argc, argv, &N, &M);

	GetLatticeParams();

	int my_i, my_j;
	int block_h, block_w;
	
	ComputeLatticeCoord(world_rank, &my_i, &my_j);
	ComputeBlockSize(N, M, &block_h, &block_w);
	
	LocalOperator op(my_i, my_j, N, M, block_h, block_w, X_MIN, X_MAX, Y_MIN, Y_MAX);
	
	/*
	auto my_w = Solve();
	SendLocalResultsTo0(my_w);
	
	if(world_rank == 0) {
		auto result = ReceiveLocalResults();
		DoIOStuff();
	}
	*/

	MPI_Finalize();
}
