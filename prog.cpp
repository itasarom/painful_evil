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

struct Matrix {
	V data;
	int M, N;
	Matrix(const V& data_, int N_, int M_):data(data_), N(N_), M(M_) {
		if (data.size() != N * M) {
			throw "asdfasdf";
		}
	}
	
	Matrix(int N_, int M_, double fill=0.0): data(N_ * M_, fill), N(N_), M(M_) {
	}
	
	Matrix() = default;
	double& operator[] (int i) {
		return data[i];
	}

	double& at(int i, int j) {
		int pos = j + i * M; 
		return data[pos];
	}

	void reset(int N_, int M_, double fill=0.0) {
		N = N_;
		M = M_;
		data.clear();
		data.resize(N_ * M_, fill);
	}
};


void PrintMatrix(Matrix& m) {

	for (int i = 0; i < m.N; ++i) {
		for(int j = 0; j < m.M; ++j) {
			cout << m.at(i, j) << " ";
		}
		cout << endl;
	}
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




void MeshGrid2(Matrix& x, Matrix& y, double x_min, double y_min, int N, int M, double h1, double h2) {

	x.reset(N, M);
	y.reset(N, M);

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < M; ++j) {
			x.at(i, j) = x_min + i * h1;
		}
	}

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < M; ++j) {
			y.at(i, j) = y_min + j * h2;
		}
	}



}



class LocalOperator {
public:
	const int N, M;
protected:
	int row_pos, col_pos;
	int block_h, block_w;
	double h1, h2;
	double my_x_min, my_y_min;
	Matrix x_, y_;
public:
	LocalOperator(int my_i, int my_j, int block_h_, int block_w_, int N_, int M_, double x_min, double x_max, double y_min, double y_max):
		N(N_), M(M_)
	 {
		row_pos = my_i;
		col_pos = my_j;
		block_h = block_h_;
		block_w = block_w_;

		h1 = (x_max - x_min)/(block_h - 1);
		h2 = (y_max - y_min)/(block_w - 1);
		
		my_x_min = x_min + row_pos * h1;
		my_y_min = y_min + col_pos * h2;
		
		MeshGrid2(x_, y_, my_x_min, my_y_min, block_h, block_w, h1, h2);
		


		BuildLhs();
		BuildRhs();
	 }

	Matrix w_ij_coefs, w_ip1j_coefs, w_im1j_coefs, w_ijp1_coefs, w_ijm1_coefs, rhs, phi;

	void BuildLhs() {
		w_ij_coefs.reset(block_h, block_w);
		w_ip1j_coefs.reset(block_h, block_w);
		w_im1j_coefs.reset(block_h, block_w);
		w_ijp1_coefs.reset(block_h, block_w);
		w_ijm1_coefs.reset(block_h, block_w);

		cout << h1 << " " << h2 << endl;
		cout << "=========\n";
		

		for (int i = 0; i < block_h; ++i) {
			for (int j = 0; j < block_w; ++j) {
				w_ij_coefs.at(i, j) = 
(func_k(x_.at(i, j) + 0.5 * h1, y_.at(i, j)) +  func_k(x_.at(i, j) - 0.5 * h1, y_.at(i, j)))/sqr(h1) + 
(func_k(x_.at(i, j), y_.at(i, j) + 0.5 * h2) +  func_k(x_.at(i, j), y_.at(i, j) - 0.5 * h2))/sqr(h2) + func_q(x_.at(i, j), y_.at(i, j));
		
				w_ip1j_coefs.at(i, j) = -func_k(x_.at(i, j) + 0.5 * h1, y_.at(i, j))/sqr(h1);
				w_im1j_coefs.at(i, j) = -func_k(x_.at(i, j) - 0.5 * h1, y_.at(i, j))/sqr(h1);
				w_ijp1_coefs.at(i, j) = -func_k(x_.at(i, j), y_.at(i, j) + 0.5*h2)/sqr(h2);
				w_ijm1_coefs.at(i, j) = -func_k(x_.at(i, j), y_.at(i, j) - 0.5*h2)/sqr(h2);
		
			}
		} 

 
	}
	
		

	void BuildRhs() {
		rhs.reset(block_h, block_w );
		
		for (int i = 0; i < block_h; ++i) {
			for (int j = 0; j < block_w; ++j) {
				rhs.at(i, j) = func_F(x_.at(i, j), y_.at(i,j));
			}
		}

		if (row_pos == 0) {
			for (int j = 0; j < block_w; ++j) {
				int i = 0;
				rhs.at(i, j) -= func_phi(x_.at(i, j) - h1,  y_.at(i, j)) * w_im1j_coefs.at(i, j);
			}
		}

		if (row_pos + block_h == N) {
			for (int j = 0; j < block_w; ++j) {
				int i = block_h - 1;
				rhs.at(i, j) -= func_phi(x_.at(i, j) + h1,  y_.at(i, j)) * w_ip1j_coefs.at(i, j);
			}
			
		}

		if (col_pos == 0) {
			for (int i = 0; i < block_h; ++i) {
				int j = 0;
				rhs.at(i, j) -= func_phi(x_.at(i, j),  y_.at(i, j) - h2) * w_ijm1_coefs.at(i, j);
			}
		}

		if (col_pos + block_w == M) {
			for (int i = 0; i < block_h; ++i) {
				int j = block_w - 1;
				rhs.at(i, j) -= func_phi(x_.at(i, j),  y_.at(i, j) + h2) * w_ijp1_coefs.at(i, j);
			}
		}


		PrintMatrix(rhs);
		
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
