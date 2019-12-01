#include "mpi.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <functional>
#include <sstream>
#include <string>
#include <cstdlib>

using namespace std;


stringstream out;
int world_rank;
int world_size;

int lattice_n;
int lattice_m;
int my_i, my_j;

int ComputeLatticeId(int i, int j) {
	if (i < 0 || i >= lattice_n) {
		throw "DFSDFSDF";
	}

	if (j < 0 || j >= lattice_m) {
		throw "AAAAAAA";
	}

	return lattice_m * i + j;
}

void ComputeLatticeCoord(int id, int* i, int* j) {
	*j = id % lattice_m;
	*i = id / lattice_m;
}

void GetLatticeParams() {
	lattice_n = 1;
	for (int z = 1; z <= sqrt(world_size); ++z) {
		if (world_size % z == 0) {
			lattice_n = max(lattice_n, z);
		}
	}
	
	lattice_m = world_size / lattice_n;
}

void ComputeBlockSize(int N, int M, int* block_h, int* block_w) {
	if (N % lattice_n == 0) {
		*block_h = N/lattice_n;
	} else {
		if (my_i == lattice_n - 1) {
			*block_h = N -  N/lattice_n * (lattice_n - 1);
		} else {
			*block_h = N / lattice_n;
		}	
	}

	if (M % lattice_m == 0) {
		*block_w = M/lattice_m;
	} else {
		if (my_j == lattice_m - 1) {
			*block_w = M -  M/lattice_m * (lattice_m - 1);
		} else {
			*block_w = M / lattice_m;
		}	
	}


}

typedef vector<double> V;

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
	
	Matrix() {
	}
	double& operator[] (int i) {
		return data[i];
	}

	double& at(int i, int j) {
		int pos = j + i * M; 
		return data[pos];
	}

	double cat(int i, int j) const {
		int pos = j + i * M; 
		return data[pos];
	}

	void reset(int N_, int M_, double fill=0.0) {
		N = N_;
		M = M_;
		data.clear();
		data.resize(N_ * M_, fill);
	}


	Matrix operator - (const Matrix& other) const {
		Matrix result(N, M);
		
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < M; ++j) {
				result.at(i, j) = cat(i, j) - other.cat(i, j);
			}
		}

		return result;
	}

	Matrix operator + (const Matrix& other) const {
		Matrix result(N, M);
		
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < M; ++j) {
				result.at(i, j) = cat(i, j) + other.cat(i, j);
			}
		}

		return result;
	}

	
	Matrix operator * (double num) const {
		Matrix result(N, M);
		
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < M; ++j) {
				result.at(i, j) = cat(i, j) * num;
			}
		}

		return result;
	}


	
};


void PrintMatrix(const Matrix& m) {

	stringstream out;
	for (int i = 0; i < m.N; ++i) {
		for(int j = 0; j < m.M; ++j) {
			out << m.cat(i, j) << " ";
		}
		out << endl;
	}
		
	cout << out.str();
}
double sqr(double x) {
	return x*x;
}

double func_k(double x, double y) {
	return x + 4;
}

double func_q(double x, double y) {

	double tmp = x + y;
	return (tmp >= 0) * sqr(x + y);
}

double func_p(double x, double y) {
	return -2 * (x + y);
}

double func_u(double x, double y) {
	return exp(1 - sqr(x + y));
}

double func_phi(double x, double y) {
	return func_u(x, y);
}

double func_F(double x, double y) {
	double inner = (-2 * (sqr(func_p(x, y)) - 2) * func_k(x, y) - func_p(x, y)  + func_q(x, y));

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
	int block_h, block_w;
protected:
	int row_pos, col_pos;
	double h1, h2;
	double my_x_min, my_y_min;
public:
	Matrix x_, y_;
	LocalOperator(int row_pos_, int col_pos_, int N_, int M_,int block_h_, int block_w_,  double x_min, double x_max, double y_min, double y_max):
		N(N_), M(M_)
	 {
		row_pos = row_pos_;
		col_pos = col_pos_;
		block_h = block_h_;
		block_w = block_w_;

		h1 = (x_max - x_min)/(N - 1);
		h2 = (y_max - y_min)/(M - 1);
		
		
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
		phi.reset(block_h, block_w);

		for (int i = 0; i < block_h; ++i) {
			for (int j = 0; j < block_w; ++j) {
				phi.at(i, j) = func_phi(x_.at(i, j),  y_.at(i, j));
			}
		}

		
		rhs.reset(block_h, block_w);
		
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


	
	}


	Matrix Call(const Matrix& w, const V& border_x0, const V& border_x1, const V& border_y0, const V& border_y1) {
		if (w.N != w_ij_coefs.N || w.M != w_ij_coefs.M) {
			throw "WTF";
		}	
		
		Matrix result(w.N, w.M);
		
		for (int i = 0; i < block_h; ++i) {
			for (int j = 0; j < block_w; ++j) {
				result.at(i, j) = w_ij_coefs.cat(i, j) * w.cat(i, j);
			}
		}
		
		
		//==================
		
		for (int i = 0; i < block_h - 1; ++i) {
			for (int j = 0; j < block_w; ++j) {
				result.at(i, j) += w_ip1j_coefs.cat(i, j) * w.cat(i + 1, j);
			}
		}
		if (border_x1.size() != block_w) {
			throw "asdasdfa";
		}
		
		for (int j = 0; j < block_w; ++j) {
			int i = block_h - 1;
			result.at(i, j) += w_ip1j_coefs.cat(i, j) * border_x1[j];
		}
		
		//==================

		for (int i = 1; i < block_h; ++i) {
			for (int j = 0; j < block_w; ++j) {
				result.at(i, j) += w_im1j_coefs.cat(i, j) * w.cat(i - 1, j);
			}
		}
		if (border_x0.size() != block_w) {
			throw "asdasdfa";
		}

		for (int j = 0; j < block_w; ++j) {
			int i = 0;
			result.at(i, j) += w_im1j_coefs.cat(i, j) * border_x0[j];
		}

		//==================
		for (int i = 0; i < block_h; ++i) {
			for (int j = 0; j < block_w - 1; ++j) {
				result.at(i, j) += w_ijp1_coefs.cat(i, j) * w.cat(i, j + 1);
			}
		}
		if (border_y1.size() != block_h) {
			throw "asdasdfa";
		}

		for (int i = 0; i < block_h; ++i) {
			int j = block_w - 1;
			result.at(i, j) += w_ijp1_coefs.cat(i, j) * border_y1[i];
		}

		//==================
		for (int i = 0; i < block_h; ++i) {
			for (int j = 1; j < block_w; ++j) {
				result.at(i, j) += w_ijm1_coefs.cat(i, j) * w.cat(i, j - 1);
			}
		}
		
		if (border_y0.size() != block_h) {
			throw "asdasdfa";
		}

		for (int i = 0; i < block_h; ++i) {
			int j = 0;
			result.at(i, j) += w_ijm1_coefs.cat(i, j) * border_y0[i];
		}
		

		return result;
	}
		
	double Dot(const Matrix &a, const Matrix& b) const {
		double result = 0.0;
		if (a.N != block_h || a.M != block_w) {
			throw "asdfasd";
		}
		if (b.N != block_h || b.M != block_w) {
			throw "asdfasd";
		}
		for (int i = 0; i < block_h; ++i) {
			for (int j = 0; j < block_w; ++j) {
				result += a.cat(i,j) * b.cat(i, j);
			}
		}

		
		return result * (h1 * h2);
	}

	double Norm2(const Matrix& a) const {
		return Dot(a, a);
	}

	double Norm(const Matrix& a) const {
		return sqrt(Norm2(a));
	}
	
};


void ParseArgs(int argc, char **argv, int* N, int* M) {
	*N = atoi(argv[1]);
	*M = atoi(argv[2]);
}


void GetBorders(const Matrix& m, 
	V& border_x0,
	V& border_x1,
	V& border_y0,
	V& border_y1
	) {
		for (int j = 0; j < m.M; ++j) {
			border_x0[j] = m.cat(0, j);
		}		

		for (int j = 0; j < m.M; ++j) {
			border_x1[j] = m.cat(m.N - 1, j);
		}		

		for (int i = 0; i < m.N; ++i) {
			border_y0[i] = m.cat(i, 0);
		}		

		for (int i = 0; i < m.N; ++i) {
			border_y1[i] = m.cat(i, m.M - 1);
		}		
}


void SyncBorder(const V& s, V& r, int other_i, int other_j) {
	
	if (other_i < 0 || other_i >= lattice_n) {
		r.resize(r.size(), 0.0);
		return;
	}	

	if (other_j < 0 || other_j >= lattice_m) {
		r.resize(r.size(), 0.0);
		return;
	}	



	int me = world_rank;
	int dst = ComputeLatticeId(other_i, other_j);
	

	MPI_Send(s.data(), s.size(), MPI_DOUBLE, dst, 0, MPI_COMM_WORLD);
	MPI_Status status;
	
        MPI_Recv(r.data(), r.size(), MPI_DOUBLE, dst, 0, MPI_COMM_WORLD, &status);
	
} 

double SyncDouble(double value) {
	double result;
	MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
	return result;
}


Matrix Solve(LocalOperator &op, int max_iter, double eps=1e-6) {
	Matrix w(op.block_h, op.block_w, 0.0);

	V border_x0_s(op.block_w, 0);
	V border_x1_s(op.block_w, 0);

	V border_y0_s(op.block_h, 0);
	V border_y1_s(op.block_h, 0);

	V border_x0_r(op.block_w, 0);
	V border_x1_r(op.block_w, 0);

	V border_y0_r(op.block_h, 0);
	V border_y1_r(op.block_h, 0);
	
	for (int iter = 0; iter < max_iter; ++iter) {
		GetBorders(w, border_x0_s, border_x1_s, border_y0_s, border_y1_s);
	
		SyncBorder(border_x0_s, border_x0_r, my_i - 1, my_j);
		SyncBorder(border_x1_s, border_x1_r, my_i + 1, my_j);
		SyncBorder(border_y0_s, border_y0_r, my_i, my_j - 1);
		SyncBorder(border_y1_s, border_y1_r, my_i, my_j + 1);
		

		Matrix Aw = op.Call(w, border_x0_r, border_x1_r, border_y0_r, border_y1_r);
		
		
	
		Matrix r = Aw - op.rhs;
		GetBorders(r, border_x0_s, border_x1_s, border_y0_s, border_y1_s);
	
		SyncBorder(border_x0_s, border_x0_r, my_i - 1, my_j);
		SyncBorder(border_x1_s, border_x1_r, my_i + 1, my_j);
		SyncBorder(border_y0_s, border_y0_r, my_i, my_j - 1);
		SyncBorder(border_y1_s, border_y1_r, my_i, my_j + 1);
		Matrix Ar = op.Call(r, border_x0_r, border_x1_r, border_y0_r, border_y1_r);

		
		
		
		double rAr = SyncDouble(op.Dot(r, Ar));
		double ArAr = SyncDouble(op.Dot(Ar, Ar));
		double rr = SyncDouble(op.Dot(r, r));
		
		double tau = rAr/ArAr;
		
		
		w = w - r * tau;
		double d = sqrt(rr) * tau;
		double delta = SyncDouble(op.Norm2(w - op.phi));
		
		if (world_rank == 0) {
			out << iter << "\t";
			out << sqrt(rr) << "\t";
			out << d << "\t";
			out << sqrt(delta) << "\t";
			out << tau;
			out << endl;
			cout << out.str();
			out.str("");
			out.clear();
		}


		if (d  < eps) {
			break;
		}

		
		
	}
	MPI_Barrier(MPI_COMM_WORLD);

	return w;
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

	int block_h, block_w;
	
	ComputeLatticeCoord(world_rank, &my_i, &my_j);

	ComputeBlockSize(N, M, &block_h, &block_w);
	int row_pos, col_pos;
	
	row_pos = N / lattice_n * my_i;
	col_pos = M / lattice_m * my_j;

	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	LocalOperator op(row_pos, col_pos, N, M, block_h, block_w, X_MIN, X_MAX, Y_MIN, Y_MAX);
	
	Matrix my_w = Solve(op, 10000000, 1e-6); 

	MPI_Finalize();
}
