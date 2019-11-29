#include "mpi.h"
#include <math.h>
#include <iostream>
#include <vector>

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

class LocalOperator {
public:
	const int N, M;
protected:
	int row_pos, col_pos;
	int block_h, block_w;
public:
	LocalOperator(int my_i, int my_j, int block_h_, int block_w_, int N_, int M_):
		N(N_), M(M_)
	 {
		row_pos = my_i;
		col_pos = my_j;
		block_h = block_h_;
		block_w = block_w_;
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
	*N = 10;
	*M = 10;
}

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
	
	LocalOperator op(my_i, my_j, N, M, block_h, block_w);
	
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
