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

void SendBorder(const vector<double>& border, int dst, int me) {
	MPI_Send(border.data(), border.size(), MPI_DOUBLE, dst, me, COMM_BORDER_SYNC);
}

void RecvBorder(vector<double>& border, int src) {
	MPI_Status status;
	MPI_Recv(border.data(), border.size(), MPI_DOUBLE, src, src, COMM_BORDER_SYNC, &status);
}

vector<double> SyncBorder(int my_i, int my_j, int other_i, int other_j, vector<double> my_border) {
	vector<double> other_border(my_border.size());

	
	SendBorder(my_border, ComputeLatticeId(other_i, other_j), ComputeLatticeId(my_i, my_j));
	RecvBorder(other_border, ComputeLatticeId(other_i, other_j));

	return other_border;
}

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


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(COMM_BORDER_SYNC, &world_rank);
	MPI_Comm_size(COMM_BORDER_SYNC, &world_size);
	
	InitLocalOperator();
	auto my_w = Solve();
	SendLocalResultsTo0(my_w);

	if(world_rank == 0) {
		auto result = ReceiveLocalResults();
		DoIOStuff();
	}

	MPI_Finalize();
}
