#include "TestParallelUnionFind.h"

TestParallelUnionFind::TestParallelUnionFind(int argc, char **argv)
{
    mNumOfProc = -1;
    mMyRank = -1;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &mNumOfProc);
    MPI_Comm_rank (MPI_COMM_WORLD, &mMyRank);
}

TestParallelUnionFind::~TestParallelUnionFind(void)
{
    // Cleanup MPI.
    MPI_Finalize();
}

void TestParallelUnionFind::runTests()
{
    std::cout << "Processor " << mMyRank << ": number of processors is "
        << mNumOfProc << std::endl;

    // DEBUG: force the debug window to remain.
    if(0 == mMyRank)
    {
        std::cin >> mMyRank;
    }
    // DEBUG end
}
