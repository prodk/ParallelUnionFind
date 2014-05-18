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
    test1();
#ifdef _DEBUG
    forceWindowToStay();
#endif
}

void TestParallelUnionFind::test1()
{
    ParallelUnionFind puf("2DStripes", defineDecomposition());
    puf.analyze();
    puf.printClusterStatistics("fileName.txt");
}

DecompositionInfo TestParallelUnionFind::defineDecomposition()
{
    DecompositionInfo info;

    info.domainWidth = 256;
    info.domainHeight = 256;
    info.myRank = mMyRank;
    info.numOfProc = mNumOfProc;

    return info;
}

#ifdef _DEBUG
void TestParallelUnionFind::forceWindowToStay() const
{
    if(0 == mMyRank)
    {
        int dummy;
        std::cin >> dummy;
    }
}
#endif
