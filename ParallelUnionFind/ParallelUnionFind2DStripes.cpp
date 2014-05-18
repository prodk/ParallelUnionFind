#include "ParallelUnionFind2DStripes.h"

ParallelUnionFind2DStripes::ParallelUnionFind2DStripes(std::size_t N, int numOfProc, int myRank):
ParallelUnionFindImpl(numOfProc, myRank),
wuf(new WeightedUnionFind(N))
{
    if(mNumOfProc <= 0)
    {
        std::cerr << "0 number of processors!" << std::endl;
        std::cerr << "Check whether the MPI has been initialized!" << std::endl;
    }
    else
    {
        std::cout << "Processor " << mMyRank << " of " << mNumOfProc << std::endl;
        std::cout << "PUF of type 2DStripes created with " << N << " elements" << std::endl;
    }
}

ParallelUnionFind2DStripes::~ParallelUnionFind2DStripes(void)
{
}

void ParallelUnionFind2DStripes::runLocalUnionFind(void)
{
}

void ParallelUnionFind2DStripes::constructGlobalLabeling(void)
{
}

void ParallelUnionFind2DStripes::mergeLabelsAcrossProcessors(void)
{
}

void ParallelUnionFind2DStripes::performFinalLabelingOfClusters(void)
{
}
