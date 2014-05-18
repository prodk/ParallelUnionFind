// ParallelUnionFind.cpp - implementation of ParalleUnionFind class

#include "ParallelUnionFind.h"
#include "ParallelUnionFindFactory.h"

ParallelUnionFind::ParallelUnionFind(const std::string& spatialConfiguration,
                                     const std::size_t N,
                                     int numOfProc,
                                     int myRank):
mNumOfProc(numOfProc),
mMyRank(myRank),
mParallelUf(ParallelUnionFindFactory::makeParallelUnionFind(spatialConfiguration, N, numOfProc, myRank))
{
}

ParallelUnionFind::~ParallelUnionFind(void)
{
}

void ParallelUnionFind::analyze(void)
{
#ifdef _DEBUG
    logMsg(std::cout, "PUF::analyze()");
#endif

    mParallelUf->analyze();
}

void ParallelUnionFind::printClusterStatistics(const std::string& fileName) const
{
#ifdef _DEBUG
    logMsg(std::cout, "PUF::printClusterStatistics()");
#endif
}

void ParallelUnionFind::printClusterSizeHistogram(const std::string& fileName) const
{
#ifdef _DEBUG
    logMsg(std::cout, "PUF::printClusterSizeHistogram()");
#endif
}

#ifdef _DEBUG
void ParallelUnionFind::logMsg(std::ostream& out, const std::string& msg) const
{
    out << msg << std::endl;
}
#endif
