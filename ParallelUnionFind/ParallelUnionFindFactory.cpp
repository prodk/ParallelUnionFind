#include "ParallelUnionFindFactory.h"

ParallelUnionFindFactory::ParallelUnionFindFactory(void)
{
}

ParallelUnionFindFactory::~ParallelUnionFindFactory(void)
{
}

ParallelUnionFindImpl* 
ParallelUnionFindFactory::makeParallelUnionFind(const std::string& spatialConfiguration,
                                                const std::size_t N,
                                                int numOfProc,
                                                int myRank)
{
    if(spatialConfiguration == "2DStripes")
    {
        return new ParallelUnionFind2DStripes(N, numOfProc, myRank);
    }
    else
    {
        return NULL;
    }
}
