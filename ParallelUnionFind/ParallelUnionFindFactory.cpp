#include "ParallelUnionFindFactory.h"

ParallelUnionFindFactory::ParallelUnionFindFactory(void)
{
}

ParallelUnionFindFactory::~ParallelUnionFindFactory(void)
{
}

ParallelUnionFindImpl* ParallelUnionFindFactory::makeParallelUnionFind(const std::string& spatialConfiguration)
{
    if(spatialConfiguration == "2DStripes")
    {
        return new ParallelUnionFind2DStripes();
    }
    else
    {
        return NULL;
    }
}
