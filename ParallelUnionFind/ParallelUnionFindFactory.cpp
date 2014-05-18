#include "ParallelUnionFindFactory.h"

ParallelUnionFindFactory::ParallelUnionFindFactory(void)
{
}

ParallelUnionFindFactory::~ParallelUnionFindFactory(void)
{
}

ParallelUnionFindImpl* 
ParallelUnionFindFactory::makeParallelUnionFind(const std::string& spatialConfiguration,
                                                const DecompositionInfo& info)
{
    if(spatialConfiguration == "2DStripes")
    {
        return new ParallelUnionFind2DStripes(info);
    }
    else
    {
        return NULL;
    }
}
