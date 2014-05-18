#ifndef PARALLEL_UNION_FIND_FACTORY
#define PARALLEL_UNION_FIND_FACTORY

#include "ParallelUnionFind2DStripes.h"

class ParallelUnionFindFactory
{
public:
    ParallelUnionFindFactory(void);
    ~ParallelUnionFindFactory(void);

    static ParallelUnionFindImpl*
    makeParallelUnionFind(const std::string& spatialConfiguration, const DecompositionInfo& info);
};

#endif // PARALLEL_UNION_FIND_FACTORY
