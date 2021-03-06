// ParallelUnionFindFactory.h - a factory for creating the desired PUF implementation.

#ifndef PARALLEL_UNION_FIND_FACTORY
#define PARALLEL_UNION_FIND_FACTORY

//---------------------------------------------------------------------------
#include <string>

//---------------------------------------------------------------------------
class ParallelUnionFindImpl;
struct DecompositionInfo;

//---------------------------------------------------------------------------
class ParallelUnionFindFactory
{
public:
    ParallelUnionFindFactory(void);
    ~ParallelUnionFindFactory(void);

    static ParallelUnionFindImpl*
    makeParallelUnionFind(const std::string& spatialConfiguration, const DecompositionInfo& info);
};

#endif // PARALLEL_UNION_FIND_FACTORY