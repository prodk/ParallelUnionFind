// ParallelUnionFind.h - declaration of ParalleUnionFind class

#ifndef PARALLEL_UNION_FIND
#define PARALLEL_UNION_FIND

#include "ParallelUnionFindImpl.h"

class ParallelUnionFind
{
public:
    ParallelUnionFind(void);
    ~ParallelUnionFind(void);
    ParallelUnionFind(const std::string& spatialConfiguration);

    void analyze(void);
    void printClusterStatistics(const std::string& fileName) const;
    void printClusterSizeHistogram(const std::string& fileName) const;

private:
    std::tr1::shared_ptr<ParallelUnionFindImpl> mParallelUf;
};

#endif // PARALLEL_UNION_FIND
