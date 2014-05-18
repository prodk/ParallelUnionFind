#ifndef PARALLEL_UNION_FIND_2D_STRIPES
#define PARALLEL_UNION_FIND_2D_STRIPES

#include "parallelunionfindimpl.h"

class ParallelUnionFind2DStripes :
    public ParallelUnionFindImpl
{
public:
    ParallelUnionFind2DStripes(std::size_t N, int numOfProc, int myRank);
    ~ParallelUnionFind2DStripes(void);

    // Overriden functions participating in the Template Method pattern
    void runLocalUnionFind(void);
    void constructGlobalLabeling(void);
    void mergeLabelsAcrossProcessors(void);
    void performFinalLabelingOfClusters(void);

private:
    std::tr1::shared_ptr<WeightedUnionFind> wuf;
};

#endif // PARALLEL_UNION_FIND_2D_STRIPES
