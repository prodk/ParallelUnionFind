#ifndef PARALLEL_UNION_FIND_2D_STRIPES
#define PARALLEL_UNION_FIND_2D_STRIPES

#include "parallelunionfindimpl.h"

class ParallelUnionFind2DStripes :
    public ParallelUnionFindImpl
{
public:
    ParallelUnionFind2DStripes(void);
    ~ParallelUnionFind2DStripes(void);

    // Overriden functions participating in the Template Method pattern
    void runLocalUnionFind(void);
    void constructGlobalLabeling(void);
    void mergeLabelsAcrossProcessors(void);
    void performFinalLabelingOfClusters(void);
};

#endif // PARALLEL_UNION_FIND_2D_STRIPES
