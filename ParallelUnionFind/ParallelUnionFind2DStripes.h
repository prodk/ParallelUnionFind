#ifndef PARALLEL_UNION_FIND_2D_STRIPES
#define PARALLEL_UNION_FIND_2D_STRIPES

#include "parallelunionfindimpl.h"

class ParallelUnionFind2DStripes :
    public ParallelUnionFindImpl
{
public:
    ParallelUnionFind2DStripes(const DecompositionInfo& info);
    ~ParallelUnionFind2DStripes(void);

    // Overriden functions participating in the Template Method pattern
    void runLocalUnionFind(void);
    void constructGlobalLabeling(void);
    void mergeLabelsAcrossProcessors(void);
    void performFinalLabelingOfClusters(void);

private:
    std::size_t mNumOfPixels;
    std::tr1::shared_ptr<WeightedUnionFind> mWuf;
    std::vector<int> mPixels;

private:
    void copyPixels();
    inline int indexTo1D(int ix, int iy) const;
    inline void mergePixels(int idq, int idp);
};

//---------------------------------------------------------------------------
inline int ParallelUnionFind2DStripes::indexTo1D(int ix, int iy) const
{
    return ix * mDecompositionInfo.indexFactor * mDecompositionInfo.domainHeight + iy;
}

//---------------------------------------------------------------------------
inline void ParallelUnionFind2DStripes::mergePixels(int idq, int idp)
{
    if(mDecompositionInfo.pixelValue == mPixels[idq])
    {
        mWuf->setInitialRoot(idq);       // Specify the non-zero size (if it was 0) and init root.
        if( !mWuf->connected(idp, idq) ) // Merge verteces only if they're not yet merged.
        {
            mWuf->makeUnion(idp, idq);
        }
    }
}

#endif // PARALLEL_UNION_FIND_2D_STRIPES
