#ifndef PARALLEL_UNION_FIND_2D_STRIPES
#define PARALLEL_UNION_FIND_2D_STRIPES

#include "parallelunionfindimpl.h"

class ParallelUnionFind2DStripes : public ParallelUnionFindImpl
{
public:
    ParallelUnionFind2DStripes(const DecompositionInfo& info);
    ~ParallelUnionFind2DStripes(void);

    // Overriden functions participating in the Template Method pattern
    void runLocalUnionFind(void);
    void constructGlobalLabeling(void);
    void mergeLabelsAcrossProcessors(void);
    void performFinalLabelingOfClusters(void);

    // Overriden output.
    void printClusterSizes(const std::string& fileName) const;
    void printClusterStatistics(const std::string& fileName) const;
    void printClusterSizeHistogram(const int bins, const std::string& fileName) const;

    // Overriden helper interface functions. Virual functions CANNOT be inlined!
    void setPixelValue(const int value);

private:
    void copyPixels();
    inline int indexTo1D(int ix, int iy) const;
    inline void mergePixels(int idq, int idp);

private:
    std::size_t mNumOfPixels;
    std::tr1::shared_ptr<WeightedUnionFind> mLocalWuf; // UF DS with local (per-processor) labeling.
    std::tr1::shared_ptr<WeightedUnionFind> mGlobalWuf;// UF DS with global labeling.
    std::vector<int> mPixels;
    std::vector<int> mGlobalLabels;     // Global labels of 
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
        mLocalWuf->setInitialRoot(idq);       // Specify the non-zero size (if it was 0) and init root.
        if( !mLocalWuf->connected(idp, idq) ) // Merge verteces only if they're not yet merged.
        {
            mLocalWuf->makeUnion(idp, idq);
        }
    }
}

#endif // PARALLEL_UNION_FIND_2D_STRIPES
