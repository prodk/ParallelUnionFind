// ParallelUnionFind2DStripes.h - declaration of the class ParallelUnionFind2DStripes.

#ifndef PARALLEL_UNION_FIND_2D_STRIPES
#define PARALLEL_UNION_FIND_2D_STRIPES

//---------------------------------------------------------------------------
#include "ParallelUnionFindImpl.h"
#include "WeightedUnionFind.h"
#include <memory>

//---------------------------------------------------------------------------
class ParallelUnionFind2DStripes : public ParallelUnionFindImpl
{
public:
    ParallelUnionFind2DStripes(const DecompositionInfo& info);
    ~ParallelUnionFind2DStripes(void);

    // Overriden output.
    // TODO: rename/add functions that specify whether local or global quantities are printed.
    void printClusterSizes(const std::string& fileName) const;
    void printClusterStatistics(const std::string& fileName) const;
    void printClusterSizeHistogram(const int bins, const std::string& fileName) const;

    // Overriden helper interface functions. Virual functions CANNOT be inlined!
    void setPixelValue(const int value);

private:
    // Overriden functions participating in the Template Method pattern.
    void runLocalUnionFind(void);
    void constructGlobalLabeling(void);
    void mergeLabelsAcrossProcessors(void);
    void performFinalLabelingOfClusters(void);

private:
    void copyPixels();
    inline int indexTo1D(int ix, int iy) const;
    inline void mergePixels(int idq, int idp);

    // Implementation helpers.
    void printLocalExtendedPicture(const DecompositionInfo& info) const;

    // Stage 1 helpers.
    int getNeighborPeriodicBC(const int index, const int size) const;    // Next neighbor assuming the PBC.
    int getNeighborNonPeriodicBC(const int index, const int size) const; // Next neighbor, no PBC. Returns -1 at the border.

    // Stage 2 helpers.
    int receiveNumberOfClustersFromPreviousProcs() const;
    void sendTotalClustersToNextProcs(const int numOfClustersOnSmallerProcIds, const int numOfMyClusters) const;

    // Stage 3 helpers.
    void setLocalPartOfGloblaPixels(void);
    void copyLeftColumnAndSendToLeftNeighbor(void);
    void copyRightColumnAndSendToRightNeighbor(void);
    void copyColumn(void);
    void sendColumn(void);


private:
    std::size_t mNumOfPixels;
    std::size_t mNumOfGlobalPixels;
    std::tr1::shared_ptr<WeightedUnionFind> mLocalWuf; // UF with local (per-processor) labeling.
    std::tr1::shared_ptr<WeightedUnionFind> mGlobalWuf;// UF with global labeling.
    std::vector<int> mLocalPixels;                     // Local mesh points.
    std::map<int, int> mGlobalLabels;                  // Non-consecutive local root is a key, a consecutive global root is a value.
    std::vector<Pixel> mGlobalPixels;                  // Extended pixels.
};

//---------------------------------------------------------------------------
inline int ParallelUnionFind2DStripes::indexTo1D(int ix, int iy) const
{
    return ix*mDecompositionInfo.domainHeight + iy;
}

//---------------------------------------------------------------------------
inline void ParallelUnionFind2DStripes::mergePixels(int idq, int idp)
{
    if (mDecompositionInfo.pixelValue == mLocalPixels[idq])
    {
        mLocalWuf->setInitialRoot(idq);       // Specify the non-zero size (if it was 0) and init root.
        if (!mLocalWuf->connected(idp, idq))  // Merge vertices only if they're not yet merged.
        {
            mLocalWuf->makeUnion(idp, idq);
        }
    }
}

//---------------------------------------------------------------------------
inline int ParallelUnionFind2DStripes::getNeighborPeriodicBC(const int index, const int size) const
{
    return (index + 1) % size;
}

//---------------------------------------------------------------------------
inline int ParallelUnionFind2DStripes::getNeighborNonPeriodicBC(const int index, const int size) const
{
    return (index >= size - 1) 
           ? -1 
           : index + 1;
}

#endif // PARALLEL_UNION_FIND_2D_STRIPES
