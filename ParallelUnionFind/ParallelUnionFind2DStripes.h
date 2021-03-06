// ParallelUnionFind2DStripes.h - declaration of the class ParallelUnionFind2DStripes.

#ifndef PARALLEL_UNION_FIND_2D_STRIPES
#define PARALLEL_UNION_FIND_2D_STRIPES

//---------------------------------------------------------------------------
#include "ParallelUnionFindImpl.h"
#include "WeightedUnionFind.h"
#include <memory>

#ifndef _WIN32
#include <tr1/memory>  // tr1::shared_ptr for gcc
#endif

//---------------------------------------------------------------------------
class ParallelUnionFind2DStripes : public ParallelUnionFindImpl
{
public:
    ParallelUnionFind2DStripes(const DecompositionInfo& info);
    ~ParallelUnionFind2DStripes(void);

    // Overriden output.
    void printClusterSizes(const std::string& fileName) const;
    void printClusterStatistics(const std::string& fileName) const;
    void printClusterSizeHistogram(const int bins, const std::string& fileName) const;

    // Overriden output of local data.
    void printPerProcessorClusterSizes(const std::string& fileName) const;
    void printPerProcessorClusterStatistics(const std::string& fileName) const;
    void printPerProcessorClusterSizeHistogram(const int bins, const std::string& fileName) const;

    // Overriden helper interface functions. Virtual functions CANNOT be inlined!
    void setPixelValue(const int value);

private:
    // Overriden functions participating in the Template Method pattern.
    void runLocalUnionFind(void);
    void constructGlobalLabeling(void);
    void mergeLabelsAcrossProcessors(void);
    void performFinalLabelingOfClusters(void);

private:
    void printInputInfo() const;
    void copyPixels();
    inline std::ptrdiff_t indexTo1D(int ix, int iy) const;
    inline void mergePixels(std::ptrdiff_t idq, std::ptrdiff_t idp, std::tr1::shared_ptr<WeightedUnionFind> wuf,
                            const int pixelValue) const;

    // Stage 1 helpers.
    int getNeighborPeriodicBC(const int index, const int size) const;    // Next neighbor assuming the PBC.
    int getNeighborNonPeriodicBC(const int index, const int size) const; // Next neighbor, no PBC. Returns -1 at the border.

    // Stage 2 helpers.
    int receiveNumberOfClustersFromPreviousProcs() const;
    void sendTotalClustersToNextProcs(const std::ptrdiff_t numOfClustersOnSmallerProcIds, const std::ptrdiff_t numOfMyClusters) const;

    // Stage 3 helpers.
    void initializeGlobalPixels(void);
    void copyLeftColumnAndSendToLeftNeighbor(void);
    void copyRightColumnAndSendToRightNeighbor(void);
    void runUfOnGlobalLabels();
    void runLocalUfOnGlobalLabelsToSetInitialRoots();
    void runUfOnGlobalPixelsAndRecordGlobalMerges();

    void mergeClusterIds(std::ptrdiff_t idq, std::ptrdiff_t idp, std::tr1::shared_ptr<WeightedUnionFind> wuf) const;
    void mergePixelsAndRecordMerge(std::ptrdiff_t idq, std::ptrdiff_t idp, std::tr1::shared_ptr<WeightedUnionFind> wuf,
                                   const int pixelValue);
    void recordMerge(const std::ptrdiff_t idp, const std::ptrdiff_t idq);

    bool isPixelValid(const int pixel) const;
    bool isClusterIdValid(const std::ptrdiff_t pixel) const;
    bool isBoss() const;

    // Stage 4 helpers.
    void getMergesFromAllProcs();
    void broadcastMergeAndAddToAllMerges(const std::vector<std::ptrdiff_t> & arrayToSend, std::vector<std::ptrdiff_t> & arrayToReceive,int numOfMerges, int root);
    int calculateNumberOfGlobalLabels() const;
    void getUniqueLabelForEachComponent();

    // Data analysis.
    void getMinMaxClusterSizes();
    void lookForPercolation();
    void lookForHorizontalPercolation();
    void bossLookForHorizontalPercolation();
    void slaveLookForHorizontalPercolation();
    void lookForHorizontalPercolation1Proc();
    void lookForVerticalPercolation();
    void bossLookForVerticalPercolation(std::set<int>& topmostRoots, std::set<int>& bottommostRoots);

    std::ptrdiff_t getFinalRootOfPixel(const int pixelId);
    std::ptrdiff_t getLocalRootFromGloablRoot(const std::ptrdiff_t globalRoot) const;
    void printPercolationInfo() const;
    void printPercolationPhrase(const std::string& contact, const std::string& vh, const int size) const;

    void adjustFinalHistogram(const int bins, std::vector<double>& finalHistogram,
                              std::multimap<std::ptrdiff_t, std::ptrdiff_t>& rootsInBin) const;
    void adjustFinalHistogramIfSlave(const int bins, std::multimap<std::ptrdiff_t, std::ptrdiff_t>& rootsInBin) const;
    void adjustFinalHistogramIfBoss(const int bins,
                                    std::vector<double>& finalHistogram,
                                    std::multimap<std::ptrdiff_t, std::ptrdiff_t>& rootsInBin) const;

    void outputSizeHistogram(const int bins, const double binWidth, const std::string& fileName, const std::vector<double>& finalHistogram) const;

    void packDataAndSendIt(const std::set<int>& data, const int destinationProc) const;
    void receiveData(std::vector<std::ptrdiff_t>& data, int& numOfElements, const int sendingProc) const;

    // Implementation helpers.
    void printConsecutiveIds(const std::map<std::ptrdiff_t, std::ptrdiff_t>& consecutiveLocalIds);
    void printLocalExtendedPicture(const DecompositionInfo& info) const;
    void printReceivedGlobalLabels() const;
    void printGlobalUfRootsAfterFirstMerge() const;
    void printMerges() const;

private:
    std::ptrdiff_t mNumOfPixels;
    std::ptrdiff_t mNumOfGlobalPixels;
    std::tr1::shared_ptr<WeightedUnionFind> mLocalWuf;      // UF with local (per-processor) labeling.
    std::tr1::shared_ptr<WeightedUnionFind> mGlobalWuf;     // UF with global labeling.
    std::vector<std::ptrdiff_t> mLocalPixels;               // Local mesh points.
    std::map<std::ptrdiff_t, std::ptrdiff_t> mGlobalLabels; // Non-consecutive local root is a key, a consecutive global root is a value.
    std::vector<Pixel> mGlobalPixels;                       // Extended pixels.
    Merge mMerge;                                           // A merge of 2 clusters residing on different processors but belonging to one cluster.
    Merge mAllMerges;                                       // All merges (including the one from the current proc).
    std::tr1::shared_ptr<WeightedUnionFind> mFinalWuf;      // UF with final labeling, includes merges that span several procs.
    std::ptrdiff_t mTotalNumOfClusters;                     // Number of connected components on all the processors.
    std::ptrdiff_t mMinClusterSize;                         // Global value of the smallest cluster size.
    std::ptrdiff_t mMaxClusterSize;                         // Global value of the largest cluster size.
    std::map<std::ptrdiff_t, std::ptrdiff_t> mClusterSizes; // Cluster size is a key, root is a value. Important: use map to avoid duplicate copies.
    bool mPercolatesHorizontally;                           // True if percolation in horizontal direction takes place.
    bool mPercolatesVertically;
    std::ptrdiff_t mHorizPercolatedSize;                    // Size of the horizontally percolated cluster.
    std::ptrdiff_t mVertPercolatedSize;                     // Size of the vertically percolated cluster.

    enum {INVALID_VALUE = -1, BOSS, MSG_1};                 // BOSS is 0 by default.
};

//---------------------------------------------------------------------------
inline std::ptrdiff_t ParallelUnionFind2DStripes::indexTo1D(int ix, int iy) const
{
    return ix*mDecompositionInfo.domainHeight + iy;
}

//---------------------------------------------------------------------------
inline void ParallelUnionFind2DStripes::mergePixels(std::ptrdiff_t idq, std::ptrdiff_t idp, std::tr1::shared_ptr<WeightedUnionFind> wuf,
                                                    const int pixelValue) const
{
    if (mDecompositionInfo.pixelValue == pixelValue)
    {
        wuf->setInitialRoot(idq);       // Specify the non-zero size (if it was 0) and init root.
        if (!wuf->connected(idp, idq))  // Merge vertices only if they're not yet merged.
        {
            wuf->makeUnion(idp, idq);
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
           ? INVALID_VALUE
           : index + 1;
}

//---------------------------------------------------------------------------
inline bool ParallelUnionFind2DStripes::isPixelValid(const int pixel) const
{
    return (pixel != INVALID_VALUE);
}

//---------------------------------------------------------------------------
inline bool ParallelUnionFind2DStripes::isClusterIdValid(const std::ptrdiff_t pixel) const
{
    return (mGlobalPixels[pixel].globalClusterId > INVALID_VALUE);
}

//---------------------------------------------------------------------------
inline bool ParallelUnionFind2DStripes::isBoss() const
{
    return (BOSS == mDecompositionInfo.myRank);
}

#endif // PARALLEL_UNION_FIND_2D_STRIPES
