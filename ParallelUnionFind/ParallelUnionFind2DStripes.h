// ParallelUnionFind2DStripes.h - declaration of the class ParallelUnionFind2DStripes.

#ifndef PARALLEL_UNION_FIND_2D_STRIPES
#define PARALLEL_UNION_FIND_2D_STRIPES

//---------------------------------------------------------------------------
#include "ParallelUnionFindImpl.h"
#include "WeightedUnionFind.h"
#include <memory>

// TODO: replace all the integers that relate to the system size by ptrdiff_t (necessary for the 64-bit JUQUEEN and huge systems).

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
    void copyPixels();
    inline int indexTo1D(int ix, int iy) const;
    inline void mergePixels(int idq, int idp, std::tr1::shared_ptr<WeightedUnionFind> wuf,
                            const int pixelValue) const;

    // Stage 1 helpers.
    int getNeighborPeriodicBC(const int index, const int size) const;    // Next neighbor assuming the PBC.
    int getNeighborNonPeriodicBC(const int index, const int size) const; // Next neighbor, no PBC. Returns -1 at the border.

    // Stage 2 helpers.
    int receiveNumberOfClustersFromPreviousProcs() const;
    void sendTotalClustersToNextProcs(const int numOfClustersOnSmallerProcIds, const int numOfMyClusters) const;

    // Stage 3 helpers.
    void initializeGloblaPixels(void);
    void copyLeftColumnAndSendToLeftNeighbor(void);
    void copyRightColumnAndSendToRightNeighbor(void);
    void runUfOnGlobalLabels();
    void runLocalUfOnGlobalLabelsToSetInitialRoots();
    void runUfOnGlobalPixelsAndRecordGlobalMerges();

    void mergeClusterIds(int idq, int idp, std::tr1::shared_ptr<WeightedUnionFind> wuf) const;
    void mergePixelsAndRecordMerge(int idq, int idp, std::tr1::shared_ptr<WeightedUnionFind> wuf,
                                    const int pixelValue);
    void recordMerge(const int idp, const int idq);

    bool isPixelValid(const int pixel) const;
    bool isClusterIdValid(const int pixel) const;

    // Stage 4 helpers.
    void getMergesFromAllProcs();
    void broadcastMergeAndAddToAllMerges(const std::vector<int> & arrayToSend, std::vector<int> & arrayToReceive,int numOfMerges, int root);
    int calculateNumberOfGlobalLabels() const;
    void getUniqueLabelForEachComponent();

    // Data analysis.
    void getMinMaxClusterSizes();
    void lookForPercolation();
    void lookForHorizontalPercolation();
    void lookForVerticalPercolation();
    int getFinalRootOfPixel(const int pixelId);
    void printPercolationInfo() const;
    void printPercolationPhrase(const std::string& contact, const std::string& vh, const int size) const;

    void adjustFinalHistogram(const int bins, std::vector<double>& finalHistogram, std::multimap<int, int>& rootsInBin) const;
    void outputSizeHistogram(const int bins, const double binWidth, const std::string& fileName, const std::vector<double>& finalHistogram) const;

    void packDataAndSendIt(const std::set<int>& data, const int destinationProc) const;
    void receiveData(std::vector<int>& data, int& numOfElements, const int sendingProc) const;

    // Implementation helpers.
    void printLocalExtendedPicture(const DecompositionInfo& info) const;
    void printReceivedGlobalLabels() const;
    void printGlobalUfRootsAfterFirstMerge() const;
    void printMerges() const;

private:
    std::size_t mNumOfPixels;
    std::size_t mNumOfGlobalPixels;
    std::tr1::shared_ptr<WeightedUnionFind> mLocalWuf; // UF with local (per-processor) labeling.
    std::tr1::shared_ptr<WeightedUnionFind> mGlobalWuf;// UF with global labeling.
    std::vector<int> mLocalPixels;                     // Local mesh points.
    std::map<int, int> mGlobalLabels;                  // Non-consecutive local root is a key, a consecutive global root is a value.
    std::vector<Pixel> mGlobalPixels;                  // Extended pixels.
    Merge mMerge;                                      // A merge of 2 clusters residing on different processors but belonging to one cluster.
    Merge mAllMerges;                                  // All merges (including the one from the current proc).
    std::tr1::shared_ptr<WeightedUnionFind> mFinalWuf; // UF with final labeling, includes merges that span several procs.
    int mTotalNumOfClusters;                           // Number of connected components on all the processors.
    int mMinClusterSize;                               // Global value of the smallest cluster size.
    int mMaxClusterSize;                               // Global value of the largest cluster size.
    std::map<int, int> mClusterSizes;                  // Cluster size is a key, root is a value. Important: use map to avoid duplicate copies.
    bool mPercolatesHorizontally;                      // True if percolation in horizontal direction takes place.
    bool mPercolatesVertically;
    int mHorizPercolatedSize;                          // Size of the horizontally percolated cluster.
    int mVertPercolatedSize;                           // Size of the vertically percolated cluster.

    enum {INVALID_VALUE = -1, BOSS, MSG_1};            // BOSS is 0 by default.
};

//---------------------------------------------------------------------------
inline int ParallelUnionFind2DStripes::indexTo1D(int ix, int iy) const
{
    return ix*mDecompositionInfo.domainHeight + iy;
}

//---------------------------------------------------------------------------
inline void ParallelUnionFind2DStripes::mergePixels(int idq, int idp, std::tr1::shared_ptr<WeightedUnionFind> wuf,
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

#endif // PARALLEL_UNION_FIND_2D_STRIPES
