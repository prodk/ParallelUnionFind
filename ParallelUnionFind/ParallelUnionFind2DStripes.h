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

    // Stage 1 helpers.
    int getNeighborPeriodicBC(const int index, const int size) const;    // Next neighbor assuming the PBC.
    int getNeighborNonPeriodicBC(const int index, const int size) const; // Next neighbor, no PBC. Returns -1 at the border.

    // Stage 2 helpers.
    int receiveNumberOfClustersFromPreviousProcs() const;
    void sendTotalClustersToNextProcs(const int numOfClustersOnSmallerProcIds, const int numOfMyClusters) const;

    // Stage 3 helpers.
    void setLocalPartOfGloblaPixels(void); // TODO: rename this function (see implementation for why).

    // TODO: think about refactoring of left/right stripes to avoid code duplication.
    // TODO: encapsulate this sending/receiving code into a strategy class.
    // At the moment I decided to hard code this for performance/clarity reasons.
    void copyLeftColumnAndSendToLeftNeighbor(void);

    // A helper data structure containing attributes of the pixel stripe sent between processors.
    struct SPixelStripe
    {
        std::vector<int> pixelValue;
        std::vector<int> globalClusterId;
        std::vector<int> sizeOfCluster;

        explicit SPixelStripe(size_t size)
            : pixelValue(size)
            , globalClusterId(size)
            , sizeOfCluster(size)
        {}
    };

    void copyLeftPixelStripeToSend(SPixelStripe & stripeToSend);

    void sendLeftStripeFromEvenReceiveOnOdd(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const;
    void sendLeftStripeFromOddReceiveOnEven(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const;

    void sendLeftStripe(SPixelStripe & stripeToSend, const int msgId[], const int size) const;
    void receiveLeftStripe(SPixelStripe & stripeToReceive, const int msgId[], const int size) const;

    void saveReceivedStripeToRightStripe(const SPixelStripe & stripeToReceive);

    // Coping with the right stripe.
    void copyRightColumnAndSendToRightNeighbor(void);
    void copyRightPixelStripeToSend(SPixelStripe & stripeToSend);

    void sendRightStripeFromEvenReceiveOnOdd(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const;
    void sendRightStripeFromOddReceiveOnEven(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const;

    void sendRightStripe(SPixelStripe & stripeToSend, const int msgId[], const int size) const;
    void receiveRightStripe(SPixelStripe & stripeToReceive, const int msgId[], const int size) const;

    void saveReceivedStripeToLeftStripe(const SPixelStripe & stripeToReceive);

    // Used both for left/right stripes
    int getLeftNeighborProcessor() const;      // Periodic BCs are taken into account via DecompositionInfo.
    int getRightNeighborProcessor() const;     // Periodic BCs are taken into account via DecompositionInfo.
    bool isNeighborProcessorValid(const int rank) const;

    void runUfOnGlobalLabelsAndRecordMerges();
    void unionExtendedPixelsUsingGlobalLabels();

    // Implementation helpers.
    void printLocalExtendedPicture(const DecompositionInfo& info) const;
    void printReceivedGlobalLabels() const;

    void printGlobalUfRootsAfterFirstMerge() const;

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
