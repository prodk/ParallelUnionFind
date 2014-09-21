// SendColumnStrategy.cpp

#include "SendColumnStrategy.h"
#include "ParallelUnionFindImpl.h"

//---------------------------------------------------------------------------
SendColumnStrategy::SendColumnStrategy(const DecompositionInfo & decompositionInfo,
                                       const std::vector<std::ptrdiff_t> & localPixels,
                                       const std::vector<Pixel> & globalPixels,
                                       const std::tr1::shared_ptr<WeightedUnionFind> & localWuf,
                                       std::map<std::ptrdiff_t, std::ptrdiff_t> & globalLabels)
    : mDecompositionInfo(decompositionInfo)
    , mLocalPixels(localPixels)
    , mGlobalPixels(globalPixels)
    , mLocalWuf(localWuf)
    , mGlobalLabels(globalLabels)
{
}

//---------------------------------------------------------------------------
SendColumnStrategy::~SendColumnStrategy(void)
{
}

//---------------------------------------------------------------------------
void SendColumnStrategy::sendReceivePixelStripes(std::vector<Pixel> & globalPixels)
{
    SPixelStripe stripeToSend(mDecompositionInfo.domainHeight);
    copyPixelStripeToSend(stripeToSend);

    SPixelStripe stripeToReceive(mDecompositionInfo.domainHeight);
    sendStripeFromEvenReceiveOnOdd(stripeToSend, stripeToReceive);
    sendStripeFromOddReceiveOnEven(stripeToSend, stripeToReceive);

    saveReceivedStripe(stripeToReceive, globalPixels);
}

//---------------------------------------------------------------------------
int SendColumnStrategy::getLeftNeighborProcessor() const
{
    if (mDecompositionInfo.periodicBoundaryX)
    {
        return (mDecompositionInfo.myRank - 1) % mDecompositionInfo.numOfProc;
    }
    else
    {
        if (mDecompositionInfo.myRank > 0)
        {
            return mDecompositionInfo.myRank - 1;
        }
        else
        {
            return INVALID_VALUE;
        }
    }

    return INVALID_VALUE;
}

//---------------------------------------------------------------------------
int SendColumnStrategy::getRightNeighborProcessor() const
{
    if (mDecompositionInfo.periodicBoundaryX)
    {
        return (mDecompositionInfo.myRank + 1) % mDecompositionInfo.numOfProc;
    }
    else
    {
        if (mDecompositionInfo.myRank < mDecompositionInfo.numOfProc - 1)
        {
            return mDecompositionInfo.myRank + 1;
        }
        else
        {
            return INVALID_VALUE;
        }
    }

    return INVALID_VALUE;
}

//---------------------------------------------------------------------------
bool SendColumnStrategy::isNeighborProcessorValid(const int rank) const
{
    return (rank >= BOSS);
}
