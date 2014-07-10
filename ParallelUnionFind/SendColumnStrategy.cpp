// SendColumnStrategy.cpp

#include "SendColumnStrategy.h"
#include "ParallelUnionFindImpl.h"

SendColumnStrategy::SendColumnStrategy(const DecompositionInfo & decompositionInfo, std::vector<Pixel> & pixels)
    : mDecompositionInfo(decompositionInfo)
    , mGlobalPixels(pixels)
{
}


SendColumnStrategy::~SendColumnStrategy(void)
{
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
            return -1; // TODO: get rid of magic numbers, use enum hack instead.
        }
    }

    return -1; // TODO: get rid of magic numbers, use enum hack instead.
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
            return -1; // TODO: get rid of magic numbers, use enum hack instead.
        }
    }

    return -1; // TODO: get rid of magic numbers, use enum hack instead.
}

//---------------------------------------------------------------------------
bool SendColumnStrategy::isNeighborProcessorValid(const int rank) const
{
    return (rank >= 0);
}
