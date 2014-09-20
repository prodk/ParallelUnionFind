// SendColumnStrategy.cpp - the implementation of the SendColumnStrategy class.

#include "SendLeftColumnStrategy.h"
#include "WeightedUnionFind.h"
#include "ParallelUnionFindImpl.h"
#include <mpi.h>

//---------------------------------------------------------------------------
SendLeftColumnStrategy::SendLeftColumnStrategy(const DecompositionInfo & decompositionInfo,
                                               const std::vector<int> & localPixels,
                                               const std::vector<Pixel> & globalPixels,
                                               const std::tr1::shared_ptr<WeightedUnionFind> & localWuf,
                                               std::map<std::ptrdiff_t, std::ptrdiff_t> & globalLabels)
    : SendColumnStrategy(decompositionInfo, localPixels, globalPixels, localWuf, globalLabels)
{
}

//---------------------------------------------------------------------------
SendLeftColumnStrategy::~SendLeftColumnStrategy(void)
{
}

//---------------------------------------------------------------------------
// Copy left pixel stripe.
void SendLeftColumnStrategy::copyPixelStripeToSend(SPixelStripe & stripeToSend) const
{
    for (std::ptrdiff_t iy = 0u; iy < mDecompositionInfo.domainHeight; ++iy)
    {
        stripeToSend.pixelValue[iy] = mLocalPixels[iy];

        // Set attributes only of those pixels that have the desired value.
        if (mDecompositionInfo.pixelValue == stripeToSend.pixelValue[iy])
        {
            const std::ptrdiff_t pixelRoot = mLocalWuf->getPixelRoot(iy);
            stripeToSend.globalClusterId[iy] = mGlobalLabels[pixelRoot];
            stripeToSend.sizeOfCluster[iy] = mLocalWuf->getClusterSize(pixelRoot);
        }
        else
        {
            stripeToSend.globalClusterId[iy] = INVALID_VALUE;
            stripeToSend.sizeOfCluster[iy] = INVALID_VALUE;
        }
    }
}

//---------------------------------------------------------------------------
// sendLeftStripeFromEvenReceiveOnOdd
void SendLeftColumnStrategy::sendStripeFromEvenReceiveOnOdd(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const
{
    const int numOfSends = NUM_OF_SENDS;
    const int msgId[numOfSends] = { MSG_1, MSG_2, MSG_3 };
    if (0 == (mDecompositionInfo.myRank % 2))
    {
        sendLeftStripe(stripeToSend, msgId, numOfSends);
    }
    else
    {
        receiveLeftStripe(stripeToReceive, msgId, numOfSends);
    }
}

//---------------------------------------------------------------------------
// sendLeftStripeFromOddReceiveOnEven
void SendLeftColumnStrategy::sendStripeFromOddReceiveOnEven(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const
{
    const int numOfSends = NUM_OF_SENDS;
    const int msgId[numOfSends] = { MSG_1, MSG_2, MSG_3 };
    if (0 != (mDecompositionInfo.myRank % 2))
    {
        sendLeftStripe(stripeToSend, msgId, numOfSends);
    }
    else
    {
        receiveLeftStripe(stripeToReceive, msgId, numOfSends);
    }
}

//---------------------------------------------------------------------------
// saveReceivedStripeToRightStripe
void SendLeftColumnStrategy::saveReceivedStripe(const SPixelStripe & stripeToReceive, std::vector<Pixel> & globalPixels) const
{
    // Check whether we participated in receiving the stripe.
    const int procToReceiveFrom = getRightNeighborProcessor(); // Periodic BCs are taken into account.

    if ( isNeighborProcessorValid(procToReceiveFrom) )
    {
        const int rightStripeId = (mDecompositionInfo.domainWidth + 1) * mDecompositionInfo.domainHeight;
        for (std::ptrdiff_t iy = 0; iy < mDecompositionInfo.domainHeight; ++iy)
        {
            globalPixels[rightStripeId + iy].pixelValue = stripeToReceive.pixelValue[iy];
            globalPixels[rightStripeId + iy].globalClusterId = stripeToReceive.globalClusterId[iy];
            globalPixels[rightStripeId + iy].sizeOfCluster = stripeToReceive.sizeOfCluster[iy];
        }
    }
}

//---------------------------------------------------------------------------
void SendLeftColumnStrategy::sendLeftStripe(SPixelStripe & stripeToSend, const int msgId[], const int size) const
{
    const int procToSendTo = getLeftNeighborProcessor(); // Periodic BCs are taken into account.
    if ( isNeighborProcessorValid(procToSendTo) )
    {
        MPI_Send(&stripeToSend.pixelValue[0], mDecompositionInfo.domainHeight, MPI_INT, procToSendTo, msgId[0], MPI_COMM_WORLD);
        MPI_Send(&stripeToSend.globalClusterId[0], mDecompositionInfo.domainHeight, MPI_INT, procToSendTo, msgId[1], MPI_COMM_WORLD);
        MPI_Send(&stripeToSend.sizeOfCluster[0], mDecompositionInfo.domainHeight, MPI_INT, procToSendTo, msgId[2], MPI_COMM_WORLD);
    }
}

//---------------------------------------------------------------------------
void SendLeftColumnStrategy::receiveLeftStripe(SPixelStripe & stripeToReceive, const int msgId[], const int size) const
{
    MPI_Status mpiStatus = {0};
    
    const int procToReceiveFrom = getRightNeighborProcessor(); // Periodic BCs are taken into account.
    if ( isNeighborProcessorValid(procToReceiveFrom) )
    {
        MPI_Recv(&stripeToReceive.pixelValue[0], mDecompositionInfo.domainHeight, MPI_INT,
                 procToReceiveFrom, msgId[0], MPI_COMM_WORLD, &mpiStatus);
        MPI_Recv(&stripeToReceive.globalClusterId[0], mDecompositionInfo.domainHeight, MPI_INT,
                 procToReceiveFrom, msgId[1], MPI_COMM_WORLD, &mpiStatus);
        MPI_Recv(&stripeToReceive.sizeOfCluster[0], mDecompositionInfo.domainHeight, MPI_INT,
                 procToReceiveFrom, msgId[2], MPI_COMM_WORLD, &mpiStatus);
    }
}
