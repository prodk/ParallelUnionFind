// SendRightColumnStrategy.cpp - the implementation of the SendRightColumnStrategy class.

#include "SendRightColumnStrategy.h"
#include "WeightedUnionFind.h"
#include "ParallelUnionFindImpl.h"
#include <mpi.h>

//---------------------------------------------------------------------------
SendRightColumnStrategy::SendRightColumnStrategy(const DecompositionInfo & decompositionInfo,
                                                 const std::vector<int> & localPixels,
                                                 const std::vector<Pixel> & globalPixels,
                                                 const std::tr1::shared_ptr<WeightedUnionFind> & localWuf,
                                                 std::map<int, int> & globalLabels)
    : SendColumnStrategy(decompositionInfo, localPixels, globalPixels, localWuf, globalLabels)
{
}

//---------------------------------------------------------------------------
SendRightColumnStrategy::~SendRightColumnStrategy(void)
{
}

//---------------------------------------------------------------------------
void SendRightColumnStrategy::copyPixelStripeToSend(SPixelStripe & stripeToSend) const
{
    const int lastStripeStart = (mDecompositionInfo.domainWidth - 1)*mDecompositionInfo.domainHeight;
    for (std::size_t iy = 0; iy < mDecompositionInfo.domainHeight; ++iy)
    {
        const int lastStripeId = iy + lastStripeStart;
        stripeToSend.pixelValue[iy] = mLocalPixels[lastStripeId];

        // Set attributes only of those pixels that have the desired value.
        if (mDecompositionInfo.pixelValue == stripeToSend.pixelValue[iy])
        {
            const int pixelRoot = mLocalWuf->getPixelRoot(lastStripeId);
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
void SendRightColumnStrategy::sendStripeFromEvenReceiveOnOdd(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const
{
    const int numOfSends = NUM_OF_SENDS;
    const int msgId[numOfSends] = { MSG_1, MSG_2, MSG_3 };
    if (0 == (mDecompositionInfo.myRank % 2))
    {
        sendRightStripe(stripeToSend, msgId, numOfSends);
    }
    else
    {
        receiveRightStripe(stripeToReceive, msgId, numOfSends);
    }
}

//---------------------------------------------------------------------------
void SendRightColumnStrategy::sendStripeFromOddReceiveOnEven(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const
{
    const int numOfSends = NUM_OF_SENDS;
    const int msgId[numOfSends] = { MSG_1, MSG_2, MSG_3 };
    if (0 != (mDecompositionInfo.myRank % 2))
    {
        sendRightStripe(stripeToSend, msgId, numOfSends);
    }
    else
    {
        receiveRightStripe(stripeToReceive, msgId, numOfSends);
    }
}

//---------------------------------------------------------------------------
void SendRightColumnStrategy::saveReceivedStripe(const SPixelStripe & stripeToReceive, std::vector<Pixel> & globalPixels) const
{
    // Check whether we participated in receiving the stripe.
    const int procToReceiveFrom = getLeftNeighborProcessor(); // Periodic BCs are taken into account.
    if ( isNeighborProcessorValid(procToReceiveFrom) )
    {
        for (std::size_t iy = 0; iy < mDecompositionInfo.domainHeight; ++iy)
        {
            globalPixels[iy].pixelValue = stripeToReceive.pixelValue[iy];
            globalPixels[iy].globalClusterId = stripeToReceive.globalClusterId[iy];
            globalPixels[iy].sizeOfCluster = stripeToReceive.sizeOfCluster[iy];
        }
    }
}

//---------------------------------------------------------------------------
void SendRightColumnStrategy::sendRightStripe(SPixelStripe & stripeToSend, const int msgId[], const int size) const
{
    const int procToSendTo = getRightNeighborProcessor(); // Periodic BCs are taken into account.
    if ( isNeighborProcessorValid(procToSendTo) )
    {
        MPI_Send(&stripeToSend.pixelValue[0], mDecompositionInfo.domainHeight, MPI_INT, procToSendTo, msgId[0], MPI_COMM_WORLD);
        MPI_Send(&stripeToSend.globalClusterId[0], mDecompositionInfo.domainHeight, MPI_INT, procToSendTo, msgId[1], MPI_COMM_WORLD);
        MPI_Send(&stripeToSend.sizeOfCluster[0], mDecompositionInfo.domainHeight, MPI_INT, procToSendTo, msgId[2], MPI_COMM_WORLD);
    }
}

//---------------------------------------------------------------------------
void SendRightColumnStrategy::receiveRightStripe(SPixelStripe & stripeToReceive, const int msgId[], const int size) const
{
    MPI_Status mpiStatus = {0};
    
    const int procToReceiveFrom = getLeftNeighborProcessor(); // Periodic BCs are taken into account.
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
