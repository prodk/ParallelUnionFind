// ParallelUnionFind2DStripes.cpp - implementation of the ParallelUnionFind2DStripes class
#include "ParallelUnionFind2DStripes.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <mpi.h>

//---------------------------------------------------------------------------
ParallelUnionFind2DStripes::ParallelUnionFind2DStripes(const DecompositionInfo& info)
    : ParallelUnionFindImpl(info)
    , mNumOfPixels(info.domainWidth * info.domainHeight)
    , mNumOfGlobalPixels((info.domainWidth + 2) * info.domainHeight)
    , mLocalWuf(new WeightedUnionFind(mNumOfPixels))
{
    if (mDecompositionInfo.numOfProc <= 0)
    {
        std::cerr << "0 number of processors!" << std::endl;
        std::cerr << "Check whether the MPI has been initialized!" << std::endl;
    }
    else
    {
        std::cout << "Processor " << mDecompositionInfo.myRank << " of " << mDecompositionInfo.numOfProc << ": PUF of type 2DStripes created.";
        std::cout << " Width " << mDecompositionInfo.domainWidth << " height " << mDecompositionInfo.domainHeight << "." << std::endl;

        copyPixels();
    }
}

//---------------------------------------------------------------------------
ParallelUnionFind2DStripes::~ParallelUnionFind2DStripes(void)
{
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::copyPixels()
{
    mLocalPixels.resize(mNumOfPixels);

    if (0 != mDecompositionInfo.pixels)
    {
        for (std::size_t i = 0u; i < mNumOfPixels; ++i)
        {
            mLocalPixels[i] = mDecompositionInfo.pixels[i];
        }
    }
}

//---------------------------------------------------------------------------
// Stage 1.
void ParallelUnionFind2DStripes::runLocalUnionFind(void)
{
    if ((0 != mDecompositionInfo.pixels) && (mNumOfPixels > 0))
    {
        mLocalWuf->reset(mNumOfPixels);                // Clear the UF. Necessary if we reuse the same UF.

        const int nx = mDecompositionInfo.domainWidth;
        const int ny = mDecompositionInfo.domainHeight;

        for (int ix = 0; ix < nx; ++ix)                // Loop through the pixels, columns fastest.
        {
            for (int iy = 0; iy < ny; ++iy)
            {
                // Act only if the current pixel contains the desired value.
                if (mDecompositionInfo.pixelValue == mLocalPixels[indexTo1D(ix, iy)])
                {
                    int idp = indexTo1D(ix, iy);      // Convert 2D pixel coordinates into 1D index.
                    mLocalWuf->setInitialRoot(idp);   // Set the root and the tree size (if it was 0).

                    // See whether neighboring (in both directions) pixels should be merged.
                    const int neighbX = getNeighborNonPeriodicBC(ix, nx);   // Right neighbor without periodic boundaries.
                    if (neighbX >= 0)
                    {
                        mergePixels(indexTo1D(neighbX, iy), idp);
                    }

                    const int neighbY = getNeighborPeriodicBC(iy, ny);      // Bottom neighbor with periodic boundaries.
                    mergePixels(indexTo1D(ix, neighbY), idp);
                }
            } // End for iy.
        } // End for ix.
    } // End if.
}

//---------------------------------------------------------------------------
// Stage 2.
void ParallelUnionFind2DStripes::constructGlobalLabeling(void)
{
    // Get local consecutive labels of the clusters.
    const std::map<int, int>& consecutiveLocalIds = mLocalWuf->getConsecutiveRootIds();

    // Get the offset from the procs with smaller ids.
    const int numOfClustersOnSmallerProcIds = receiveNumberOfClustersFromPreviousProcs();

    // Send to the following proc the offset that includes the # of clusters on the current processor.
    const int numOfMyClusters = consecutiveLocalIds.size();
    sendTotalClustersToNextProcs(numOfClustersOnSmallerProcIds, numOfMyClusters);

    // Construct global consecutive ids of the clusters.
    std::map<int, int>::const_iterator iter;
    for (iter = consecutiveLocalIds.begin(); iter != consecutiveLocalIds.end(); ++iter)
    {
        mGlobalLabels[iter->first] =                     // Non-consecutive local root is a key (i.e. iter->first).
            // iter->second is a local consecutive id of the cluster.
            iter->second + numOfClustersOnSmallerProcIds;// Consecutive global root is a value stored in the map.
    }

#ifdef _DEBUG
    // Print ids.
    std::cout << "Loc  LocConsec GlobalConsec" << std::endl;
    for (iter = consecutiveLocalIds.begin(); iter != consecutiveLocalIds.end(); ++iter)
    {
        std::cout << iter->first << "\t" << iter->second << "\t" << mGlobalLabels[iter->first] << std::endl;
    }
#endif
}

//---------------------------------------------------------------------------
int ParallelUnionFind2DStripes::receiveNumberOfClustersFromPreviousProcs() const
{
    int numOfClusters = 0;

    // Receive the number of clusters located on the processors with ids smaller than ours.
    const int msgId = 1;
    MPI_Status mpiStatus;

    // Root doesn't receive anything, its offset is 0. The root initiates sending.
    if (0 != mDecompositionInfo.myRank)
    {
        const int receiveFromProc = mDecompositionInfo.myRank-1;
        MPI_Recv(&numOfClusters, 1, MPI_INT, receiveFromProc, msgId, MPI_COMM_WORLD, &mpiStatus);
    }

    return numOfClusters;
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::sendTotalClustersToNextProcs(const int numOfClustersOnSmallerProcIds, const int numOfMyClusters) const
{
    int offsetForTheNextProcessor = numOfClustersOnSmallerProcIds + numOfMyClusters;

    const int msgId = 1;
    if (mDecompositionInfo.myRank < mDecompositionInfo.numOfProc - 1) // Exclude the last processor from sending.
    {
        const int sendToProc = mDecompositionInfo.myRank+1;
        MPI_Send(&offsetForTheNextProcessor, 1, MPI_INT, sendToProc, msgId, MPI_COMM_WORLD);
    }
}

//---------------------------------------------------------------------------
// Stage 3.
void ParallelUnionFind2DStripes::mergeLabelsAcrossProcessors(void)
{
    // An array to set the data of the globalWuf. It contains 2 additional columns of data.
    mGlobalPixels.resize(mNumOfGlobalPixels);

    // Copy the data from the localWuf to the array.
    // Setup the global pixels. Note: first/last columns are not empty.
    setLocalPartOfGloblaPixels();

    copyLeftColumnAndSendToLeftNeighbor();

    copyRightColumnAndSendToRightNeighbor();

#ifdef _DEBUG
    printLocalExtendedPicture(mDecompositionInfo);
#endif

    // Initialize the global UF with the local and received data.

    // Run UF on the global UF and record the merges that happen.

}

//---------------------------------------------------------------------------
// TODO: rename this function or factor out the code for the first/last stripes.
void ParallelUnionFind2DStripes::setLocalPartOfGloblaPixels(void)
{
    if (0 != mDecompositionInfo.pixels)
    {
         // Initialize the first and the last columns of the global pixels.
        const int lastStripeStart = (mDecompositionInfo.domainWidth + 1)*mDecompositionInfo.domainHeight;
        for (std::size_t iy = 0u; iy < mDecompositionInfo.domainHeight; ++iy)
        {
            const int firstStripeId = iy;
            // TODO: get rid of the magic -1 number: introduce an enum.
            mGlobalPixels[firstStripeId].pixelValue = -1;
            mGlobalPixels[firstStripeId].globalClusterId = -1;
            mGlobalPixels[firstStripeId].sizeOfCluster = -1;

            const int lastStripeId = iy + lastStripeStart;
            // TODO: get rid of the magic -1 number: introduce an enum.
            mGlobalPixels[lastStripeId].pixelValue = -1;
            mGlobalPixels[lastStripeId].globalClusterId = -1;
            mGlobalPixels[lastStripeId].sizeOfCluster = -1;
        }

        // Copy the data from the localWuf to the inner part of the global pixels.
        for (std::size_t ix = 0u; ix < mDecompositionInfo.domainWidth; ++ix)
        {
            for (std::size_t iy = 0u; iy < mDecompositionInfo.domainHeight; ++iy)
            {
                const int pixelGlobalId = indexTo1D(ix + 1, iy);
                const int pixelLocalId = indexTo1D(ix, iy);

                mGlobalPixels[pixelGlobalId].pixelValue = mLocalPixels[pixelLocalId];

                const int pixelRoot = mLocalWuf->getPixelRoot(pixelLocalId);
                mGlobalPixels[pixelGlobalId].globalClusterId = mGlobalLabels[pixelRoot];
                mGlobalPixels[pixelGlobalId].sizeOfCluster = mLocalWuf->getClusterSize(pixelRoot);
            }
        }
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::copyLeftColumnAndSendToLeftNeighbor(void)
{
    SPixelStripe stripeToSend(mDecompositionInfo.domainHeight);
    copyLeftPixelStripeToSend(stripeToSend);

    SPixelStripe stripeToReceive(mDecompositionInfo.domainHeight);
    sendLeftStripeFromEvenReceiveOnOdd(stripeToSend, stripeToReceive);
    sendLeftStripeFromOddReceiveOnEven(stripeToSend, stripeToReceive);

    saveReceivedStripeToRightStripe(stripeToReceive);
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::copyLeftPixelStripeToSend(SPixelStripe & stripeToSend)
{
    for (std::size_t iy = 0; iy < mDecompositionInfo.domainHeight; ++iy)
    {
        stripeToSend.pixelValue[iy] = mLocalPixels[iy];

        const int pixelRoot = mLocalWuf->getPixelRoot(iy);
        stripeToSend.globalClusterId[iy] = mGlobalLabels[pixelRoot];
        stripeToSend.sizeOfCluster[iy] = mLocalWuf->getClusterSize(pixelRoot);
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::sendLeftStripeFromEvenReceiveOnOdd( SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const
{
    const int numOfSends = 3;                 // TODO: get rid of magic numbers, perhaps use sizeof(stripeToSend)/sizeof(stripeToSend.clusterId).
    const int msgId[numOfSends] = { 123, 456, 789 };   // TODO: get rid of magic numbers.
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
void ParallelUnionFind2DStripes::sendLeftStripeFromOddReceiveOnEven(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const
{
    const int numOfSends = 3;                 // TODO: get rid of magic numbers, perhaps use sizeof(stripeToSend)/sizeof(stripeToSend.clusterId).
    const int msgId[numOfSends] = { 123, 456, 789 };   // TODO: get rid of magic numbers.
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
void ParallelUnionFind2DStripes::sendLeftStripe(SPixelStripe & stripeToSend, const int msgId[], const int size) const
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
void ParallelUnionFind2DStripes::receiveLeftStripe(SPixelStripe & stripeToReceive, const int msgId[], const int size) const
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

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::saveReceivedStripeToRightStripe(const SPixelStripe & stripeToReceive)
{
    // Check whether we participated in receiving the stripe.
    const int procToReceiveFrom = getRightNeighborProcessor(); // Periodic BCs are taken into account.

    if ( isNeighborProcessorValid(procToReceiveFrom) )
    {
        const int rightStripeId = (mDecompositionInfo.domainWidth + 1) * mDecompositionInfo.domainHeight;
        for (std::size_t iy = 0; iy < mDecompositionInfo.domainHeight; ++iy)
        {
            mGlobalPixels[rightStripeId + iy].pixelValue = stripeToReceive.pixelValue[iy];
            mGlobalPixels[rightStripeId + iy].globalClusterId = stripeToReceive.globalClusterId[iy];
            mGlobalPixels[rightStripeId + iy].sizeOfCluster = stripeToReceive.sizeOfCluster[iy];
        }
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::copyRightColumnAndSendToRightNeighbor(void)
{
    SPixelStripe stripeToSend(mDecompositionInfo.domainHeight);
    copyRightPixelStripeToSend(stripeToSend);
    
    SPixelStripe stripeToReceive(mDecompositionInfo.domainHeight);
    sendRightStripeFromEvenReceiveOnOdd(stripeToSend, stripeToReceive);
    sendRightStripeFromOddReceiveOnEven(stripeToSend, stripeToReceive);

    saveReceivedStripeToLeftStripe(stripeToReceive);

}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::copyRightPixelStripeToSend(SPixelStripe & stripeToSend)
{
    const int lastStripeStart = (mDecompositionInfo.domainWidth - 1)*mDecompositionInfo.domainHeight;
    for (std::size_t iy = 0; iy < mDecompositionInfo.domainHeight; ++iy)
    {
        const int lastStripeId = iy + lastStripeStart;
        stripeToSend.pixelValue[iy] = mLocalPixels[lastStripeId];

        const int pixelRoot = mLocalWuf->getPixelRoot(lastStripeId);
        stripeToSend.globalClusterId[iy] = mGlobalLabels[pixelRoot];
        stripeToSend.sizeOfCluster[iy] = mLocalWuf->getClusterSize(pixelRoot);
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::sendRightStripeFromEvenReceiveOnOdd( SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const
{
    const int numOfSends = 3;                 // TODO: get rid of magic numbers, perhaps use sizeof(stripeToSend)/sizeof(stripeToSend.clusterId).
    const int msgId[numOfSends] = { 123, 456, 789 };   // TODO: get rid of magic numbers.
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
void ParallelUnionFind2DStripes::sendRightStripeFromOddReceiveOnEven(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const
{
    const int numOfSends = 3;                 // TODO: get rid of magic numbers, perhaps use sizeof(stripeToSend)/sizeof(stripeToSend.clusterId).
    const int msgId[numOfSends] = { 123, 456, 789 };   // TODO: get rid of magic numbers.
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
void ParallelUnionFind2DStripes::sendRightStripe(SPixelStripe & stripeToSend, const int msgId[], const int size) const
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
void ParallelUnionFind2DStripes::receiveRightStripe(SPixelStripe & stripeToReceive, const int msgId[], const int size) const
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

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::saveReceivedStripeToLeftStripe(const SPixelStripe & stripeToReceive)
{
    // Check whether we participated in receiving the stripe.
    const int procToReceiveFrom = getLeftNeighborProcessor(); // Periodic BCs are taken into account.
    if ( isNeighborProcessorValid(procToReceiveFrom) )
    {
        for (std::size_t iy = 0; iy < mDecompositionInfo.domainHeight; ++iy)
        {
            mGlobalPixels[iy].pixelValue = stripeToReceive.pixelValue[iy];
            mGlobalPixels[iy].globalClusterId = stripeToReceive.globalClusterId[iy];
            mGlobalPixels[iy].sizeOfCluster = stripeToReceive.sizeOfCluster[iy];
        }
    }
}

//---------------------------------------------------------------------------
int ParallelUnionFind2DStripes::getLeftNeighborProcessor() const
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
int ParallelUnionFind2DStripes::getRightNeighborProcessor() const
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
bool ParallelUnionFind2DStripes::isNeighborProcessorValid(const int rank) const
{
    return (rank >= 0);
}

//---------------------------------------------------------------------------
// Stage 4.
void ParallelUnionFind2DStripes::performFinalLabelingOfClusters(void)
{
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printClusterSizes(const std::string& fileName) const
{
    if ((fileName.length() > 0) && ("" != fileName))
    {
        std::fstream fileStream(fileName);
        if (fileStream.good())
        {
            mLocalWuf->printClusterSizes(fileStream);
            fileStream.close();
        }
    }
    else
    {
        mLocalWuf->printClusterSizes(std::cout);
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printClusterStatistics(const std::string& fileName) const
{
    if ((fileName.length() > 0) && ("" != fileName))
    {
        std::fstream fileStream(fileName);
        if (fileStream.good())
        {
            fileStream << "Pixel value used: " << mDecompositionInfo.pixelValue << std::endl;
            mLocalWuf->printClusterStatistics(fileStream);
            fileStream.close();
        }
    }
    else
    {
        std::cout << "Pixel value used: " << mDecompositionInfo.pixelValue << std::endl;
        mLocalWuf->printClusterStatistics(std::cout);
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printClusterSizeHistogram(const int bins, const std::string& fileName) const
{
    mLocalWuf->printClusterSizeHistogram(bins, fileName);
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::setPixelValue(const int value)
{
    mDecompositionInfo.pixelValue = value;
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printLocalExtendedPicture(const DecompositionInfo& info) const
{
    std::stringstream fileName;
    fileName << "proc_" << info.myRank << "_mGLobal_picture.dat";
    std::ofstream outFile(fileName.str());

    if (outFile.good())
    {
        const int nx = info.domainWidth + 2;
        const int ny = info.domainHeight;

        for (int iy = 0; iy < ny; ++iy)            // Loop through the pixels, rows fastest.
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                outFile << mGlobalPixels[indexTo1D(ix, iy)].pixelValue;
            }
            outFile << std::endl;
        }
        outFile.close();
    }
}
