// ParallelUnionFind2DStripes.cpp - implementation of the ParallelUnionFind2DStripes class
#include "ParallelUnionFind2DStripes.h"
#include "SendLeftColumnStrategy.h"
#include "SendRightColumnStrategy.h"
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
    , mGlobalWuf()
    , mMerge()
    , mAllMerges()
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
                int idp = indexTo1D(ix, iy);           // Convert 2D pixel coordinates into 1D index.
                if (mDecompositionInfo.pixelValue == mLocalPixels[idp])
                {
                    mLocalWuf->setInitialRoot(idp);    // Set the root and the tree size (if it was 0).

                    // See whether neighboring (in both directions) pixels should be merged.
                    const int neighbX = getNeighborNonPeriodicBC(ix, nx);   // Right neighbor without periodic boundaries.
                    if (isNeighborPixelValid(neighbX))
                    {
                        const int idx = indexTo1D(neighbX, iy);
                        mergePixels(idx, idp, mLocalWuf, mLocalPixels[idx]);
                    }

                    const int neighbY = getNeighborPeriodicBC(iy, ny);      // Bottom neighbor with periodic boundaries.
                    const int idy = indexTo1D(ix, neighbY);
                    mergePixels(idy, idp, mLocalWuf, mLocalPixels[idy]);
                }
            } // End for iy.
        } // End for ix.
    } // End if.
    else // We have no valid pixels.
    {
        if (BOSS == mDecompositionInfo.myRank)
        {
            std::cerr << "Error: there are no valid pixels to work with!" << std::endl;
        }
        const int errorCode = 123;
        MPI_Abort(MPI_COMM_WORLD, errorCode);
    }
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
    const int msgId = MSG_1;
    MPI_Status mpiStatus;

    // Root doesn't receive anything, its offset is 0. The root initiates sending.
    if (0 != mDecompositionInfo.myRank)
    {
        const int receiveFromProc = mDecompositionInfo.myRank - 1;
        MPI_Recv(&numOfClusters, 1, MPI_INT, receiveFromProc, msgId, MPI_COMM_WORLD, &mpiStatus);
    }

    return numOfClusters;
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::sendTotalClustersToNextProcs(const int numOfClustersOnSmallerProcIds, const int numOfMyClusters) const
{
    int offsetForTheNextProcessor = numOfClustersOnSmallerProcIds + numOfMyClusters;

    const int msgId = MSG_1;
    if (mDecompositionInfo.myRank < mDecompositionInfo.numOfProc - 1) // Exclude the last processor from sending.
    {
        const int sendToProc = mDecompositionInfo.myRank + 1;
        MPI_Send(&offsetForTheNextProcessor, 1, MPI_INT, sendToProc, msgId, MPI_COMM_WORLD);
    }
}

//---------------------------------------------------------------------------
// Stage 3.
void ParallelUnionFind2DStripes::mergeLabelsAcrossProcessors(void)
{
    // An array to set the data of the globalWuf. It contains 2 additional columns of data (per processor).
    mGlobalPixels.resize(mNumOfGlobalPixels);

    // Copy the data from the localWuf to the array.
    // Setup the global pixels. Note: first/last columns are not empty (i.e. they don't contain garbage).
    initializeGloblaPixels();

    copyLeftColumnAndSendToLeftNeighbor();

    copyRightColumnAndSendToRightNeighbor();

#ifdef _DEBUG
    printLocalExtendedPicture(mDecompositionInfo);
    printReceivedGlobalLabels();
#endif

    // Run UF on the global UF and record the merges that happen.
    runUfOnGlobalLabels();
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::initializeGloblaPixels(void)
{
    if (0 != mDecompositionInfo.pixels)
    {
        // Init everything to -1.
        const std::size_t numOfExtendedPixels = mGlobalPixels.size();
        for (std::size_t index = 0u; index < numOfExtendedPixels; ++index)
        {
            mGlobalPixels[index].pixelValue = INVALID_VALUE;
            mGlobalPixels[index].globalClusterId = INVALID_VALUE;
            mGlobalPixels[index].sizeOfCluster = INVALID_VALUE;
        }

        // Copy the data from the localWuf to the inner part of the global pixels.
        for (std::size_t ix = 0u; ix < mDecompositionInfo.domainWidth; ++ix)
        {
            for (std::size_t iy = 0u; iy < mDecompositionInfo.domainHeight; ++iy)
            {
                const int pixelGlobalId = indexTo1D(ix + 1, iy);
                const int pixelLocalId = indexTo1D(ix, iy);

                mGlobalPixels[pixelGlobalId].pixelValue = mLocalPixels[pixelLocalId];

                // Fill in only those pixels that have the desired value.
                if (mDecompositionInfo.pixelValue == mGlobalPixels[pixelGlobalId].pixelValue)
                {
                    const int pixelRoot = mLocalWuf->getPixelRoot(pixelLocalId);
                    mGlobalPixels[pixelGlobalId].globalClusterId = mGlobalLabels[pixelRoot];
                    mGlobalPixels[pixelGlobalId].sizeOfCluster = mLocalWuf->getClusterSize(pixelRoot);
                }
                else // TODO: perhaps remove this initialization if it is redundant.
                {
                    mGlobalPixels[pixelGlobalId].globalClusterId = INVALID_VALUE;
                    mGlobalPixels[pixelGlobalId].sizeOfCluster = INVALID_VALUE;
                }
            }
        }
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::copyLeftColumnAndSendToLeftNeighbor(void)
{
    SendLeftColumnStrategy leftColumn(mDecompositionInfo, mLocalPixels, mGlobalPixels, mLocalWuf, mGlobalLabels);
    leftColumn.sendReceivePixelStripes(mGlobalPixels);
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::copyRightColumnAndSendToRightNeighbor(void)
{
    SendRightColumnStrategy rightColumn(mDecompositionInfo, mLocalPixels, mGlobalPixels, mLocalWuf, mGlobalLabels);
    rightColumn.sendReceivePixelStripes(mGlobalPixels);
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::runUfOnGlobalLabels()
{
    // We use the same extended width even when there are no periodic BCs.
    // In the latter case the pixels of the boundary stripes contain -1 and are don't contribute to clusters.
    const int numOfExtendedPixels = mDecompositionInfo.domainHeight * (mDecompositionInfo.domainWidth + 2);

    mGlobalWuf = static_cast<std::tr1::shared_ptr<WeightedUnionFind> >(new WeightedUnionFind(numOfExtendedPixels));

    // At first run the local UF, but merge the pixels based on their global cluster id, not the pixel value.
    runLocalUfOnGlobalLabelsToSetInitialRoots();

#ifdef _DEBUG
    printGlobalUfRootsAfterFirstMerge();
#endif

    runUfOnGlobalPixelsAndRecordGlobalMerges();

    printMerges();
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::runLocalUfOnGlobalLabelsToSetInitialRoots()
{
    // We run a local UF but on all the pixels, including the received stripe(s).
    // But clusters are identified based on the global cluster id, and not on the pixel value.
    // In this way we will prepare the UF tree of merged clusters.

    const int numOfExtendedPixels = mDecompositionInfo.domainHeight * (mDecompositionInfo.domainWidth + 2);
    mGlobalWuf->reset(numOfExtendedPixels);                // Clear the UF. Necessary if we reuse the same UF.

    const int nx = mDecompositionInfo.domainWidth + 2;
    const int ny = mDecompositionInfo.domainHeight;

    for (int ix = 0; ix < nx; ++ix)                      // Loop through the pixels, columns fastest.
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            // Act only if the cluster id is valid.
            int idp = indexTo1D(ix, iy);                 // Convert 2D pixel coordinates into 1D index.
            if ( isClusterIdValid(idp) )
            {
                mGlobalWuf->setInitialRoot(idp);         // Set the root and the tree size (if it was 0).

                // See whether neighboring (in both directions) pixels should be merged.
                const int neighbX = getNeighborNonPeriodicBC(ix, nx);   // Right neighbor without periodic boundaries.
                // Merge pixels only if they belong to the same cluster.
                if (neighbX >= 0)
                {
                    const int idx = indexTo1D(neighbX, iy);
                    mergeClusterIds(idx, idp, mGlobalWuf);
                }

                const int neighbY = getNeighborPeriodicBC(iy, ny);      // Bottom neighbor with periodic boundaries.
                const int idy = indexTo1D(ix, neighbY);
                mergeClusterIds(idy, idp, mGlobalWuf);
            }
        } // End for iy.
    } // End for ix.
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::mergeClusterIds(int idq, int idp, std::tr1::shared_ptr<WeightedUnionFind> wuf) const
{
    if (mGlobalPixels[idq].globalClusterId == mGlobalPixels[idp].globalClusterId)
    {
        wuf->setInitialRoot(idq);       // Specify the non-zero size (if it was 0) and init root.
        if (!wuf->connected(idp, idq))  // Merge vertices only if they're not yet merged.
        {
            wuf->makeUnion(idp, idq);
        }
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::runUfOnGlobalPixelsAndRecordGlobalMerges()
{
    // Run the local UF again on the extended pixels.
    // But now merge clusters based on the pixel values and not on the global cluster's ids.
    // Record the merges that occur.
    const int numOfExtendedPixels = mDecompositionInfo.domainHeight * (mDecompositionInfo.domainWidth + 2);
    mGlobalWuf->reset(numOfExtendedPixels);                // Clear the UF. Necessary if we reuse the same UF.

    const int nx = mDecompositionInfo.domainWidth + 2;
    const int ny = mDecompositionInfo.domainHeight;

    for (int ix = 0; ix < nx; ++ix)                // Loop through the pixels, columns fastest.
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            // TODO: recheck this logic because still the same merges are recorded several times.
            // Act only if the cluster id is valid.
            int idp = indexTo1D(ix, iy);      // Convert 2D pixel coordinates into 1D index.
            if (mDecompositionInfo.pixelValue == mGlobalPixels[idp].pixelValue)
            {
                // TODO: perhaps delete this line, we assume that mGlobalWuf has alredy been run at least once.
                mGlobalWuf->setInitialRoot(idp);   // Set the root and the tree size (if it was 0).

                // See whether neighboring (in both directions) pixels should be merged.
                const int neighbX = getNeighborNonPeriodicBC(ix, nx);   // Right neighbor without periodic boundaries.
                // Merge pixels only if they belong to the same cluster.
                if (neighbX >= 0)
                {
                    const int idx = indexTo1D(neighbX, iy);
                    mergePixelsAndRecordMerge(idx, idp, mGlobalWuf, mGlobalPixels[idx].pixelValue);
                }

                const int neighbY = getNeighborPeriodicBC(iy, ny);      // Bottom neighbor with periodic boundaries.
                const int idy = indexTo1D(ix, neighbY);
                mergePixelsAndRecordMerge(idy, idp, mGlobalWuf, mGlobalPixels[idy].pixelValue);
            }
        } // End for iy.
    } // End for ix.
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::mergePixelsAndRecordMerge(int idq, int idp, std::tr1::shared_ptr<WeightedUnionFind> wuf,
                                                            const int pixelValue)
{
    if ( (mDecompositionInfo.pixelValue == pixelValue) && 
         (mGlobalPixels[idq].globalClusterId != mGlobalPixels[idp].globalClusterId) ) // Avoid redundant merges of the same cluster.
    {
        wuf->setInitialRoot(idq);       // Specify the non-zero size (if it was 0) and init root.
        if (!wuf->connected(idp, idq))  // Merge vertices only if they're not yet merged.
        {
            wuf->makeUnion(idp, idq);
            recordMerge(idp, idq);
        }
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::recordMerge(const int idp, const int idq)
{
    mMerge.p.push_back(mGlobalPixels[idp].globalClusterId);
    mMerge.q.push_back(mGlobalPixels[idq].globalClusterId);
    mMerge.pClusterSize.push_back(mGlobalPixels[idp].sizeOfCluster);
    mMerge.qClusterSize.push_back(mGlobalPixels[idq].sizeOfCluster);
    mMerge.clusterSize.push_back(mGlobalPixels[idp].sizeOfCluster + mGlobalPixels[idq].sizeOfCluster);
}

//---------------------------------------------------------------------------
bool ParallelUnionFind2DStripes::isNeighborPixelValid(const int pixel) const
{
    return (pixel != INVALID_VALUE);
}

//---------------------------------------------------------------------------
bool ParallelUnionFind2DStripes::isClusterIdValid(const int pixel) const
{
    return (mGlobalPixels[pixel].globalClusterId > INVALID_VALUE);
}

//---------------------------------------------------------------------------
// Stage 4.
void ParallelUnionFind2DStripes::performFinalLabelingOfClusters(void)
{
    // Broadcast each processor's merges that we've just recorded.
    getMergesFromAllProcs();

    // Replay all the unions from the global union list
    getUniqueLabelForEachComponent();

    // Transform the local labels of islands to the global ones using the global UF

    // Transform the final labeling so that the labels range from 0 to Nc-1
    // where Nc - the total number of connected components

}

void ParallelUnionFind2DStripes::getMergesFromAllProcs()
{
    for (int root = 0; root < mDecompositionInfo.numOfProc; ++root)
    {
        int numOfMerges = 0;
        if (mDecompositionInfo.myRank == root)
        {
            numOfMerges = static_cast<int>(mMerge.p.size()); // Number of merges to send for the current root.
        }

        // Send/receive the number of merges.
        MPI_Bcast(&numOfMerges, 1, MPI_INT, root, MPI_COMM_WORLD);

#ifdef _DEBUG
        std::cout << "proc " << root << " numOfMerges " << numOfMerges << std::endl;
#endif

        // Having sent/received the number of merges, we can send/receive the actual info about them
        broadcastMergeAndAddToAllMerges(mMerge.clusterSize, mAllMerges.clusterSize, numOfMerges, root);
        broadcastMergeAndAddToAllMerges(mMerge.p, mAllMerges.p, numOfMerges, root);
        broadcastMergeAndAddToAllMerges(mMerge.pClusterSize, mAllMerges.pClusterSize, numOfMerges, root);
        broadcastMergeAndAddToAllMerges(mMerge.q, mAllMerges.q, numOfMerges, root);
        broadcastMergeAndAddToAllMerges(mMerge.qClusterSize, mAllMerges.qClusterSize, numOfMerges, root);
    } // end for (root = 0;
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::broadcastMergeAndAddToAllMerges(const std::vector<int> & arrayToSend,
                                                                 std::vector<int> & arrayToReceive,
                                                                 int numOfMerges, int root)
{
    std::vector<int> buffer(numOfMerges); // Allocate a buffer where we put the received data.

    // Pack the data if we are sending it, otherwise set the buffer to 0s.
    if (mDecompositionInfo.myRank == root)
    {
        buffer = arrayToSend;
    }
    else
    {
        std::fill(buffer.begin(), buffer.end(), 0);
    }

    MPI_Bcast(&buffer[0], numOfMerges, MPI_INT, root, MPI_COMM_WORLD);
    Common::addOneVectorToAnother<int>(arrayToReceive, buffer);
}

//---------------------------------------------------------------------------
int ParallelUnionFind2DStripes::calculateNumberOfGlobalLabels() const
{
    int numberOfGlobalLabels = 0;

    int myNumberOfLables = mGlobalLabels.size();

    MPI_Reduce(&myNumberOfLables, &numberOfGlobalLabels, 1, MPI_INT, MPI_SUM, BOSS, MPI_COMM_WORLD);

    // Give the just calculated value to all the processors.
    MPI_Bcast(&numberOfGlobalLabels, 1, MPI_INT, BOSS, MPI_COMM_WORLD);

    return numberOfGlobalLabels;
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::getUniqueLabelForEachComponent()
{
    // Get the number of unmerged clusters (some of them will be merged further).
    const int numberOfGlobalLabels = calculateNumberOfGlobalLabels();
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

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printReceivedGlobalLabels() const
{
    std::stringstream fileName;
    fileName << "proc_" << mDecompositionInfo.myRank << "_mGLobal_labels.dat";
    std::ofstream outFile(fileName.str());

    if (outFile.good())
    {
        const int nx = mDecompositionInfo.domainWidth + 2;
        const int ny = mDecompositionInfo.domainHeight;

        for (int iy = 0; iy < ny; ++iy)            // Loop through the pixels, rows fastest.
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                outFile << mGlobalPixels[indexTo1D(ix, iy)].globalClusterId << " ";
            }
            outFile << std::endl;
        }
        outFile.close();
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printGlobalUfRootsAfterFirstMerge() const
{
    std::stringstream fileName;
    fileName << "proc_" << mDecompositionInfo.myRank << "_globalUf_roots.dat";
    std::ofstream outFile(fileName.str());

    if (outFile.good())
    {
        const int nx = mDecompositionInfo.domainWidth + 2;
        const int ny = mDecompositionInfo.domainHeight;

        for (int iy = 0; iy < ny; ++iy)            // Loop through the pixels, rows fastest.
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                outFile << mGlobalWuf->getPixelRoot(indexTo1D(ix, iy)) << " ";
            }
            outFile << std::endl;
        }
        outFile.close();
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printMerges() const
{
    const std::size_t numOfMerges = mMerge.p.size();

    for (std::size_t i = 0u; i < numOfMerges; ++i)
    {
        std::cout << "merge " << i << " p " << mMerge.p[i] << " q " << mMerge.q[i] << std::endl;
    }
}
