// ParallelUnionFind2DStripes.cpp - implementation of the ParallelUnionFind2DStripes class
#include "ParallelUnionFind2DStripes.h"
#include "SendLeftColumnStrategy.h"
#include "SendRightColumnStrategy.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <mpi.h>
#include <iterator>

//---------------------------------------------------------------------------
ParallelUnionFind2DStripes::ParallelUnionFind2DStripes(const DecompositionInfo& info)
    : ParallelUnionFindImpl(info)
    , mNumOfPixels(info.domainWidth * info.domainHeight)
    , mNumOfGlobalPixels((info.domainWidth + 2) * info.domainHeight)
    , mLocalWuf(new WeightedUnionFind(mNumOfPixels))
    , mGlobalWuf()
    , mMerge()
    , mAllMerges()
    , mFinalWuf()
    , mTotalNumOfClusters(0)
    , mMinClusterSize(1)
    , mMaxClusterSize(1)
    , mClusterSizes()
    , mPercolatesHorizontally(false)
    , mPercolatesVertically(false)
    , mHorizPercolatedSize(0)
    , mVertPercolatedSize(0)
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
                    if (isPixelValid(neighbX))
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
    // TODO: put this to a separate method.
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
    if (BOSS != mDecompositionInfo.myRank)
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
    // In the latter case the pixels of the boundary stripes contain -1 and don't contribute to clusters.
    const int numOfExtendedPixels = mDecompositionInfo.domainHeight * (mDecompositionInfo.domainWidth + 2);

    mGlobalWuf = std::tr1::shared_ptr<WeightedUnionFind>(new WeightedUnionFind(numOfExtendedPixels));

    // At first run the local UF, but merge the pixels based on their global cluster id, not the pixel value.
    runLocalUfOnGlobalLabelsToSetInitialRoots();

#ifdef _DEBUG
    printGlobalUfRootsAfterFirstMerge();
#endif

    runUfOnGlobalPixelsAndRecordGlobalMerges();

#ifdef _DEBUG
    printMerges();
#endif
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::runLocalUfOnGlobalLabelsToSetInitialRoots()
{
    // We run a local UF but on all the pixels, including the received stripe(s).
    // But clusters are identified based on the global cluster id, and not on the pixel value.
    // In this way we will prepare the UF tree of merged clusters.

    const int numOfExtendedPixels = mDecompositionInfo.domainHeight * (mDecompositionInfo.domainWidth + 2);
    mGlobalWuf->reset(numOfExtendedPixels);              // Clear the UF. Necessary if we reuse the same UF.

    const int nx = mDecompositionInfo.domainWidth + 2;   // Two additional columns.
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
            int idp = indexTo1D(ix, iy);           // Convert 2D pixel coordinates into 1D index.
            if (mDecompositionInfo.pixelValue == mGlobalPixels[idp].pixelValue)
            {
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
}

//---------------------------------------------------------------------------
bool ParallelUnionFind2DStripes::isPixelValid(const int pixel) const
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

    // Replay all the unions from the global union list.
    getUniqueLabelForEachComponent();

    // I skip the two 2 final steps of the algorithm:
    // i) Transform the global labels of islands to the final ones using the final UF.
    // ii) Transform the final labeling so that the labels range from 0 to Nc-1, Nc - the total number of connected components.
    // Reasons for skipping:
    // 1) we can perform the transformation i) on the fly using mFinalWuf->getPixelRoot(globalRoot).
    // 2) transform ii) is not required for our goals.

    // Get the desired statistics.

    // Define min/max cluster sizes.
    getMinMaxClusterSizes();

    lookForPercolation();

    // Size histogram is computed in the corresponding 'print' method.
}

//---------------------------------------------------------------------------
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

        // Having sent/received the number of merges, we can send/receive the actual info about them.
        if (numOfMerges > 0)
        {
            broadcastMergeAndAddToAllMerges(mMerge.p, mAllMerges.p, numOfMerges, root);
            broadcastMergeAndAddToAllMerges(mMerge.pClusterSize, mAllMerges.pClusterSize, numOfMerges, root);
            broadcastMergeAndAddToAllMerges(mMerge.q, mAllMerges.q, numOfMerges, root);
            broadcastMergeAndAddToAllMerges(mMerge.qClusterSize, mAllMerges.qClusterSize, numOfMerges, root);
        }
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
    // Get the number of unmerged (across processors) clusters (some of them will be merged further).
    const int numOfGlobalLabels = calculateNumberOfGlobalLabels();

    // Create a UF with the 'numberOfGlobalLabels' entries.
    mFinalWuf = std::tr1::shared_ptr<WeightedUnionFind>(new WeightedUnionFind(numOfGlobalLabels));

    // Set the default roots for all the entries (necessary to get the consecutive labeling in the end).
    for (int index = 0; index < numOfGlobalLabels; ++index)
    {
        mFinalWuf->setInitialRoot(index, 0); // Initial cluster size is 0, not 1 as in the non-overloaded version.
    }

    // Replay all the unions that happened accross the processors' boundaries.
    const std::size_t numOfAllMerges = mAllMerges.p.size();
    for (std::size_t index = 0u; index < numOfAllMerges; ++index)
    {
        // Set initial roots/cluster sizes for the entries that are to be merged.
        const int idp = mAllMerges.p[index];
        const int idpClusterSize = mAllMerges.pClusterSize[index];
        mFinalWuf->setInitialRoot(idp, idpClusterSize);

        const int idq = mAllMerges.q[index];
        const int idqClusterSize = mAllMerges.qClusterSize[index];
        mFinalWuf->setInitialRoot(idq, idqClusterSize);

        // Connect the globla labels.
        if (!mFinalWuf->connected(idp, idq))  // Merge vertices only if they're not merged yet.
        {
            mFinalWuf->makeUnion(idp, idq);
        }
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::getMinMaxClusterSizes()
{
    // Define total number of clusters.
    // Map: first - nonconsecutive final label; second - consecutive final label.
    const std::map<int, int>& consecutiveFinalIds = mFinalWuf->getConsecutiveRootIds();
    mTotalNumOfClusters = consecutiveFinalIds.size();

    // Loop over all the (local) roots residing on the current proc.
    const std::map<int, int>& consecutiveLocalIds = mLocalWuf->getConsecutiveRootIds();
    std::map<int, int>::const_iterator iter = consecutiveLocalIds.begin();

    // Important: use the final UF for initialization.
    int minClusterSize = mDecompositionInfo.domainHeight * mDecompositionInfo.domainWidth;
    int maxClusterSize = INVALID_VALUE;

    for (; iter != consecutiveLocalIds.end(); ++iter)
    {
        const int localRoot = iter->first;
        const int globalRoot = mGlobalLabels[iter->first];
        // Map the global root id to the final id. For local clusters finalRoot == globalRoot.

        // Getting the final root using the global root (and not the pixel id as the function assumes)
        // is correct in this case, because mFinalWuf is formed based on the root ids and not pixel ids.
        const int finalRoot = mFinalWuf->getPixelRoot(globalRoot); // TODO: recheck, this might be wrong.

        // Assume at first the cluster spans several procs, hence, use mFinalWuf and final cluster size.
        int clusterSize = mFinalWuf->getClusterSize(finalRoot);

        // If cluster size is 0 then it is a (global label of a) local cluster, hence use the local cluster size.
        if (clusterSize <=0)
        {
            clusterSize = mLocalWuf->getClusterSize(localRoot);
        }

        minClusterSize = std::min(clusterSize, minClusterSize);
        maxClusterSize = std::max(clusterSize, maxClusterSize);

        // Save the cluster size and its final root (will be used for the size histogram).
        // Important: finalRoot should be the key as it is unique, clusterSize can have the same value for different clusters.
        mClusterSizes.insert( std::pair<int, int>(finalRoot, clusterSize) );
    } // end for root

    // Get the final min/max values.
    MPI_Reduce(&minClusterSize, &mMinClusterSize, 1, MPI_INT, MPI_MIN, BOSS, MPI_COMM_WORLD);
    MPI_Reduce(&maxClusterSize, &mMaxClusterSize, 1, MPI_INT, MPI_MAX, BOSS, MPI_COMM_WORLD);

    // Broadcast the values to others for calculating the cluster size histogram.
    MPI_Bcast(&mMinClusterSize, 1, MPI_INT, BOSS, MPI_COMM_WORLD);
    MPI_Bcast(&mMaxClusterSize, 1, MPI_INT, BOSS, MPI_COMM_WORLD);
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printPerProcessorClusterSizes(const std::string& fileName) const
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
void ParallelUnionFind2DStripes::printPerProcessorClusterStatistics(const std::string& fileName) const
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
void ParallelUnionFind2DStripes::printPerProcessorClusterSizeHistogram(const int bins, const std::string& fileName) const
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

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printClusterSizes(const std::string& fileName) const
{
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printClusterStatistics(const std::string& fileName) const
{
    if (BOSS == mDecompositionInfo.myRank)
    {
        std::cout << "Processor " << mDecompositionInfo.myRank << " report:" << std::endl;
        std::cout << "Found clusters: " << mTotalNumOfClusters << std::endl;
        std::cout << "Min cluster: " << mMinClusterSize << std::endl;
        std::cout << "Max cluster: " << mMaxClusterSize << std::endl;

        printPercolationInfo();
    }
}

//---------------------------------------------------------------------------
// TODO: introduce the version for 1 processor.
void ParallelUnionFind2DStripes::printClusterSizeHistogram(const int bins, const std::string& fileName) const
{
    if (bins > 1)
    {
        const double binWidth = static_cast<double>(mMaxClusterSize - mMinClusterSize)/static_cast<double>(bins - 1);

        // Avoid division by 0.
        if (binWidth < std::numeric_limits<double>::epsilon())
        {
            std::cerr << " Warning: processor " << mDecompositionInfo.myRank
                      << " binWidth is 0! Using binWidth == 1.0 to avoid division by 0." << std::endl;
        }
        
        // Calculate local size histogram for the current processor. Save the corresponding cluster root.
        std::vector<double> localHistogram(bins);
        std::multimap<int, int> rootsInBin; // Bin id is a key (duplicates are assumed), cluster root is value.

        std::map<int, int>::const_iterator iter;
        for (iter = mClusterSizes.begin(); iter != mClusterSizes.end(); ++iter)
        {
            int iChannel = static_cast<int>( ( (*iter).second - mMinClusterSize)*1./binWidth + 0.5 );
            rootsInBin.insert( std::pair<int, int>(iChannel, (*iter).first) );
            ++localHistogram[iChannel];
        }

        // Calculate the incorrect global histogram.
        // It is incorrect because those clusters that span several processors 
        // are taken into account for several times (depending on the number of procs they span).
        std::vector<double> finalHistogram(bins);
        MPI_Reduce(&localHistogram[0], &finalHistogram[0], bins, MPI_DOUBLE, MPI_SUM, BOSS, MPI_COMM_WORLD);

        // Ensure that clusters that span several processors are taken into account only once.
        adjustFinalHistogram(bins, finalHistogram, rootsInBin);

        outputSizeHistogram(bins, binWidth, fileName, finalHistogram);
    }
    else 
    {
        if (BOSS == mDecompositionInfo.myRank)
        {
            std::cerr << "Error: we need at least 2 bins for the cluster size histogram. The histogram has not been computed." << std::endl;
        }
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::adjustFinalHistogram(const int bins,
                                                      std::vector<double>& finalHistogram,
                                                      std::multimap<int, int>& rootsInBin) const
{
    // TODO: split this function into 2 functions.
    // Adjust the values in every bin.
    // Non-bosses pack the cluster roots to the vector and send them to the boss.
    typedef std::multimap<int, int>::iterator mapIter;
    if (BOSS != mDecompositionInfo.myRank)
    {
        // Put the roots sequentially in the array, separate bins with the INVALID_VALUE.
        std::vector<int> roots;
        int binId = 0;
        while (binId < bins)
        {
            // Get the range that contains the roots of the specified bin.
            std::pair<mapIter, mapIter> range = rootsInBin.equal_range(binId);

            // Pack the roots into the array.
            for (mapIter rangeIter = range.first; rangeIter != range.second; ++rangeIter)
            {
                roots.push_back( (*rangeIter).second );
            }

            // Separate the bin with the INVALID_VALUE divider.
            roots.push_back(INVALID_VALUE);

            ++binId;
        }

        // TODO: put sending/receiving a vector of the unknown size to a separate method.
        // Send the packed values to the BOSS.
        int numOfElements = roots.size();
        MPI_Send(&numOfElements, 1, MPI_INT, BOSS, MSG_1, MPI_COMM_WORLD);

        MPI_Send(&roots[0], numOfElements, MPI_INT, BOSS, MSG_1, MPI_COMM_WORLD);
    }
    else
    {
        // Boss receives the cluster roots, unpacks them and uses to adjust the histogram.
        for (int procId = 1; procId < mDecompositionInfo.numOfProc; ++procId)
        {
            // TODO: put sending/receiving a vector of the unknown size to a separate method.
            MPI_Status mpiStatus;
            int numOfElements = 0;
            // Receive the number of roots.
            MPI_Recv(&numOfElements, 1, MPI_INT, procId, MSG_1, MPI_COMM_WORLD, &mpiStatus);

            // Allocate the buffer.
            std::vector<int> buffer(numOfElements);
            // Receive the data.
            MPI_Recv(&buffer[0], numOfElements, MPI_INT, procId, MSG_1, MPI_COMM_WORLD, &mpiStatus);

            // Process the data: eliminate those clusters from the histogram that has been counted more than once.
            int binId = 0;
            int index = 0;
label:
            while (index < numOfElements)
            {
                // Skip all the empty bins.
                while (INVALID_VALUE == buffer[index])
                {
                    ++binId;
                    ++index;
                    goto label; // Use goto to be sure that other loops are not affected.
                }

                // Get roots that are known to the BOSS in the current bin.
                std::pair<mapIter, mapIter> range = rootsInBin.equal_range(binId);

                // Needed to prevent looping over the newly added element (this would happen if rootsInBin was used instead).
                std::multimap<int, int> bufferMap;

                // Loop over the received roots in the current bin and see whether BOSS knows about them.
                while (INVALID_VALUE != buffer[index])
                {
                    // If the rootsInBins bin is empty then insert a new bin.
                    if (range.first == range.second)
                    {
                        bufferMap.insert(std::pair<int, int> (binId, buffer[index]));
                    }

                    for (mapIter rangeIter = range.first; rangeIter != range.second; ++rangeIter)
                    {
                        if (buffer[index] == (*rangeIter).second)
                        {
                            --finalHistogram[binId];
                        }
                        else
                        {
                            bufferMap.insert(std::pair<int, int> (binId, buffer[index]));
                        }
                    }
                    ++index;
                }
                ++index; // Go to the 1st element of the next bin.
                ++binId;

                // Copy the elements to the original map.
                for (mapIter iter = bufferMap.begin(); iter != bufferMap.end(); ++iter)
                {
                    rootsInBin.insert(std::pair<int, int> ( (*iter).first, (*iter).second ));
                }
            } // End while (index < numOfElements).
        } // End for (int procId.
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::outputSizeHistogram(const int bins,
                                                     const double binWidth,
                                                     const std::string& fileName,
                                                     const std::vector<double>& finalHistogram) const
{
    // The final correct histogram is on the BOSS. Print it to a file.
    if (BOSS == mDecompositionInfo.myRank)
    {
        std::ofstream histFile(fileName);

        if (histFile.good() && (mTotalNumOfClusters > 0))
        {
            histFile.setf(std::ios::fixed,std::ios::floatfield); //histFile.width(10);
            double histSum = 0.;
            for (int i = 0; i < bins; ++i)
            {
                histSum += static_cast<double>(finalHistogram[i])/mTotalNumOfClusters;
                histFile << i*binWidth + mMinClusterSize << "\t";
                // Normalize only by the total number of islands (without delta).

                // TODO: in the end normalize.
                //histFile << static_cast<double>(finalHistogram[i])/mTotalNumOfClusters << "\t" << histSum << std::endl;
                histFile << static_cast<double>(finalHistogram[i]) << "\t" << histSum << std::endl;
            }
            histFile.close();
        }
        else
        {
            std::cerr << "Error: Bad histogram file." << std::endl;
        }
    }
}

//---------------------------------------------------------------------------
// TODO: introduce the version for 1 processor.
void ParallelUnionFind2DStripes::lookForPercolation()
{
    lookForHorizontalPercolation();
    lookForVerticalPercolation();
}

//---------------------------------------------------------------------------
// TODO: introduce the version for 1 processor.
void ParallelUnionFind2DStripes::lookForHorizontalPercolation()
{
    // TODO: split this function into 2 functions.
    if (BOSS == mDecompositionInfo.myRank)
    {
        // BOSS loops over the leftmost stripe of pixels and gets their final roots.
        std::set<int> leftMostVericalPixelRoots; // Use set to eliminate duplicates.
        for (std::size_t iy = 0u; iy < mDecompositionInfo.domainHeight; ++iy)
        {
            // TODO: recheck, getting final root using the globalRoot (not the pixel id) might be wrong.
            const int finalRoot = getFinalRootOfPixel(iy);
            if (INVALID_VALUE != finalRoot)
            {
                leftMostVericalPixelRoots.insert(finalRoot);
            }
        }

        // Receive the right-most roots from the last processor.
        int numOfReceivedRoots = 0;
        std::vector<int> rightMostVerticalPixelRoots;
        receiveData(rightMostVerticalPixelRoots, numOfReceivedRoots, mDecompositionInfo.numOfProc - 1);

        // Loop over the received roots and see whether one of them is in the BOSS roots.
        for (int id = 0; id < numOfReceivedRoots; ++id)
        {
            const int currentRoot = rightMostVerticalPixelRoots[id];
            if ( leftMostVericalPixelRoots.find(currentRoot) != leftMostVericalPixelRoots.end() )
            {
                mPercolatesHorizontally = true;
                mHorizPercolatedSize = mFinalWuf->getClusterSize(currentRoot);
                break;
            }
        }
    }
    else if (mDecompositionInfo.numOfProc - 1 == mDecompositionInfo.myRank)
    {
        // The last processor loops over the rightmost pixel stripe and gets its roots.
        std::set<int> rightMostVerticalPixelRoots;
        const int lastStripeStart = (mDecompositionInfo.domainWidth - 1)*mDecompositionInfo.domainHeight;
        for (std::size_t iy = 0; iy < mDecompositionInfo.domainHeight; ++iy)
        {
            const int finalRoot = getFinalRootOfPixel(iy + lastStripeStart);
            if (INVALID_VALUE != finalRoot)
            {
                rightMostVerticalPixelRoots.insert(finalRoot);
            }
        }

        // Pack the roots for sending.
        packDataAndSendIt(rightMostVerticalPixelRoots, BOSS);
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::lookForVerticalPercolation()
{
    // Every proc gets roots of the top-/bottom-most pixel stripe.
    std::set<int> topmostRoots;
    std::set<int> bottommostRoots;
    for (int ix = 0; ix < mDecompositionInfo.domainWidth; ++ix)
    {
        const int topRoot = getFinalRootOfPixel(ix * mDecompositionInfo.domainHeight);
        if (INVALID_VALUE != topRoot)
        {
            topmostRoots.insert(topRoot);
        }

        const int bottomRoot = getFinalRootOfPixel(mDecompositionInfo.domainHeight * (ix + 1) - 1);
        if (INVALID_VALUE != bottomRoot)
        {
            bottommostRoots.insert(bottomRoot);
        }
    }

    // Non-bosses send these data to the BOSS.
    if (BOSS != mDecompositionInfo.myRank)
    {
        packDataAndSendIt(topmostRoots, BOSS);
        packDataAndSendIt(bottommostRoots, BOSS);
    }
    else
    {
        // BOSS receives the roots.
        for (int procId = 1; procId < mDecompositionInfo.numOfProc; ++procId)
        {
            // Receive and save the topmost roots.
            int numOfReceivedRoots = 0;
            std::vector<int> dataToReceive;
            receiveData(dataToReceive, numOfReceivedRoots, procId);
            // Save the data to the BOSS's set.
            if (numOfReceivedRoots > 0)
            {
                std::copy( dataToReceive.begin(), dataToReceive.end(), std::inserter( topmostRoots, topmostRoots.end() ) );

            }

            // Receive and save the bottommost roots.
            receiveData(dataToReceive, numOfReceivedRoots, procId);
            if (numOfReceivedRoots > 0)
            {
                std::copy( dataToReceive.begin(), dataToReceive.end(), std::inserter( bottommostRoots, bottommostRoots.end() ) );
            }
        }

        // BOSS searches for at least one coincidence of the top and bottom roots.
        std::set<int>::iterator iter = topmostRoots.begin();
        for (; iter != topmostRoots.end(); ++iter)
        {
            const int currentRoot = *iter;
            if ( bottommostRoots.find(currentRoot) != bottommostRoots.end() )
            {
                mPercolatesVertically = true;
                mVertPercolatedSize = mFinalWuf->getClusterSize(currentRoot);
                break;
            }
        }
    }
}

//---------------------------------------------------------------------------
int ParallelUnionFind2DStripes::getFinalRootOfPixel(const int pixelId)
{
    int finalRoot = INVALID_VALUE;

    const int localRoot = mLocalWuf->getPixelRoot(pixelId);
    if (mLocalWuf->getClusterSize(localRoot) > 0) // Look only at the pixels of interest.
    {
        const int globalRoot = mGlobalLabels[localRoot];
        finalRoot = mFinalWuf->getPixelRoot(globalRoot);
    }

    return finalRoot;
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printPercolationInfo() const
{
    if (BOSS == mDecompositionInfo.myRank)
    {
        std::string h = " horizontal ";
        std::string v = " vertical ";

        std::string contact = " contact ";
        if (0 == mDecompositionInfo.pixelValue)
        {
            contact = " non-contact ";
        }

        if (mPercolatesHorizontally)
        {
            printPercolationPhrase(contact, h, mHorizPercolatedSize);
        }

        if (mPercolatesVertically)
        {
            printPercolationPhrase(contact, v, mVertPercolatedSize);
        }
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printPercolationPhrase(const std::string& contact, const std::string& vh, const int size) const
{
    std::cout << "# There is at least one" << contact << "percolated cluster in" << vh << "direction" << std:: endl;
    std::cout << "Cluster size is " << size << std::endl;
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::packDataAndSendIt(const std::set<int>& data, const int destinationProc) const
{
    std::vector<int> dataToSend;
    int numOfRoots = 0;

    std::set<int>::iterator iter = data.begin();
    for (; iter != data.end(); ++iter)
    {
        dataToSend.push_back( *iter );
        ++numOfRoots;
    }

    MPI_Send(&numOfRoots, 1, MPI_INT, destinationProc, MSG_1, MPI_COMM_WORLD);

    if (numOfRoots > 0)
    {
        MPI_Send(&dataToSend[0], numOfRoots, MPI_INT, destinationProc, MSG_1, MPI_COMM_WORLD);
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::receiveData(std::vector<int>& data, int& numOfElements, const int sendingProc) const
{
    MPI_Status mpiStatus;
    numOfElements = 0;
    MPI_Recv(&numOfElements, 1, MPI_INT, sendingProc, MSG_1, MPI_COMM_WORLD, &mpiStatus);

    if (numOfElements > 0)
    {
        data.resize(numOfElements);
        MPI_Recv(&data[0], numOfElements, MPI_INT, sendingProc, MSG_1, MPI_COMM_WORLD, &mpiStatus);
    }
}
