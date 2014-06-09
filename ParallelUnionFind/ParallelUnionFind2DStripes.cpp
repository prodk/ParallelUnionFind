// ParallelUnionFind2DStripes.cpp - implementation of the ParallelUnionFind2DStripes class
#include "ParallelUnionFind2DStripes.h"

//---------------------------------------------------------------------------
ParallelUnionFind2DStripes::ParallelUnionFind2DStripes(const DecompositionInfo& info) :
                                ParallelUnionFindImpl(info),
                                mNumOfPixels(info.domainWidth * info.domainHeight),
                                mNumOfGlobalPixels((info.domainWidth + 2) * info.domainHeight),
                                mLocalWuf(new WeightedUnionFind(mNumOfPixels))
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
    mPixels.resize(mNumOfPixels);

    if (0 != mDecompositionInfo.pixels)
    {
        for (std::size_t i = 0u; i < mNumOfPixels; ++i)
        {
            mPixels[i] = mDecompositionInfo.pixels[i];
        }
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::runLocalUnionFind(void)
{
    if ((0 != mDecompositionInfo.pixels) && (mNumOfPixels > 0))
    {
        mLocalWuf->reset(mNumOfPixels);                // Clear the UF. Necessary if we reuse the same UF.

        const int nx = mDecompositionInfo.domainWidth;
        const int ny = mDecompositionInfo.domainHeight;

        for (int ix = 0; ix < nx; ++ix)                // Loop through the pixels.
        {
            for (int iy = 0; iy < ny; ++iy)
            {
                // TODO: reconsider periodic boundaries here!
                int neighbX = (ix + 1) % nx;     // Right neighbor with periodic boundaries.
                int neighbY = (iy + 1) % ny;     // Bottom neighbor with periodic boundaries.

                // Act only if the current pixel contains value.
                if (mDecompositionInfo.pixelValue == mPixels[indexTo1D(ix, iy)])
                {
                    int idp = indexTo1D(ix, iy); // Convert 2D pixel coordinates into 1D index.
                    mLocalWuf->setInitialRoot(idp);   // Set the root and the tree size (if it was 0).

                    // See whether neighboring (in both directions) pixels should be merged.
                    mergePixels(indexTo1D(neighbX, iy), idp); // Right neighbor.
                    mergePixels(indexTo1D(ix, neighbY), idp); // Bottom neighbor.
                }
            } // End for iy.
        } // End for ix.
    } // End if.
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::constructGlobalLabeling(void)
{
    // Get local consecutive labels of the clusters.
    const std::map<int, int>& consecutiveLocalIds = mLocalWuf->getConsecutiveRootIds();

    //----------------
    // TODO: put this part of code into a separate function.

    // Convert them to global labels.
    int myOffset = 0;
    // At first receive the number of clusters located on the processors with smaller ids.
    int msgId = 1;
    MPI_Status mpiStatus;
    if (0 != mDecompositionInfo.myRank) // Root doesn't receive anything, its offset is 0. The root initiates sending.
    {
        MPI_Recv(&myOffset, 1, MPI_INT, mDecompositionInfo.myRank-1, msgId, MPI_COMM_WORLD, &mpiStatus);
    }

    // Then send to the following proc the offset that includes the # of clusters on the current processor.
    // Offset for the next neighboring proc (i.e. the number of clusters to add on the next processor).
    const int numOfMyClusters = consecutiveLocalIds.size();
    int offsetForTheNextProcessor = myOffset + numOfMyClusters;

    if (mDecompositionInfo.myRank < mDecompositionInfo.numOfProc - 1) // Exclude the last processor from sending.
    {
        MPI_Send(&offsetForTheNextProcessor, 1, MPI_INT, mDecompositionInfo.myRank+1, msgId, MPI_COMM_WORLD);
    }

    // Construct global ids.
    std::map<int, int>::const_iterator iter;
    for (iter = consecutiveLocalIds.begin(); iter != consecutiveLocalIds.end(); ++iter)
    {
        // Non-consecutive local root is a key, a consecutive global root is a value.
        mGlobalLabels[iter->first] = iter->second + myOffset;
    }

    // TODO END.
    //----------------

#ifdef _DEBUG
    // Print ids.
    std::cout << "Loc \t LCsc \t Global" << std::endl;
    for (iter = consecutiveLocalIds.begin(); iter != consecutiveLocalIds.end(); ++iter)
    {
        std::cout << iter->first << "\t" << iter->second << "\t" << mGlobalLabels[iter->first] << std::endl;
    }
#endif
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::mergeLabelsAcrossProcessors(void)
{
    //----------------
    //TODO: split this function into several smaller ones.

    // An array to set the data of the globalWuf. It contains 2 additional columns of data.
    mGlobalPixels.resize(mNumOfGlobalPixels);

    // Copy the data from the localWuf to the array. Leave the first and the last columns empty.
    if (0 != mDecompositionInfo.pixels)
    {
        for (std::size_t ix = 0u; ix < mDecompositionInfo.domainWidth; ++ix)
        {
            for (std::size_t iy = 0u; iy < mDecompositionInfo.domainHeight; ++iy)
            {
                const int pixelGlobalId = indexTo1D(ix + 1, iy);
                const int pixelLocalId = indexTo1D(ix, iy);

                mGlobalPixels[pixelGlobalId].pixelValue = mPixels[pixelLocalId];

                const int pixelRoot = mLocalWuf->getPixelRoot(pixelLocalId);
                mGlobalPixels[pixelGlobalId].globalClusterId = mGlobalLabels[pixelRoot];
                mGlobalPixels[pixelGlobalId].sizeOfCluster = mLocalWuf->getClusterSize(pixelRoot);
            }
        }
    }

    // Copy the Merge data of the left stripe and send it to the left neighbor.
    std::vector<int> stripeToSend(mDecompositionInfo.domainHeight);
    std::vector<int> stripeToReceive(mDecompositionInfo.domainHeight);
    for (std::size_t iy = 0; iy < mDecompositionInfo.domainHeight; ++iy)
    {
        stripeToSend[iy] = mPixels[iy];

       /* const int pixelRoot = mLocalWuf->getPixelRoot(iy);
        stripe[iy].globalClusterId = mGlobalLabels[pixelRoot];
        stripe[iy].sizeOfCluster = mLocalWuf->getClusterSize(pixelRoot);*/
    }
    // Send from even procs, receive from odd
    const int msgId = 123;
    if (0 == (mDecompositionInfo.myRank % 2))
    {
        // TODO: reconsider periodic boundaries!
        if (mDecompositionInfo.myRank > 0)
        {
            const int procToSendTo = mDecompositionInfo.myRank - 1;
            MPI_Send(&stripeToSend[0], mDecompositionInfo.domainHeight, MPI_INT, procToSendTo, msgId, MPI_COMM_WORLD);
        }
    }
    else
    {
        if (mDecompositionInfo.myRank != mDecompositionInfo.numOfProc - 1)
        {
            MPI_Status mpiStatus;
            // TODO: reconsider periodic boundaries!
            const int procToReceiveFrom = mDecompositionInfo.myRank + 1;
            MPI_Recv(&stripeToReceive[0], mDecompositionInfo.domainHeight, MPI_INT, procToReceiveFrom, msgId, MPI_COMM_WORLD, &mpiStatus);
        }
    }

    // Then send form odd, receive from even
    if (0 != (mDecompositionInfo.myRank % 2))
    {
        // TODO: reconsider periodic boundaries!
        const int procToSendTo = (mDecompositionInfo.myRank - 1);
        MPI_Send(&stripeToSend[0], mDecompositionInfo.domainHeight, MPI_INT, procToSendTo, msgId, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Status mpiStatus;
        // TODO: reconsider periodic boundaries!
        const int procToReceiveFrom = (mDecompositionInfo.myRank + 1);
        MPI_Recv(&stripeToReceive[0], mDecompositionInfo.domainHeight, MPI_INT, procToReceiveFrom, msgId, MPI_COMM_WORLD, &mpiStatus);
    }

    // Save the values of the right stripe, received from the left neighbor.
    for (std::size_t iy = 0; iy < mDecompositionInfo.domainHeight; ++iy)
    {
        const int rightStripeId = (mDecompositionInfo.domainWidth + 1) * mDecompositionInfo.domainHeight;
        mGlobalPixels[rightStripeId + iy].pixelValue = stripeToReceive[iy];

       /* const int pixelRoot = mLocalWuf->getPixelRoot(iy);
        stripe[iy].globalClusterId = mGlobalLabels[pixelRoot];
        stripe[iy].sizeOfCluster = mLocalWuf->getClusterSize(pixelRoot);*/
    }

#ifdef _DEBUG
    printLocalExtendedPicture(mDecompositionInfo);
#endif

    // Receive the stripe from the left neighbor.

    // Copy the Merge data of the right stripe and send it to the right neighbor

    // Receive the stripe from the right neighbor.

    // Initialize the global UF with the data.

    // TODO END.
    //----------------
}

//---------------------------------------------------------------------------
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
