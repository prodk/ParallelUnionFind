// ParallelUnionFind2DStripes.cpp - implementation of the ParallelUnionFind2DStripes class
#include "ParallelUnionFind2DStripes.h"

//---------------------------------------------------------------------------
ParallelUnionFind2DStripes::ParallelUnionFind2DStripes(const DecompositionInfo& info) :
                                ParallelUnionFindImpl(info),
                                mNumOfPixels(info.domainWidth * info.domainHeight),
                                mLocalWuf(new WeightedUnionFind(mNumOfPixels))
{
    if(mDecompositionInfo.numOfProc <= 0)
    {
        std::cerr << "0 number of processors!" << std::endl;
        std::cerr << "Check whether the MPI has been initialized!" << std::endl;
    }
    else
    {
        std::cout << "Processor " << mDecompositionInfo.myRank << " of " << mDecompositionInfo.numOfProc << std::endl;
        std::cout << "PUF of type 2DStripes created" << std::endl;
        std::cout << "width " << mDecompositionInfo.domainWidth << " height " << mDecompositionInfo.domainHeight << std::endl;

        std::cout << "Copying pixels to the internal array ..." << std::endl << std::endl;
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

    if(0 != mDecompositionInfo.pixels)
    {
        for(std::size_t i = 0u; i < mNumOfPixels; ++i)
        {
            mPixels[i] = mDecompositionInfo.pixels[i];
        }
    }
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::runLocalUnionFind(void)
{
    if((0 != mDecompositionInfo.pixels) && (mNumOfPixels > 0))
    {
        mLocalWuf->reset(mNumOfPixels);                // Clear the UF. Necessary if we run UF several times.

        const int nx = mDecompositionInfo.domainWidth;
        const int ny = mDecompositionInfo.domainHeight;

        for(int ix = 0; ix < nx; ++ix)            // Loop through the pixels.
        {
            for(int iy = 0; iy < ny; ++iy)
            {
                int neighbX = (ix + 1) % nx;     // Right neighbor with periodic boundaries.
                int neighbY = (iy + 1) % ny;     // Bottom neighbor with periodic boundaries.

                // Act only if the current pixel contains value.
                if(mDecompositionInfo.pixelValue == mPixels[indexTo1D(ix, iy)]){
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
    mGlobalLabels.resize(numOfMyClusters);
    std::map<int, int>::const_iterator iter;
    int count = 0;
    for (iter = consecutiveLocalIds.begin(); iter != consecutiveLocalIds.end(); ++iter)
    {
        mGlobalLabels[count] = iter->second + myOffset;
        ++count;
    }

#ifdef _DEBUG
    // Print ids.
    std::cout << " Loc \t LCsc \t Global" << std::endl;
    int i = 0;
    for (iter = consecutiveLocalIds.begin(); iter != consecutiveLocalIds.end(); ++iter)
    {
        std::cout << iter->first << "\t" << iter->second << "\t" << mGlobalLabels[i] << std::endl;
        ++i;
    }
#endif
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::mergeLabelsAcrossProcessors(void)
{
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::performFinalLabelingOfClusters(void)
{
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::printClusterSizes(const std::string& fileName) const
{
    if ( (fileName.length() > 0) && ("" != fileName) )
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
    if ( (fileName.length() > 0) && ("" != fileName) )
    {
        std::fstream fileStream(fileName);
        if (fileStream.good())
        {
            fileStream << std::endl << "Pixel value used: " << mDecompositionInfo.pixelValue << std::endl;
            mLocalWuf->printClusterStatistics(fileStream);
            fileStream.close();
        }
    }
    else
    {
        std::cout << std::endl << "Pixel value used: " << mDecompositionInfo.pixelValue << std::endl;
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
