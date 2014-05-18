// ParallelUnionFind2DStripes.cpp - implementation of the ParallelUnionFind2DStripes class
#include "ParallelUnionFind2DStripes.h"

//---------------------------------------------------------------------------
ParallelUnionFind2DStripes::ParallelUnionFind2DStripes(const DecompositionInfo& info) :
                                ParallelUnionFindImpl(info),
                                mNumOfPixels(info.domainWidth * info.domainHeight),
                                mWuf(new WeightedUnionFind(mNumOfPixels))
{
    if(mDecompositionInfo.numOfProc <= 0)
    {
        std::cerr << "0 number of processors!" << std::endl;
        std::cerr << "Check whether the MPI has been initialized!" << std::endl;
    }
    else
    {
        std::cout << "Processor " << mDecompositionInfo.myRank << " of " << mDecompositionInfo.numOfProc << std::endl;
        std::cout << "PUF of type 2DStripes created, width " 
            << mDecompositionInfo.domainWidth << " height " << mDecompositionInfo.domainHeight << std::endl;

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
        mWuf->reset(mNumOfPixels);                // Clear the UF. Necessary if we reuse a Pixels object.

        const std::size_t nx = mDecompositionInfo.domainWidth;
        const std::size_t ny = mDecompositionInfo.domainHeight;

        for(int ix = 0; ix < nx; ++ix)           // Loop through the pixels.
        {
            for(int iy = 0; iy < ny; ++iy)
            {
                int neighbX = (ix + 1) % nx;     // Right neighbor with periodic boundaries.
                int neighbY = (iy + 1) % ny;     // Bottom neighbor with periodic boundaries.

                // Act only if the current pixel contains value.
                if(mDecompositionInfo.pixelValue == mPixels[indexTo1D(ix, iy)]){
                    int idp = indexTo1D(ix, iy); // Convert 2D pixel coordinates into 1D index.
                    mWuf->setInitialRoot(idp);   // Set the root and the tree size (if it was 0).

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
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::mergeLabelsAcrossProcessors(void)
{
}

//---------------------------------------------------------------------------
void ParallelUnionFind2DStripes::performFinalLabelingOfClusters(void)
{
}
