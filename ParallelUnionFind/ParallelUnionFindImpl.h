// ParallelUnionFindImpl.h - declaration of the ParallelUnionFindImpl class which is a part of the Bridge pattern.
// Through the Template Method pattern defines the general sceleton
// of the parallel union-find algorithm described in:
// C.Harrison et al. Eurographics Symposium on Parallel Graphics and Visualization, 2011.

#ifndef PARALLEL_UNION_FIND_IMPL
#define PARALLEL_UNION_FIND_IMPL

#include <mpi.h>
#include <string>
#include "WeightedUnionFind.h"

struct DecompositionInfo
{
    DecompositionInfo() : 
    numOfProc(0), myRank(-1), domainWidth(1), domainHeight(1), pixels(0), indexFactor(1), pixelValue(1) {}

    int numOfProc;
    int myRank;
    std::size_t domainWidth; // Width of the domain per processor, in grid points.
    std::size_t domainHeight;// Height of the domain per processor, in grid points.
    const int* pixels;       // A pointer to the array of grid points in 1D format.
    int indexFactor;         // A factor to convert 2D indeces into 1D form. 
                             // It is necessary, because in the parallel FFTW it is 2, not 1.
    int pixelValue;          // Value that is used for uniting the points.
                             // For the contact area we use 1, which corresponds to contact.
};

class ParallelUnionFindImpl
{
public:
    ParallelUnionFindImpl(const DecompositionInfo& info);
    virtual ~ParallelUnionFindImpl(void);

    virtual void analyze(void);  // Template method pattern: specify a skeleton for the algorithm.

    // Output functions, obligatory to override.
    virtual void printClusterSizes(const std::string& fileName) const = 0;
    virtual void printClusterStatistics(const std::string& fileName) const = 0;
    virtual void printClusterSizeHistogram(const int bins, const std::string& fileName) const = 0;

    // Helper interface functions.
    virtual void setPixelValue(const int value) = 0;

private:
    // These are the fixed steps of the algorithm called from the analyze().
    // Private functions can be overriden by children, but can't be invoked. They're invoked only in the base class.
    virtual void runLocalUnionFind(void) = 0;
    virtual void constructGlobalLabeling(void) = 0;
    virtual void mergeLabelsAcrossProcessors(void) = 0;
    virtual void performFinalLabelingOfClusters(void) = 0;

protected:
    DecompositionInfo mDecompositionInfo;
};

#endif // PARALLEL_UNION_FIND_IMPL
