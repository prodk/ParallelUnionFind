// ParallelUnionFindImpl.h - declaration of the ParallelUnionFindImpl class 
// which is a part of the Bridge pattern.

// Using the Template Method pattern this class defines the general sceleton
// of the parallel union-find algorithm described in:
// C.Harrison et al. Eurographics Symposium on Parallel Graphics and Visualization, 2011.

#ifndef PARALLEL_UNION_FIND_IMPL
#define PARALLEL_UNION_FIND_IMPL

//---------------------------------------------------------------------------
#include <string>
#include <vector>
#include <cstddef> // For ptrdiff_t.
                   // ptrdiff_t is necessary for huge systems. It works even when the integer value exceeds INT_MAX = 2^32.

//---------------------------------------------------------------------------
class WeightedUnionFind;

//---------------------------------------------------------------------------
struct DecompositionInfo
{
    DecompositionInfo() :
    numOfProc(0), myRank(-1), domainWidth(1), domainHeight(1), pixels(0), pixelValue(1), periodicBoundaryX(false) {}

    int numOfProc;
    int myRank;
    std::ptrdiff_t domainWidth;  // Width of the domain per processor, in grid points.
    std::ptrdiff_t domainHeight; // Height of the domain per processor, in grid points.
    const int* pixels;           // A pointer to the array of grid points in 1D format.
    int pixelValue;              // Value that is used for uniting the points.
                                 // 1 corresponds to contact, 0 to non-contact.
    bool periodicBoundaryX;      // True if use pbc in the x direction. False by default.
};

//---------------------------------------------------------------------------
// Records a single merge.
struct Merge
{
    // We put the arrays into the structure to simplify using Merge structure in MPI
    // We do not want to bother with creating new MPI data type,
    // as it seems easier to just send arrays of MPI_INTs.
    std::vector<std::ptrdiff_t> p;                   // Root id of the 1st cluster.
    std::vector<std::ptrdiff_t> q;                   // Root id of the 2nd cluster.
    std::vector<std::ptrdiff_t> pClusterSize;        // Size of the cluster of the p root.
    std::vector<std::ptrdiff_t> qClusterSize;        // Size of the cluster of the q root.
};

//---------------------------------------------------------------------------
class ParallelUnionFindImpl
{
public:
    ParallelUnionFindImpl(const DecompositionInfo& info);
    virtual ~ParallelUnionFindImpl(void);

    void analyze(void);         // TM pattern: specify a skeleton for the algorithm.

    // Output functions, obligatory to override.
    virtual void printClusterSizes(const std::string& fileName) const = 0;
    virtual void printClusterStatistics(const std::string& fileName) const = 0;
    virtual void printClusterSizeHistogram(const int bins, const std::string& fileName) const = 0;

    // Output local (i.e. related to the current processors) quantities (usually for debugging purposes).
    virtual void printPerProcessorClusterSizes(const std::string& fileName) const = 0;
    virtual void printPerProcessorClusterStatistics(const std::string& fileName) const = 0;
    virtual void printPerProcessorClusterSizeHistogram(const int bins, const std::string& fileName) const = 0;

    // Helper interface functions.
    virtual void setPixelValue(const int value) = 0;

private:
    // These are the fixed steps of the algorithm called from the analyze().
    // Private functions can be overriden by children, but can't be explicitly invoked. They're called only in the base class.
    virtual void runLocalUnionFind(void) = 0;
    virtual void constructGlobalLabeling(void) = 0;
    virtual void mergeLabelsAcrossProcessors(void) = 0;
    virtual void performFinalLabelingOfClusters(void) = 0;

protected:
    DecompositionInfo mDecompositionInfo;
};

#endif // PARALLEL_UNION_FIND_IMPL