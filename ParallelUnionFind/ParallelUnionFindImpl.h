// ParallelUnionFindImpl.h - declaration of the ParallelUnionFindImpl class 
// which is a part of the Bridge pattern.
// Using the Template Method pattern this class defines the general sceleton
// of the parallel union-find algorithm described in:
// C.Harrison et al. Eurographics Symposium on Parallel Graphics and Visualization, 2011.

#ifndef PARALLEL_UNION_FIND_IMPL
#define PARALLEL_UNION_FIND_IMPL

//---------------------------------------------------------------------------
#include <string>

//---------------------------------------------------------------------------
class WeightedUnionFind;

//---------------------------------------------------------------------------
struct DecompositionInfo
{
    DecompositionInfo() :
    numOfProc(0), myRank(-1), domainWidth(1), domainHeight(1), pixels(0), pixelValue(1) {}

    int numOfProc;
    int myRank;
    std::size_t domainWidth; // Width of the domain per processor, in grid points.
    std::size_t domainHeight;// Height of the domain per processor, in grid points.
    const int* pixels;       // A pointer to the array of grid points in 1D format.
    int pixelValue;          // Value that is used for uniting the points.
                             // 1 corresponds to contact, 0 to non-contact.
};

//---------------------------------------------------------------------------
struct Pixel
{
    int pixelValue;
    int globalClusterId;
    int sizeOfCluster;
};

//---------------------------------------------------------------------------
// Records a single merge.
struct Merge
{
    int p;
    int q;
    int rootId;
    int clusterSize;
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