// ParallelUnionFind.h - declaration of ParalleUnionFind class

#ifndef PARALLEL_UNION_FIND
#define PARALLEL_UNION_FIND

#include "ParallelUnionFindImpl.h"

class ParallelUnionFind
{
public:
    ~ParallelUnionFind(void);
    ParallelUnionFind(const std::string& spatialConfiguration, // How the domain decomposition is done.
                      const DecompositionInfo& info);

    // Interface the client must work with.
    void analyze(void);
    void printClusterStatistics(const std::string& fileName) const;
    void printClusterSizeHistogram(const int bins, const std::string& fileName) const;
    void setPixelValue(const int value);

    // Log messages in the debug mode.
#ifdef _DEBUG
private:
    void logMsg(std::ostream& out, const std::string& msg) const;
#endif

private:
    std::tr1::shared_ptr<ParallelUnionFindImpl> mParallelUf; // Use the Bridge to separate the interface from the impl.
};

#endif // PARALLEL_UNION_FIND
