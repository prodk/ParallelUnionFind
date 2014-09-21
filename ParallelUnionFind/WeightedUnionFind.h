// WeightedUnionFind.h - WeightedUnionFind class declaration.

#ifndef WEIGHTED_UNION_FIND
#define WEIGHTED_UNION_FIND

//---------------------------------------------------------------------------
#include <vector>
#include <set>
#include <map>
#include <cstddef>    // For ptrdiff_t.
                      // ptrdiff_t is necessary for huge systems. It works even when the integer value exceeds INT_MAX = 2^32.
                      // This might happen for the 64-bit JUQUEEN and huge systems.
#include <iostream>
#include <fstream>
#include <string>
#include <limits>     // std::numeric_limits<> for gcc

//---------------------------------------------------------------------------
struct Pixel
{
    int pixelValue;
    std::ptrdiff_t globalClusterId;
    std::ptrdiff_t sizeOfCluster;
};

//---------------------------------------------------------------------------
namespace Common
{
    // Non-member non-friend function to add one vector to the end of another.
    template <typename T>
    void addOneVectorToAnother (std::vector<T> & mainVector, std::vector<T> & toAdd)
    {
        mainVector.insert(mainVector.end(), toAdd.begin(), toAdd.end());
    }
}

//---------------------------------------------------------------------------
class WeightedUnionFind
{
public:
    WeightedUnionFind(const ptrdiff_t N);
    ~WeightedUnionFind(void);

    // Interface for the client.
    bool connected(std::ptrdiff_t p, std::ptrdiff_t q);        // True if both vertices are in the same cluster.
    void makeUnion(std::ptrdiff_t p, std::ptrdiff_t q);        // Put the vertices to the same cluster.
    void setInitialRoot(std::ptrdiff_t id);                    // Set the initial root of the vertex.
    void setInitialRoot(std::ptrdiff_t id, std::ptrdiff_t clusterSize); // Overloaded version.
    void reset(std::ptrdiff_t N);                              // Set the union find to a state as if it had been newly created.

    // Change local labeling of the clusters such that it is consecutive.
    const std::map<std::ptrdiff_t, std::ptrdiff_t>& getConsecutiveRootIds();

    // Output the results.
    void printClusterSizes(std::ostream& out);                 // Print sizes of clusters that correspond to tree vertices.
    void printClusterStatistics(std::ostream& out);            // Cluster statistics.

    void setPixelRoot(std::ptrdiff_t pixelId, std::ptrdiff_t rootId);
    void setClusterSize(std::ptrdiff_t rootId, std::ptrdiff_t clusterSize);

    std::ptrdiff_t getPixelRoot(const std::ptrdiff_t id) const
    {
        return root(id);
    }

    std::ptrdiff_t getClusterSize(const std::ptrdiff_t rootId) const
    {
        return mSize[rootId];
    }

    // Output cluster histogram to a file.
    int printClusterSizeHistogram(const int bins, const std::string& fileOut);

private:
    std::ptrdiff_t root(std::ptrdiff_t i) const;                                     // Find the root of the vertex i.
    void getMinMaxClusterSize(std::ptrdiff_t *min, std::ptrdiff_t *max);  // Define the minimum and maximum cluster sizes.
    void buildAndPrintSizeHistogram(const int bins, const std::string& fileOut) const;

private:
    static const int defaultInt;
    std::vector<std::ptrdiff_t> mId;                            // Ids of the vertices.
    std::vector<std::ptrdiff_t> mSize;                          // Number of vertices in each tree (i. e. cluster).
    std::set<std::ptrdiff_t> mRoots;                            // A set of ids of the roots of the trees.
    std::map<std::ptrdiff_t, std::ptrdiff_t> mConsecutiveRoots; // Relabeled cluster indices.
    std::ptrdiff_t mMinClusterSize;
    std::ptrdiff_t mMaxClusterSize;
};

#endif // WEIGHTED_UNION_FIND