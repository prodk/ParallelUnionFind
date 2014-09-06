// WeightedUnionFind.h - WeightedUnionFind class declaration.

#ifndef WEIGHTED_UNION_FIND
#define WEIGHTED_UNION_FIND

//---------------------------------------------------------------------------
#include <vector>
#include <set>
#include <map>

//---------------------------------------------------------------------------
struct Pixel
{
    int pixelValue;
    int globalClusterId;
    int sizeOfCluster;
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
    WeightedUnionFind(const std::size_t N);
    ~WeightedUnionFind(void);

    // Interface for the client.
    bool connected(int p, int q);        // True if both vertices are in the same cluster.
    void makeUnion(int p, int q);        // Put the vertices to the same cluster.
    void setInitialRoot(int id);         // Set the initial root of the vertex.
    void reset(int N);                   // Set the union find to a state as if it had been newly created.

    // Change local labeling of the clusters such that it is consecutive.
    const std::map<int, int>& getConsecutiveRootIds();

    // Output the results.
    void printClusterSizes(std::ostream& out);     // Prints sizes of clusters that correspond to tree vertices.
    void printClusterStatistics(std::ostream& out);// Cluster statistics.

    void setPixelRoot(int pixelId, int rootId);
    void setClusterSize(int rootId, int clusterSize);

    int getPixelRoot(const int id) const
    {
        return root(id);
    }

    int getClusterSize(const int rootId) const
    {
        return mSize[rootId];
    }

    // Output cluster histogram to a file.
    int printClusterSizeHistogram(const int bins, const std::string& fileOut);

private:
    int root(int i) const;                          // Find the root of the vertex i.
    void getMinMaxClusterSize(int *min, int *max);  // Define clusters with the minimum and maximum sizes.
    void buildAndPrintSizeHistogram(const int bins, const std::string& fileOut) const;

private:
    static const int defaultInt;
    std::vector<int> mId;          // Ids of the vertices.
    std::vector<int> mSize;        // Number of vertices in each tree (i. e. cluster).
    std::set<int> mRoots;          // A set of ids of the roots of the trees.
    std::map<int, int> mConsecutiveRoots; // Relabeled cluster indices.
    int mMinClusterSize;
    int mMaxClusterSize;
};

#endif // WEIGHTED_UNION_FIND