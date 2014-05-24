// WeightedUnionFind.h - WeightedUnionFind class declaration.

#ifndef WEIGHTED_UNION_FIND
#define WEIGHTED_UNION_FIND

#include <iostream>
#include <vector>
#include <set>
#include <fstream>

class WeightedUnionFind
{
public:
    WeightedUnionFind(const std::size_t N);
    ~WeightedUnionFind(void);

    // Interface the client is restricted to use.
    bool connected(int p, int q);        // True if both verteces are in the same cluster.
    void makeUnion(int p, int q);        // Put the verteces to the same cluster.
    void setInitialRoot(int id);         // Set the initial root of the vertex.
    void reset(int N);                   // Set the union find to a state as if it had been newly created.

    // Output the results.
    //void printId(std::ostream out);      // Prints ids of all the vertices.
    void printClusterSizes(std::ostream& out);     // Prints sizes of clusters that correspond to tree vertices.
    void printClusterStatistics(std::ostream& out);// Cluster statistics.

    // Output cluster histogram to a file.
    int printClusterSizeHistogram(const int bins, const std::string& fileOut);

private:
    int root(int i);               // Find the root of the vertex i.
    void getMinMaxClusterSize(int *min, int *max);  // Define clusters with the minimum and maximum sizes.
    void buldAndPrintSizeHistogram(const int bins, const std::string& fileOut) const;

private:
    std::vector<int> mId;          // Ids of the verteces.
    std::vector<int> mSize;        // Number of vertices in each tree (i. e. cluster).
    std::set<int> mRoots;          // A set of ids of the roots of the trees.
    int mMinClusterSize;
    int mMaxClusterSize;
};

#endif // WEIGHTED_UNION_FIND