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

    // TODO: add find method that returns the root of the component a point belongs to.
    bool Connected(int p, int q);// True if both verteces are in the same cluster.
    void Union(int p, int q);    // Put the verteces to the same cluster.
    void SetInitialRoot(int id); // Set the initial root of the vertex.
    void Reset(int N);           // Set the union find to a state as if it had been newly created.

    // Testing.
    bool CompareResults(int resultMinClusterSize, int resultMaxClusterSize, int numOfClusters);

    // I/O.
    void PrintId();             // Prints ids of all the vertices.
    void PrintSize();           // Prints sizes of clusters that correspond to tree vertices.
    void PrintClusters();       // Cluster statistics.
    int PrintClusterSizeHistogram(const int bins, 
        const std::string &fileOut);  // Output cluster histogram to a file.

private:
    std::vector<int> id;        // Ids of the verteces.
    std::vector<int> size;      // Number of vertices in each tree (or cluster).
    std::set<int> roots;        // A set of ids of the roots of the trees.
    int minClusterSize;
    int maxClusterSize;

    int Root(int i);            // Find the root of the vertex i.
    void GetMinMaxClusterSize(int *min, int *max);  // Define clusters with the minimum and maximum sizes.
};

#endif // WEIGHTED_UNION_FIND