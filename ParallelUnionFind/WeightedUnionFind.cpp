// WeightedUnionFind.cpp - WeightedUnionFind class implementation.
#include "WeightedUnionFind.h"

WeightedUnionFind::WeightedUnionFind(int N)
{
    id.resize(N);
    size.resize(N);

    for(int i = 0; i < N; ++i)
    {
        id[i] = i;		// Assign the id to the sequence number of the vertex.
        size[i] = 0;	// No elements in all the trees at the beginning.
    }

    minClusterSize = 0;
    maxClusterSize = 0;
}

WeightedUnionFind::~WeightedUnionFind(void)
{
}

bool WeightedUnionFind::Connected(int p, int q)
{
    return Root(p) == Root(q);
}

void WeightedUnionFind::Union(int p, int q)
{
    int i = Root(p);
    int j = Root(q);

    // There is no need to add either i or j to the 'roots' set because they must be there already.
    if(size[i] < size[j])	// Attach the smaller tree to the larger one.
    {
        id[i] = j;			// Specify the larger cluster as a root, to keep the logarithmic tree height.
        size[j] += size[i];	// Increase the larger tree size by the merged amount.
        roots.erase(i);		// Erase i from the roots, as it is not root any more.
    }
    else
    {
        id[j] = i;
        size[i] += size[j];
        roots.erase(j);		// Erase j from the roots, as it is not root any more.
    }
}

// Return the root of the tree the vertex i belongs to.
int WeightedUnionFind::Root(int i)
{
    while(i != id[i])	// i is == to id[i] only for roots.
        i = id[i];		// Get the next id.

    return i;
}

// Set the initial root of the vertex.
void WeightedUnionFind::SetInitialRoot(int idp)
{
    if(size[idp] == 0)
    {
        size[idp] = 1;		// Specify the size and the first root.
        roots.insert(idp);	// Put the root into the set.
    }
}

void WeightedUnionFind::Reset(int N)
{
    id.resize(N);
    size.resize(N);

    for(int i = 0; i < N; ++i)
    {
        id[i] = i;		// Assign the id to the sequence number of the vertex.
        size[i] = 0;	// No elements in all the trees at the beginning.
    }

    roots.clear();		// Removes all the roots from the set.
}

//// Prints ids of all the vertices.
//void WeightedUnionFind::PrintId()
//{
//    std::cout << "Ids:" << id.size() << std::endl;
//
//    for(std::size_t i = 0; i < id.size(); ++i)
//        std::cout << i << " ";
//
//    std::cout << std::endl;
//
//    for(std::size_t i = 0; i < id.size(); ++i)
//        std::cout << id[i] << " ";
//
//    std::cout << std::endl;
//}
//
//// Prints sizes of clusters that correspond to tree vertices.
//void WeightedUnionFind::PrintSize()
//{
//    std::cout << "Tree sizes:" << std::endl;
//    for(std::size_t i = 0; i < size.size(); ++i)
//        std::cout << size[i] << " ";
//    std::cout << std::endl;
//}
//
//// Cluster statistics.
//void WeightedUnionFind::PrintClusters()
//{
//    std::cout << std::endl << "Statistics:" << std::endl;
//    std::cout << "Found clusters: " << roots.size() << std::endl;
//
//    // Print root of the trees and the corresponding tree sizes.
//    /*if(roots.size() > 0)
//    {
//    std::cout << "Root \t Size" << std::endl;
//    std::set<int>::iterator iter;
//    for(iter = roots.begin(); iter != roots.end(); iter++)
//    std::cout << *iter << " \t " << size[*iter] << std::endl;
//    }*/
//}
//
//// Output cluster histogram to a file.
//int WeightedUnionFind::PrintClusterSizeHistogram(const int bins, const std::string &fileOut)
//{
//    // Check whether there are clusters at all.
//    if(0 == roots.size())
//    {
//        std::cout << "No clusters found!" << std::endl;
//        return -1;
//    }
//
//    // Get sizes of the smallest/largest clusters.
//    GetMinMaxClusterSize(&minClusterSize, &maxClusterSize);
//
//    std::cout << "Min cluster: " << minClusterSize << std::endl;
//    std::cout << "Max cluster: " << maxClusterSize << std::endl;
//
//    // Width of the bin.
//    double delta = static_cast<double>(maxClusterSize - minClusterSize)/(bins - 1);
//    std::vector<int> sizeHistogram;
//    sizeHistogram.resize(bins);
//
//    if(delta > 1.e-07)	// Prevent from division by 0.
//    {
//        // Build the histogram.
//        std::set<int>::iterator iter;
//        for(iter = roots.begin(); iter != roots.end(); iter++)	// Loop over the clusters.
//        {
//            if(size[*iter] > 0)
//            {
//                int iChannel = static_cast<int>( (size[*iter] - minClusterSize)*1./delta + 0.5 );
//                ++sizeHistogram[iChannel];
//            }
//        }
//
//        // Print the histogram to the file.
//        std::ofstream histFile(fileOut);
//        histFile.setf(std::ios::fixed,std::ios::floatfield); //histFile.width(10);
//        double histSum = 0.;
//        for(int i = 0; i < bins; ++i)
//        {
//            histSum += static_cast<double>(sizeHistogram[i])/roots.size();			
//            histFile << i*delta + minClusterSize << "\t";
//            // Normalize only by the total number of islands (without delta).
//            histFile << static_cast<double>(sizeHistogram[i])/roots.size() << "\t" << histSum << std::endl;
//        }
//        histFile.close();
//    } // End if delta > 0.
//    else
//        std::cerr << "Warning: Don't have enough different cluster sizes for the histogram!" << std::endl;
//
//    return 0;
//}

// Define clusters with the minimum and maximum sizes.
void WeightedUnionFind::GetMinMaxClusterSize(int *min, int *max)
{
    std::set<int>::iterator iter;
    iter = roots.begin();
    int minCluster = size[*iter];
    int maxCluster = size[*iter];
    for(; iter != roots.end(); iter++)
    {
        if(minCluster > size[*iter])
            minCluster = size[*iter];
        if(maxCluster < size[*iter])
            maxCluster = size[*iter];
    }

    *min = minCluster;
    *max = maxCluster;
}

//// Compare min/max cluster sizes and their total number with the assumed values.
//bool WeightedUnionFind::CompareResults(int resultMinClusterSize, int resultMaxClusterSize, int numOfClusters)
//{
//    if(resultMinClusterSize != minClusterSize)
//        return false;
//
//    if(resultMaxClusterSize != maxClusterSize)
//        return false;
//
//    if(roots.size() != numOfClusters)
//        return false;
//
//    return true;
//}
