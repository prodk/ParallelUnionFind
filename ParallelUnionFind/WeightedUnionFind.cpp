// WeightedUnionFind.cpp - WeightedUnionFind class implementation.
#include "WeightedUnionFind.h"

//---------------------------------------------------------------------------
WeightedUnionFind::WeightedUnionFind(const std::size_t N)
{
    mId.resize(N);
    mSize.resize(N);

    for(std::size_t i = 0; i < N; ++i)
    {
        mId[i] = i;      // Assign the id to the sequence number of the vertex.
        mSize[i] = 0;    // No elements in all the trees at the beginning.
    }

    mMinClusterSize = 0;
    mMaxClusterSize = 0;
}

//---------------------------------------------------------------------------
WeightedUnionFind::~WeightedUnionFind(void)
{
}

//---------------------------------------------------------------------------
bool WeightedUnionFind::connected(int p, int q)
{
    return root(p) == root(q);
}

//---------------------------------------------------------------------------
void WeightedUnionFind::makeUnion(int p, int q)
{
    int i = root(p);
    int j = root(q);

    // There is no need to add either i or j to the 'roots' set because they must be there already.
    if(mSize[i] < mSize[j])   // Attach the smaller tree to the larger one.
    {
        mId[i] = j;          // Specify the larger cluster as a root, to keep the logarithmic tree height.
        mSize[j] += mSize[i]; // Increase the larger tree size by the merged amount.
        mRoots.erase(i);     // Erase i from the roots, as it is not root any more.
    }
    else
    {
        mId[j] = i;
        mSize[i] += mSize[j];
        mRoots.erase(j);     // Erase j from the roots, as it is not root any more.
    }
}

//---------------------------------------------------------------------------
// Return the root of the tree the vertex i belongs to.
int WeightedUnionFind::root(int i)
{
    while(i != mId[i])      // i is == to id[i] only for roots.
    {
        i = mId[i];          // Get the next id.
    }

    return i;
}

//---------------------------------------------------------------------------
// Set the initial root of the vertex.
void WeightedUnionFind::setInitialRoot(int idp)
{
    if(0 == mSize[idp])
    {
        mSize[idp] = 1;      // Specify the size and the first root.
        mRoots.insert(idp);   // Put the root into the set.
    }
}

//---------------------------------------------------------------------------
void WeightedUnionFind::reset(int N)
{
    mId.resize(N);
    mSize.resize(N);

    for(int i = 0; i < N; ++i)
    {
        mId[i] = i;      // Assign the id to the sequence number of the vertex.
        mSize[i] = 0;    // No elements in all the trees at the beginning.
    }

    mRoots.clear();      // Removes all the roots from the set.
}

//---------------------------------------------------------------------------
// Prints ids of all the vertices.
void WeightedUnionFind::printId()
{
    std::cout << "Ids:" << mId.size() << std::endl;

    for(std::size_t i = 0; i < mId.size(); ++i)
        std::cout << i << " ";

    std::cout << std::endl;

    for(std::size_t i = 0; i < mId.size(); ++i)
    {
        std::cout << mId[i] << " ";
    }

    std::cout << std::endl;
}

//---------------------------------------------------------------------------
// Prints sizes of clusters that correspond to tree vertices.
void WeightedUnionFind::printSize()
{
    std::cout << "Tree sizes:" << std::endl;
    for(std::size_t i = 0; i < mSize.size(); ++i)
        std::cout << mSize[i] << " ";
    std::cout << std::endl;
}

//---------------------------------------------------------------------------
// Cluster statistics.
void WeightedUnionFind::printClusters()
{
    std::cout << std::endl << "Statistics:" << std::endl;
    std::cout << "Found clusters: " << mRoots.size() << std::endl;

    // Print root of the trees and the corresponding tree sizes.
    if(mRoots.size() > 0)
    {
        std::cout << "Root \t Size" << std::endl;
        std::set<int>::iterator iter;
        for(iter = mRoots.begin(); iter != mRoots.end(); iter++)
        {
            std::cout << *iter << " \t " << mSize[*iter] << std::endl;
        }
    }
}

//---------------------------------------------------------------------------
// Output cluster histogram to a file.
int WeightedUnionFind::printClusterSizeHistogram(const int bins, const std::string &fileOut)
{
    // Check whether there are clusters at all.
    if(0 == mRoots.size())
    {
        std::cout << "No clusters found!" << std::endl;
        return -1;
    }

    // Get sizes of the smallest/largest clusters.
    getMinMaxClusterSize(&mMinClusterSize, &mMaxClusterSize);

    std::cout << "Min cluster: " << mMinClusterSize << std::endl;
    std::cout << "Max cluster: " << mMaxClusterSize << std::endl;

    // Width of the bin.
    double delta = static_cast<double>(mMaxClusterSize - mMinClusterSize)/(bins - 1);
    std::vector<int> sizeHistogram;
    sizeHistogram.resize(bins);

    if(delta > 1.e-07)	// Prevent from division by 0.
    {
        // Build the histogram.
        std::set<int>::iterator iter;
        for(iter = mRoots.begin(); iter != mRoots.end(); iter++)	// Loop over the clusters.
        {
            if(mSize[*iter] > 0)
            {
                int iChannel = static_cast<int>( (mSize[*iter] - mMinClusterSize)*1./delta + 0.5 );
                ++sizeHistogram[iChannel];
            }
        }

        // Print the histogram to the file.
        std::ofstream histFile(fileOut);
        histFile.setf(std::ios::fixed,std::ios::floatfield); //histFile.width(10);
        double histSum = 0.;
        for(int i = 0; i < bins; ++i)
        {
            histSum += static_cast<double>(sizeHistogram[i])/mRoots.size();			
            histFile << i*delta + mMinClusterSize << "\t";
            // Normalize only by the total number of islands (without delta).
            histFile << static_cast<double>(sizeHistogram[i])/mRoots.size() << "\t" << histSum << std::endl;
        }
        histFile.close();
    } // End if delta > 0.
    else
        std::cerr << "Warning: Don't have enough different cluster sizes for the histogram!" << std::endl;

    return 0;
}

//---------------------------------------------------------------------------
// Define clusters with the minimum and maximum sizes.
void WeightedUnionFind::getMinMaxClusterSize(int *min, int *max)
{
    std::set<int>::iterator iter;
    iter = mRoots.begin();
    int minCluster = mSize[*iter];
    int maxCluster = mSize[*iter];
    for(; iter != mRoots.end(); iter++)
    {
        if(minCluster > mSize[*iter])
        {
            minCluster = mSize[*iter];
        }
        if(maxCluster < mSize[*iter])
        {
            maxCluster = mSize[*iter];
        }
    }

    *min = minCluster;
    *max = maxCluster;
}

//---------------------------------------------------------------------------
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
