// WeightedUnionFind.cpp - WeightedUnionFind class implementation.
#include "WeightedUnionFind.h"
#include <iostream>
#include <fstream>

//---------------------------------------------------------------------------
// Avoid magic numbers.
const int WeightedUnionFind::defaultInt = -1;
//---------------------------------------------------------------------------
WeightedUnionFind::WeightedUnionFind(const std::ptrdiff_t N)
    : mRoots()
    , mConsecutiveRoots()
{
    mId.resize(N);
    mSize.resize(N);

    for (std::ptrdiff_t i = 0; i < N; ++i)
    {
        mId[i] = i;      // Assign the id to the sequence number of the vertex.
        mSize[i] = 0;    // No elements in all the trees at the beginning.
    }

    mMinClusterSize = defaultInt;
    mMaxClusterSize = defaultInt;
}

//---------------------------------------------------------------------------
WeightedUnionFind::~WeightedUnionFind(void)
{
}

//---------------------------------------------------------------------------
bool WeightedUnionFind::connected(std::ptrdiff_t p, std::ptrdiff_t q)
{
    return root(p) == root(q);
}

//---------------------------------------------------------------------------
void WeightedUnionFind::makeUnion(std::ptrdiff_t p, std::ptrdiff_t q)
{
    std::ptrdiff_t i = root(p);
    std::ptrdiff_t j = root(q);

    // There is no need to add either i or j to the 'roots' set because they must be there already.
    if (mSize[i] < mSize[j])   // Attach the smaller tree to the larger one.
    {
        mId[i] = j;           // Specify the larger cluster as a root, to keep the logarithmic tree height.
        mSize[j] += mSize[i]; // Increase the larger tree size by the merged amount.
        mRoots.erase(i);      // Erase i from the roots, as it is not root any more.
    }
    else
    {
        mId[j] = i;
        mSize[i] += mSize[j];
        mRoots.erase(j);      // Erase j from the roots, as it is not root any more.
    }
}

//---------------------------------------------------------------------------
// Return the root of the tree the vertex i belongs to.
std::ptrdiff_t WeightedUnionFind::root(std::ptrdiff_t i) const
{
    while (i != mId[i])      // i is == to id[i] only for roots.
    {
        i = mId[i];          // Get the next id.
    }

    return i;
}

//---------------------------------------------------------------------------
// Set the initial root of the vertex.
void WeightedUnionFind::setInitialRoot(std::ptrdiff_t idp)
{
    if (mSize[idp] <= 0)
    {
        mSize[idp] = 1;        // Specify the size and the first root.
        mRoots.insert(idp);    // Put the root into the set.
    }
}

//---------------------------------------------------------------------------
void WeightedUnionFind::setInitialRoot(std::ptrdiff_t idp, std::ptrdiff_t clusterSize)
{
    if (mSize[idp] <= 0)
    {
        mSize[idp] = clusterSize; // Specify the size and the first root.
        mRoots.insert(idp);       // Put the root into the set.
    }
}

//---------------------------------------------------------------------------
void WeightedUnionFind::reset(std::ptrdiff_t N)
{
    mId.resize(N);
    mSize.resize(N);

    for (std::ptrdiff_t i = 0; i < N; ++i)
    {
        mId[i] = i;            // Assign the id to the sequence number of the vertex.
        mSize[i] = 0;          // No elements in all the trees at the beginning.
    }

    mRoots.clear();            // Remove all the roots from the set.
    mConsecutiveRoots.clear(); // Clear consecutive roots.

    // Important when WUF is reused: reset cluster sizes.
    mMinClusterSize = defaultInt;
    mMaxClusterSize = defaultInt;
}

//---------------------------------------------------------------------------
const std::map<std::ptrdiff_t, std::ptrdiff_t>& WeightedUnionFind::getConsecutiveRootIds()
{
    mConsecutiveRoots.clear();

    const std::ptrdiff_t numOfClusters = mRoots.size();
    std::ptrdiff_t count = 0;
    std::set<std::ptrdiff_t>::iterator iter;
    for (iter = mRoots.begin(); iter != mRoots.end(); ++iter)
    {
        mConsecutiveRoots[*iter] = count; // Root id is a key, consecutive id is a value.
        ++count;
    }

    return mConsecutiveRoots;
}

//---------------------------------------------------------------------------
// Prints sizes of clusters that correspond to tree vertices.
void WeightedUnionFind::printClusterSizes(std::ostream& out)
{
    // Print root of the trees and the corresponding tree sizes.
    if (mRoots.size() > 0)
    {
        out << "Loc  LocConsec  Size" << std::endl;
        std::set<std::ptrdiff_t>::iterator iter;
        for (iter = mRoots.begin(); iter != mRoots.end(); ++iter)
        {
            out << *iter << " \t " << mConsecutiveRoots[*iter] << " \t " << mSize[*iter] << std::endl;
        }
    }
}

//---------------------------------------------------------------------------
// Cluster statistics.
void WeightedUnionFind::printClusterStatistics(std::ostream& out)
{
    out << "Found clusters: " << mRoots.size() << std::endl;

    if ((mMinClusterSize < 0) || (mMaxClusterSize < 0))
    {
        getMinMaxClusterSize(&mMinClusterSize, &mMaxClusterSize);
    }

    out << "Min cluster: " << mMinClusterSize << std::endl;
    out << "Max cluster: " << mMaxClusterSize << std::endl;
}

//---------------------------------------------------------------------------
// Output cluster histogram to a file.
int WeightedUnionFind::printClusterSizeHistogram(const int bins, const std::string &fileOut)
{
    // Check whether there are clusters at all.
    if (0 == mRoots.size())
    {
        std::cout << "No clusters found!" << std::endl;
        return defaultInt;
    }

    // Get sizes of the smallest/largest clusters.
    if ((mMinClusterSize < 0) || (mMaxClusterSize < 0))
    {
        getMinMaxClusterSize(&mMinClusterSize, &mMaxClusterSize);
    }

    buildAndPrintSizeHistogram(bins, fileOut);

    return 0;
}

//---------------------------------------------------------------------------
// Define clusters with the minimum and maximum sizes.
void WeightedUnionFind::getMinMaxClusterSize(std::ptrdiff_t *min, std::ptrdiff_t *max)
{
    std::set<std::ptrdiff_t>::iterator iter;
    iter = mRoots.begin();
    std::ptrdiff_t minCluster = mSize[*iter];
    std::ptrdiff_t maxCluster = mSize[*iter];
    for (; iter != mRoots.end(); ++iter)
    {
        if (minCluster > mSize[*iter])
        {
            minCluster = mSize[*iter];
        }
        if (maxCluster < mSize[*iter])
        {
            maxCluster = mSize[*iter];
        }
    }

    *min = minCluster;
    *max = maxCluster;
}

//---------------------------------------------------------------------------
void WeightedUnionFind::buildAndPrintSizeHistogram(const int bins, const std::string &fileOut) const
{
    // Width of the bin.
    const double delta = static_cast<double>(mMaxClusterSize - mMinClusterSize)/(bins - 1);
    std::vector<int> sizeHistogram;
    sizeHistogram.resize(bins);

    if (delta > std::numeric_limits<double>::epsilon())  // Prevent from division by 0.
    {
        // Build the histogram.
        std::set<std::ptrdiff_t>::iterator iter;
        for (iter = mRoots.begin(); iter != mRoots.end(); ++iter)  // Loop over the clusters.
        {
            if (mSize[*iter] > 0)
            {
                int iChannel = static_cast<int>( (mSize[*iter] - mMinClusterSize)*1./delta + 0.5 );
                ++sizeHistogram[iChannel];
            }
        }

        // Print the histogram to the file.
        std::ofstream histFile(fileOut);
        const std::ptrdiff_t numOfRoots = static_cast<std::ptrdiff_t>(mRoots.size());

        if (histFile.good() && (numOfRoots > 0))
        {
            histFile.setf(std::ios::fixed,std::ios::floatfield); //histFile.width(10);
            double histSum = 0.;
            for (int i = 0; i < bins; ++i)
            {
                histSum += static_cast<double>(sizeHistogram[i])/numOfRoots;
                histFile << i*delta + mMinClusterSize << "\t";
                // Normalize only by the total number of islands (without delta).
                histFile << static_cast<double>(sizeHistogram[i])/numOfRoots << "\t" << histSum << std::endl;
            }

            histFile.close();
        }
        else
        {
            std::cout << "Error: Bad histogram file or 0 number of clusters." << std::endl;
        }
    } // End if delta > 0.
    else
    {
        std::cout << "Warning: Don't have enough different cluster sizes for the histogram!" << std::endl;
    }
}
