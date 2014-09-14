// TestParallelUnionFind.h - TestParallelUnionFind class declaration.
// Provides kind of unit tests.

#ifndef TEST_PARALLEL_UNION_FIND
#define TEST_PARALLEL_UNION_FIND

//---------------------------------------------------------------------------
#include "ParallelUnionFindImpl.h"
#include <vector>
#include <iostream>

//---------------------------------------------------------------------------
class TestParallelUnionFind
{
public:
    TestParallelUnionFind(int argc, char **argv);
    ~TestParallelUnionFind(void);
    void runTests();

private:
    void analyze(DecompositionInfo& info, const std::string& fileIn);
    void test8();
    void test256();
    void test1k();
    void test4k();
    void test8k();
    DecompositionInfo defineDecomposition8();
    DecompositionInfo defineDecomposition256();
    DecompositionInfo defineDecomposition1k();
    DecompositionInfo defineDecomposition4k();
    DecompositionInfo defineDecomposition8k();
    void fillInDecompositionInfo(DecompositionInfo& info, const int size);

    int readPixels(std::vector<int>& pixels, const std::string& filePictureIn, const DecompositionInfo& info);
    int readPixelsInParallel(const std::string& filePictureIn, const DecompositionInfo& info);
    void copyRelevantData(const std::vector<int>& allPixels, const DecompositionInfo& info);
    int indexTo1D(int ix, int iy, const DecompositionInfo& info) const;

#ifdef _DEBUG
    void forceWindowToStay() const;
    void printPartOfThePicture(const DecompositionInfo& info) const;
#endif

private:
    int mNumOfProc;
    int mMyRank;
    std::vector<int> mPixels;
    int mNumOfBins;
};

//---------------------------------------------------------------------------
inline int TestParallelUnionFind::indexTo1D(int ix, int iy, const DecompositionInfo& info) const
{
    return ix*info.domainHeight + iy;
}

//---------------------------------------------------------------------------
inline void TestParallelUnionFind::fillInDecompositionInfo(DecompositionInfo& info, const int size)
{
    if (mNumOfProc > 0)
    {
        info.domainWidth = size/mNumOfProc;
    }
    else
    {
        std::cout << "We need at least 1 processor to perform the task!" << std::endl;
        return;
    }
    info.domainHeight = size;
    info.myRank = mMyRank;
    info.numOfProc = mNumOfProc;
    info.pixels = 0;
    info.periodicBoundaryX = false;
}

#endif // TEST_PARALLEL_UNION_FIND