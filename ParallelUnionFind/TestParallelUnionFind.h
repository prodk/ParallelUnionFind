// TestParallelUnionFind.h - TestParallelUnionFind class declaration.
// Provides kind of unit test.

#ifndef TEST_PARALLEL_UNION_FIND
#define TEST_PARALLEL_UNION_FIND

#include "ParallelUnionFind.h"

class TestParallelUnionFind
{
public:
    TestParallelUnionFind(int argc, char **argv);
    ~TestParallelUnionFind(void);
    void runTests();

private:
    void analyze(DecompositionInfo& info, const std::string fileIn);
    void test256();
    void test1k();
    void test4k();
    void test8k();
    DecompositionInfo defineDecomposition256();
    DecompositionInfo defineDecomposition1k();
    DecompositionInfo defineDecomposition4k();
    DecompositionInfo defineDecomposition8k();

    int readPixels(const std::string& filePictureIn, const DecompositionInfo& info);
    inline int indexTo1D(int ix, int iy, const DecompositionInfo& info) const;

#ifdef _DEBUG
    void forceWindowToStay() const;
#endif

private:
    int mNumOfProc;
    int mMyRank;
    std::vector<int> mPixels;
};

//---------------------------------------------------------------------------
inline int TestParallelUnionFind::indexTo1D(int ix, int iy, const DecompositionInfo& info) const
{
    return ix*info.indexFactor*info.domainHeight + iy;
}

#endif // TEST_PARALLEL_UNION_FIND
