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
    int mNumOfProc;
    int mMyRank;
};

#endif // TEST_PARALLEL_UNION_FIND
