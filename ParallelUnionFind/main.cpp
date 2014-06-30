// Test client for the parallel union-find algorithm
// (c) Mykola Prodanov, 2014, Odessa, Ukraine.

//---------------------------------------------------------------------------
#include "TestParallelUnionFind.h"

//---------------------------------------------------------------------------
int main(int argc, char **argv)
{
    TestParallelUnionFind testUF(argc, argv);

    testUF.runTests();

    return 0;
}
