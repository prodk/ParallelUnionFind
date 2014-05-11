// Test client for the parallel union-find algorithm
// Nikolay Prodanov, 2014, Odessa

#include "TestParallelUnionFind.h"

int main(int argc, char **argv)
{
    TestParallelUnionFind testUF(argc, argv);

    testUF.runTests();

    return 0;
}
