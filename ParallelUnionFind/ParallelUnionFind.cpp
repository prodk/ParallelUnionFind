// ParallelUnionFind.cpp - implementation of ParalleUnionFind class

#include "ParallelUnionFind.h"
#include "ParallelUnionFindFactory.h"

ParallelUnionFind::ParallelUnionFind(void):
mParallelUf(NULL)
{
}

ParallelUnionFind::ParallelUnionFind(const std::string& spatialConfiguration):
mParallelUf(ParallelUnionFindFactory::makeParallelUnionFind(spatialConfiguration))
{
}

ParallelUnionFind::~ParallelUnionFind(void)
{
}
