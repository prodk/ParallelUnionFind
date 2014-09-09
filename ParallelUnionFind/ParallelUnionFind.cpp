// ParallelUnionFind.cpp - implementation of ParalleUnionFind class

#include "ParallelUnionFind.h"
#include "ParallelUnionFindFactory.h"
#include <iostream>

//---------------------------------------------------------------------------
ParallelUnionFind::ParallelUnionFind(const std::string& spatialConfiguration,
                                     const DecompositionInfo& info) :
mParallelUf(ParallelUnionFindFactory::makeParallelUnionFind(spatialConfiguration, info))
{
}

//---------------------------------------------------------------------------
ParallelUnionFind::~ParallelUnionFind(void)
{
}

//---------------------------------------------------------------------------
void ParallelUnionFind::analyze(void)
{
#ifdef _DEBUG
    logMsg(std::cout, "PUF::analyze()");
#endif

    mParallelUf->analyze();
}

//---------------------------------------------------------------------------
void ParallelUnionFind::printClusterStatistics(const std::string& fileName) const
{
#ifdef _DEBUG
    logMsg(std::cout, "\nPUF::printClusterStatistics()");
#endif

    mParallelUf->printClusterStatistics(fileName);
}

//---------------------------------------------------------------------------
void ParallelUnionFind::printClusterSizes(const std::string& fileName) const
{
#ifdef _DEBUG
    logMsg(std::cout, "PUF::printClusterSizes()");
#endif

    mParallelUf->printClusterSizes(fileName);
}

//---------------------------------------------------------------------------
void ParallelUnionFind::printClusterSizeHistogram(const int bins, const std::string& fileName) const
{
#ifdef _DEBUG
    logMsg(std::cout, "PUF::printClusterSizeHistogram()");
#endif

    mParallelUf->printClusterSizeHistogram(bins, fileName);
}

//---------------------------------------------------------------------------
void ParallelUnionFind::setPixelValue(const int value)
{
#ifdef _DEBUG
    logMsg(std::cout, "PUF::setPixelValue()");
#endif

    mParallelUf->setPixelValue(value);
}

//---------------------------------------------------------------------------
void ParallelUnionFind::printPerProcessorClusterStatistics(const std::string& fileName) const
{
#ifdef _DEBUG
    logMsg(std::cout, "\nPUF::printPerProcessorClusterStatistics()");
#endif

    mParallelUf->printPerProcessorClusterStatistics(fileName);
}

//---------------------------------------------------------------------------
void ParallelUnionFind::printPerProcessorClusterSizes(const std::string& fileName) const
{
#ifdef _DEBUG
    logMsg(std::cout, "\nPUF::printPerProcessorClusterSizes()");
#endif

    mParallelUf->printPerProcessorClusterSizes(fileName);
}

//---------------------------------------------------------------------------
void ParallelUnionFind::printPerProcessorClusterSizeHistogram(const int bins, const std::string& fileName) const
{
#ifdef _DEBUG
    logMsg(std::cout, "\nPUF::printPerProcessorClusterSizeHistogram()");
#endif

    mParallelUf->printPerProcessorClusterSizeHistogram(bins, fileName);
}

//---------------------------------------------------------------------------
#ifdef _DEBUG
void ParallelUnionFind::logMsg(std::ostream& out, const std::string& msg) const
{
    out << msg << std::endl;
}
#endif
