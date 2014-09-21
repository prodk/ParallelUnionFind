// TestParallelUnionFind.h - TestParallelUnionFind class declaration.
// Provides kind of unit tests.

#ifndef TEST_PARALLEL_UNION_FIND
#define TEST_PARALLEL_UNION_FIND

//---------------------------------------------------------------------------
#include "ParallelUnionFindImpl.h"
#include <vector>
#include <iostream>
#include <ctype.h>
#include <mpi.h>

//---------------------------------------------------------------------------
class TestParallelUnionFind
{
public:
    TestParallelUnionFind(int argc, char **argv);
    ~TestParallelUnionFind(void);
    void runTests();

private:
    void readInputParameters();
    void analyze(DecompositionInfo& info, const std::string& pictureFilePath);
    void testTheSystem();
    void test8();
    void test256();
    void test1k();
    void test4k();
    void test8k();
    DecompositionInfo defineDecompositionInfo();
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
    void eraseComments(std::string& line, const std::string& comment) const;
    std::string trimString(const std::string& str, const std::string& whitespace);
    void saveInteger(int &number, const std::string& text);
    void printPartOfThePicture(const DecompositionInfo& info) const;
    void getFileNameFromPath(const std::string& filePath, std::string& fileName) const;

#ifdef _DEBUG
    void forceWindowToStay() const;
#endif

private:
    int mNumOfProc;
    int mMyRank;
    std::vector<int> mPixels;
    std::string mPictureFile;
    int mSystemSize;
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
        info.domainWidth = static_cast<std::ptrdiff_t>(size/mNumOfProc);
    }
    else
    {
        std::cout << "We need at least 1 processor to perform the task!" << std::endl;
        return;
    }
    info.domainHeight = static_cast<std::ptrdiff_t>(size);
    info.myRank = mMyRank;
    info.numOfProc = mNumOfProc;
    info.pixels = 0;
    info.periodicBoundaryX = false;
}

//---------------------------------------------------------------------------
inline std::string TestParallelUnionFind::trimString(const std::string& str, const std::string& whitespace = " \t")
{
    const std::size_t strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
    {
        return ""; // no content
    }

    const std::size_t strEnd = str.find_last_not_of(whitespace);
    const std::size_t strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

//---------------------------------------------------------------------------
inline void TestParallelUnionFind::saveInteger(int &number, const std::string& text)
{
    if (isdigit(text[0]))
    {
        number = atoi(text.c_str());
    }
    else
    {
        std::cout << "Cannot convert this string " << text << " to integer" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

//---------------------------------------------------------------------------
inline void TestParallelUnionFind::eraseComments(std::string& line, const std::string& comment) const
{
    std::size_t found = line.find_first_of(comment);
    if (std::string::npos != found)
    {
        line.erase(found);
    }
}

#endif // TEST_PARALLEL_UNION_FIND