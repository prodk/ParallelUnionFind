// TestParallelUnionFind.cpp - implementation of the TestParallelUnionFind class.

#include "TestParallelUnionFind.h"
#include "ParallelUnionFind.h"
#include <fstream>
#include <sstream>
#include <mpi.h>

//---------------------------------------------------------------------------
TestParallelUnionFind::TestParallelUnionFind(int argc, char **argv)
    : mNumOfProc(-1)
    , mMyRank(-1)
    , mPixels()
    , mNumOfBins(6)
{
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &mNumOfProc);
    MPI_Comm_rank (MPI_COMM_WORLD, &mMyRank);
}

//---------------------------------------------------------------------------
TestParallelUnionFind::~TestParallelUnionFind(void)
{
    // Cleanup MPI.
    MPI_Finalize();
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::runTests()
{
    test8();
    //test256();
    //test1k();
    //test4k();
    //test8k();
#ifdef _DEBUG
    forceWindowToStay();
#endif
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::analyze(DecompositionInfo& info, const std::string& fileIn)
{
    if (readPixelsInParallel(fileIn, info) >= 0)
    {
        // Specify pixels just read from the file.
        info.pixels = &mPixels[0];
        ParallelUnionFind puf("2DStripes", info);

        // Analyze contact (pixel value 1 by default).
        puf.analyze();

#ifdef _DEBUG
        puf.printPerProcessorClusterStatistics("");
        puf.printPerProcessorClusterSizes("");
#endif

        puf.printClusterStatistics("");
        puf.printClusterSizeHistogram(mNumOfBins, "cont" + fileIn);

        //int bins = 1000;
        //puf.printClusterSizeHistogram(bins, "cont"+fileIn);

        // Analyze non-contact.
        /*std::cout << std::endl;
        const int pixelvalue = 0;
        puf.setPixelValue(pixelvalue);
        puf.analyze();
        puf.printPerProcessorClusterStatistics("");
        puf.printPerProcessorClusterSizes("");*/
        //puf.printPerProcessorClusterSizeHistogram(bins, "ncont"+fileIn);
    }
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::test8()
{
    DecompositionInfo info = defineDecomposition8();
    analyze(info, "8by8.dat");
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::test256()
{
    DecompositionInfo info = defineDecomposition256();
    analyze(info, "picture_H08_4a_256_p_0_1.dat");
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::test1k()
{
    DecompositionInfo info = defineDecomposition1k();
    analyze(info, "picture_H08_4a_1k_p_0_1.dat");
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::test4k()
{
    DecompositionInfo info = defineDecomposition4k();
    analyze(info, "picture_H03_2a_4k_p_0_035.dat");
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::test8k()
{
    DecompositionInfo info = defineDecomposition8k();
    analyze(info, "picture_H08_4a_8k_p_0_1.dat");
}

//---------------------------------------------------------------------------
DecompositionInfo TestParallelUnionFind::defineDecomposition8()
{
    DecompositionInfo info;

    const int size = 8;
    fillInDecompositionInfo(info, size);

    return info;
}

//---------------------------------------------------------------------------
DecompositionInfo TestParallelUnionFind::defineDecomposition256()
{
    DecompositionInfo info;

    const int size = 256;
    fillInDecompositionInfo(info, size);

    return info;
}

//---------------------------------------------------------------------------
DecompositionInfo TestParallelUnionFind::defineDecomposition1k()
{
    DecompositionInfo info;

    const int size = 1024;
    fillInDecompositionInfo(info, size);

    return info;
}

//---------------------------------------------------------------------------
DecompositionInfo TestParallelUnionFind::defineDecomposition4k()
{
    DecompositionInfo info;

    const int size = 4096;
    fillInDecompositionInfo(info, size);

    return info;
}

//---------------------------------------------------------------------------
DecompositionInfo TestParallelUnionFind::defineDecomposition8k()
{
    DecompositionInfo info;

    const int size = 8192;
    fillInDecompositionInfo(info, size);

    return info;
}

//---------------------------------------------------------------------------
int TestParallelUnionFind::readPixels(std::vector<int>& pixels,
                                      const std::string& filePictureIn,
                                      const DecompositionInfo& info)
{
    std::cout << "________________" << std::endl;
    std::cout << "Reading    " << filePictureIn << std::endl;

    std::ifstream fileIn(filePictureIn);
    if (fileIn.good())
    {
        // Read pixels from the file.
        std::string line;
        std::size_t lineCount = 0;
        while (std::getline(fileIn, line))
        {
            if (lineCount >= mNumOfProc*info.domainWidth)
            {
                std::cerr << std::endl;
                std::cerr << "Error: Too many lines in the picture file!" << std::endl;
                std::cerr << std::endl;
                return -1;
            }
            if (line.size() != info.domainHeight)
            {
                std::cerr << std::endl;
                std::cerr << "Error: Number of columns in the picture file is wrong!" << std::endl;
                std::cerr << std::endl;
                return -1;
            }

            // Copy pixels of the current line.
            char tmp[2] = "0";
            for (std::size_t i = 0; i < line.size(); ++i)
            {
                tmp[0] = line[i];
                pixels[indexTo1D(i, lineCount, info)] = std::atoi(tmp);
            }

            ++lineCount;
        } // End while eof.

        fileIn.close();
    }
    return 0;
}

//---------------------------------------------------------------------------
int TestParallelUnionFind::readPixelsInParallel(const std::string& filePictureIn, const DecompositionInfo& info)
{
    // Set the size of the chunk of data residing on the current processor.
    mPixels.resize(info.domainWidth*info.domainHeight);

    // For testing purposes only: read the whole file and then choose the data that the current proc needs.
    std::vector<int> allPixels;
    allPixels.resize(info.domainHeight*info.domainWidth*mNumOfProc);

    readPixels(allPixels, filePictureIn, info);

    copyRelevantData(allPixels, info);

    // Print the part of the picture residing on the current processor.
    printPartOfThePicture(info);

    return 0;
}

void TestParallelUnionFind::copyRelevantData(const std::vector<int>& allPixels, const DecompositionInfo& info)
{
    // Choose the data that corresponds to the current processor.
    const int nx = info.domainWidth;
    const int ny = info.domainHeight;
    for (int ix = 0; ix < nx; ++ix)            // Loop through the pixels.
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            int globalX = ix + nx*mMyRank;
            int globalIndex = indexTo1D(globalX, iy, info);
            int localIndex = indexTo1D(ix, iy, info);
            mPixels[localIndex] = allPixels[globalIndex];
        }
    }
}

//---------------------------------------------------------------------------
#ifdef _DEBUG
void TestParallelUnionFind::forceWindowToStay() const
{
    if (0 == mMyRank)
    {
        int dummy;
        std::cin >> dummy;
    }
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::printPartOfThePicture(const DecompositionInfo& info) const
{
    std::stringstream fileName;
    fileName << "proc_" << mMyRank << "_picture.dat";
    std::ofstream outFile(fileName.str());

    if (outFile.good())
    {
        const int nx = info.domainWidth;
        const int ny = info.domainHeight;

        for (int iy = 0; iy < ny; ++iy)            // Loop through the pixels, rows fastest.
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                outFile << mPixels[indexTo1D(ix, iy, info)];
            }
            outFile << std::endl;
        }
        outFile.close();
    }
}
#endif