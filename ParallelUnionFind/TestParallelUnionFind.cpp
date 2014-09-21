// TestParallelUnionFind.cpp - implementation of the TestParallelUnionFind class.
// Note: this test client is assumed to be run on a small number of processors (order of 1 - 10)
// and small systems (linear size ~ 10^1 - 10 ^3)
// because it does not use any parallel I/O which will be a bottleneck if the mentioned conditions are violated.

#include "TestParallelUnionFind.h"
#include "ParallelUnionFind.h"
#include <fstream>
#include <sstream>

//---------------------------------------------------------------------------
TestParallelUnionFind::TestParallelUnionFind(int argc, char **argv)
    : mNumOfProc(-1)
    , mMyRank(-1)
    , mPixels()
    , mPictureFile("8by8.dat")
    , mSystemSize(8)
    , mNumOfBins(6)
{
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &mNumOfProc);
    MPI_Comm_rank (MPI_COMM_WORLD, &mMyRank);

    readInputParameters();
}

//---------------------------------------------------------------------------
TestParallelUnionFind::~TestParallelUnionFind(void)
{
    // Cleanup MPI.
    MPI_Finalize();
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::readInputParameters()
{
    // Every processor reads the input parameters.
    // This will be a bottleneck when we have a big number of processors.
    // But we assume that this test client will not run on many processors.

    // We assume a file 'input.txt' that contains the picture file path, system size (square), number of bins, BCs.
    const std::string inputFile = "input.txt";

    // Read and parse the file.
    std::ifstream fileIn(inputFile);
    if (fileIn.good())
    {
        const std::string comment = "#";   // Symbol for comments.
        const size_t numOfParams = 3;      // Number of parameters in the input.txt
        std::string line;
        std::size_t paramCount = 0;

        while ((paramCount < numOfParams) && std::getline(fileIn, line))
        {
            // Look for comments.
            eraseComments(line, comment);

            // Trim the white spaces and save the parameter.
            std::string paramLine = trimString(line);
            if ("" != paramLine)
            {
                switch (paramCount)
                {
                case 0:       // mPictureFile
                    mPictureFile = paramLine;
                    break;
                case 1:       // mSystemSize
                    saveInteger(mSystemSize, paramLine);
                    break;
                case 2:
                    break;    // mNumOfBins
                    saveInteger(mNumOfBins, paramLine);
                }
                ++paramCount;
            }
        } // End while eof.

        fileIn.close();
    }
    else
    {
        std::cerr << "Cannot read the input file." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::runTests()
{
    testTheSystem();
    // We can call here the built-in methods such as test8() etc.
#ifdef _DEBUG
    forceWindowToStay();
#endif
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::analyze(DecompositionInfo& info, const std::string& fileIn)
{
    if (readPixelsInParallel(fileIn, info) >= 0)
    {
        if (0 == info.myRank)
        {
            std::cout << "__________________________" << std::endl;
        }
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

        // Analyze non-contact.
        info.pixelValue = 0;
        if (0 == info.myRank)
        {
            std::cout << "__________________________" << std::endl;
        }
        ParallelUnionFind nonContPuf("2DStripes", info);
        nonContPuf.analyze();
        nonContPuf.printClusterStatistics("");
        nonContPuf.printClusterSizeHistogram(mNumOfBins, "noncont" + fileIn);
    }
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::testTheSystem()
{
    DecompositionInfo info = defineDecompositionInfo();
    analyze(info, mPictureFile);
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
DecompositionInfo TestParallelUnionFind::defineDecompositionInfo()
{
    DecompositionInfo info;

    const int size = mSystemSize;
    fillInDecompositionInfo(info, size);

    return info;
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
    if (0 == mMyRank)
    {
        std::cout << "Reading    " << filePictureIn << std::endl;
    }

    std::ifstream fileIn(filePictureIn);
    if (fileIn.good())
    {
        // Read pixels from the file.
        std::string line;
        std::ptrdiff_t lineCount = 0;
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
    else
    {
        std::cerr << "Error: Cannot read the picture file." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    return 0;
}

//---------------------------------------------------------------------------
int TestParallelUnionFind::readPixelsInParallel(const std::string& filePictureIn, const DecompositionInfo& info)
{
    // Set the size of the chunk of data residing on the current processor.
    mPixels.resize(info.domainWidth*info.domainHeight);

    // For testing purposes only: read the whole file and then choose the data that the current proc needs.
    // This will be a bottleneck when we have a big number of processors.
    // But we assume that this test client will not run on many processors.
    std::vector<int> allPixels;
    allPixels.resize(info.domainHeight * info.domainWidth * mNumOfProc);

    readPixels(allPixels, filePictureIn, info);

    copyRelevantData(allPixels, info);

#ifdef _DEBUG
    // Print the part of the picture residing on the current processor.
    printPartOfThePicture(info);
#endif

    return 0;
}

//---------------------------------------------------------------------------
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
#endif