#include "TestParallelUnionFind.h"

//---------------------------------------------------------------------------
TestParallelUnionFind::TestParallelUnionFind(int argc, char **argv) :
mNumOfProc(-1), mMyRank(-1), mPixels()
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
    test1();
#ifdef _DEBUG
    forceWindowToStay();
#endif
}

//---------------------------------------------------------------------------
void TestParallelUnionFind::test1()
{
    DecompositionInfo info = defineDecomposition();

    if(readPixels("picture_H08_4a_256_p_0_1.dat", info) >=0 )
    {
        // Specify pixels just read from the file.
        info.pixels = &mPixels[0];
        ParallelUnionFind puf("2DStripes", info);
        puf.analyze();
        puf.printClusterStatistics("fileName.txt");
    }
}

//---------------------------------------------------------------------------
DecompositionInfo TestParallelUnionFind::defineDecomposition()
{
    DecompositionInfo info;

    info.domainWidth = 256;
    info.domainHeight = 256;
    info.myRank = mMyRank;
    info.numOfProc = mNumOfProc;
    info.pixels = 0;

    return info;
}

//---------------------------------------------------------------------------
int TestParallelUnionFind::readPixels(const std::string& filePictureIn, const DecompositionInfo& info)
{
    std::cout << std::endl;
    std::cout << "________________" << std::endl;
    std::cout << "Reading    " << filePictureIn << std::endl;
    std::cout << "Creating " << info.domainWidth << " by " << info.domainHeight << " picture" << std::endl;

    mPixels.resize(info.domainWidth*info.domainHeight);

    std::ifstream fileIn(filePictureIn);
    // Read pixels from the file.
    std::string line;
    int lineCount = 0;
    while(std::getline(fileIn, line))
    {
        if(lineCount >= info.domainWidth)
        {
            std::cerr << std::endl;
            std::cerr << "Error: Too many lines in the picture file!" << std::endl;
            std::cerr << std::endl;
            return -1;
        }
        if(line.size() != info.domainHeight)
        {
            std::cerr << std::endl;
            std::cerr << "Error: Number of columns in the picture file is wrong!" << std::endl;
            std::cerr << std::endl;
            return -1;
        }

        // Copy pixels of the current line.
        char tmp[2] = "0";
        for(std::size_t i = 0; i < line.size(); ++i)
        {
            tmp[0] = line[i];
            mPixels[indexTo1D(lineCount, i, info)] = std::atoi(tmp);
        }

        ++lineCount;
    } // End while eof.

    fileIn.close();
    return 0;
}

//---------------------------------------------------------------------------
#ifdef _DEBUG
void TestParallelUnionFind::forceWindowToStay() const
{
    if(0 == mMyRank)
    {
        int dummy;
        std::cin >> dummy;
    }
}
#endif
