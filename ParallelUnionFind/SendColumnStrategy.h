// SendColumnStrategy.h - the implementation of the class with this name.

#ifndef SEND_COLUMN_STRATEGY
#define SEND_COLUMN_STRATEGY

#include <vector>

struct Pixel;
struct DecompositionInfo;

class SendColumnStrategy
{
public:
    SendColumnStrategy(const DecompositionInfo & decompositionInfo, std::vector<Pixel> & pixels);
    virtual ~SendColumnStrategy(void);

protected:
    int getLeftNeighborProcessor() const;      // Periodic BCs are taken into account via DecompositionInfo.
    int getRightNeighborProcessor() const;     // Periodic BCs are taken into account via DecompositionInfo.
    bool isNeighborProcessorValid(const int rank) const;

protected:
    const DecompositionInfo & mDecompositionInfo;
    std::vector<Pixel> & mGlobalPixels;
};

#endif // SEND_COLUMN_STRATEGY
