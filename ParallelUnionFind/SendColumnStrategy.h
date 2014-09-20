// SendColumnStrategy.h - the declaration of the SendColumnStrategy class.

#ifndef SEND_COLUMN_STRATEGY
#define SEND_COLUMN_STRATEGY

#include <vector>
#include <map>
#include <cstddef> // For ptrdiff_t.
                   // ptrdiff_t is necessary for huge systems. It works even when the integer value exceeds INT_MAX = 2^32.

struct Pixel;
struct DecompositionInfo;
class WeightedUnionFind;

//---------------------------------------------------------------------------
// A helper data structure containing attributes of the pixel stripe sent between processors.
struct SPixelStripe
{
    std::vector<std::ptrdiff_t> pixelValue;
    std::vector<std::ptrdiff_t> globalClusterId;
    std::vector<std::ptrdiff_t> sizeOfCluster;

    explicit SPixelStripe(std::ptrdiff_t size)
        : pixelValue(size)
        , globalClusterId(size)
        , sizeOfCluster(size)
    {}
};

//---------------------------------------------------------------------------
class SendColumnStrategy
{
public:
    SendColumnStrategy(const DecompositionInfo & decompositionInfo,
                       const std::vector<int> & localPixels,
                       const std::vector<Pixel> & globalPixels,
                       const std::tr1::shared_ptr<WeightedUnionFind> & localWuf,
                       std::map<std::ptrdiff_t, std::ptrdiff_t> & globalLabels);
    virtual ~SendColumnStrategy(void);

    void sendReceivePixelStripes(std::vector<Pixel> & globalPixels);

protected:
    // Methods of the template method pattern to override.
    virtual void copyPixelStripeToSend(SPixelStripe & stripeToSend) const = 0;
    virtual void sendStripeFromEvenReceiveOnOdd(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const = 0;
    virtual void sendStripeFromOddReceiveOnEven(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const = 0;
    virtual void saveReceivedStripe(const SPixelStripe & stripeToReceive, std::vector<Pixel> & globalPixels) const = 0;

protected:
    int getLeftNeighborProcessor() const;      // Periodic BCs are taken into account via DecompositionInfo.
    int getRightNeighborProcessor() const;     // Periodic BCs are taken into account via DecompositionInfo.
    bool isNeighborProcessorValid(const int rank) const;

protected:
    const DecompositionInfo & mDecompositionInfo;
    const std::vector<int> & mLocalPixels;
    const std::vector<Pixel> & mGlobalPixels;
    const std::tr1::shared_ptr<WeightedUnionFind> & mLocalWuf;

    std::map<std::ptrdiff_t, std::ptrdiff_t> & mGlobalLabels;// This reference is not modified, so should have been const.
                                                             // But if it is const we cannot use [] operator. So I left it non-const.

    enum {INVALID_VALUE = -1, BOSS = 0, NUM_OF_SENDS = 3, MSG_1 =123, MSG_2 = 456, MSG_3 = 789};
};

#endif // SEND_COLUMN_STRATEGY
