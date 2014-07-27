// SendLeftColumnStrategy.h - the declaration of the SendLeftColumnStrategy class.

#ifndef SEND_LEFT_COLUMN_STRATEGY
#define SEND_LEFT_COLUMN_STRATEGY

#include "SendColumnStrategy.h"

//---------------------------------------------------------------------------
class SendLeftColumnStrategy : public SendColumnStrategy
{
public:
    SendLeftColumnStrategy(const DecompositionInfo & decompositionInfo,
                           const std::vector<int> & localPixels,
                           const std::vector<Pixel> & globalPixels,
                           const std::tr1::shared_ptr<WeightedUnionFind> & localWuf,
                           std::map<int, int> & globalLabels);
    ~SendLeftColumnStrategy(void);

private:
    // Overriden methods of the TM pattern.
    void copyPixelStripeToSend(SPixelStripe & stripeToSend) const;
    void sendStripeFromEvenReceiveOnOdd(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const;
    void sendStripeFromOddReceiveOnEven(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const;
    void saveReceivedStripe(const SPixelStripe & stripeToReceive, std::vector<Pixel> & globalPixels) const;

    // Helpers.
    void sendLeftStripe(SPixelStripe & stripeToSend, const int msgId[], const int size) const;
    void receiveLeftStripe(SPixelStripe & stripeToReceive, const int msgId[], const int size) const;
};

#endif // SEND_LEFT_COLUMN_STRATEGY