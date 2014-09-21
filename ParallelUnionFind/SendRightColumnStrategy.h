// SendRightColumnStrategy.h - the declaration of the SendRightColumnStrategy class.

#ifndef SEND_RIGHT_COLUMN_STRATEGY
#define SEND_RIGHT_COLUMN_STRATEGY

#include "SendColumnStrategy.h"

class SendRightColumnStrategy : public SendColumnStrategy
{
public:
    SendRightColumnStrategy(const DecompositionInfo & decompositionInfo,
                            const std::vector<std::ptrdiff_t> & localPixels,
                            const std::vector<Pixel> & globalPixels,
                            const std::tr1::shared_ptr<WeightedUnionFind> & localWuf,
                            std::map<std::ptrdiff_t, std::ptrdiff_t> & globalLabels);
    ~SendRightColumnStrategy(void);

private:
    // Overriden methods of the TM pattern.
    void copyPixelStripeToSend(SPixelStripe & stripeToSend) const;
    void sendStripeFromEvenReceiveOnOdd(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const;
    void sendStripeFromOddReceiveOnEven(SPixelStripe & stripeToSend, SPixelStripe & stripeToReceive) const;
    void saveReceivedStripe(const SPixelStripe & stripeToReceive, std::vector<Pixel> & globalPixels) const;

    // Helpers.
    void sendRightStripe(SPixelStripe & stripeToSend, const int msgId[], const int size) const;
    void receiveRightStripe(SPixelStripe & stripeToReceive, const int msgId[], const int size) const;
};

#endif // SEND_RIGHT_COLUMN_STRATEGY
