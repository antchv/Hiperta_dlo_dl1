#include <iostream>
#include <vector>
#include <numeric>

#include <execution>
#include <algorithm>

#include "micro_benchmark.h"

#include "AlignedAllocator.h"




///Defines a vector of index
typedef std::vector<int> VecIndex;

///Defines a vector of data
typedef std::vector<float, AlignedAllocator<float> > VecData;


///Compute the calibration
/** @param[out] tabCalibSignal : table of signal (input tabADC)
 *  @param tabPedHigh : pointer to the table of pedestal for high threshold
 *  @param tabPedLow : pointer to the table of pedestal for low threshold
 *  @param tabGainHigh : pointer to the table of gain for high threshold
 *  @param tabGainLow : pointer to the table of gain for low threshold
 *  @param tabGain : pointer to the table of gain which will be used for calibration
 *  @param tabPed : pointer to the table of pedestal which will be used for calibration
 *  @param threshold : threshold for selecting high or low values
 *  @param tabEventIdx : vector of index of events
 *  @param tabSliceIdx : vector of index of slices
 *  @param nbPixel : number of pixels
 *  @param nbSlice : number of slices
 *  @param nbEvent : number of events
*/
void compute_calibration(float* tabCalibSignal, float* tabPedHigh, float* tabPedLow, 
                         float* tabGainHigh, float* tabGainLow, float* tabGain, float* tabPed, float* tabtemppedlow, float* tabtempgainlow, float threshold, int* tabEventIdx, int* tabSliceIdx, size_t nbPixel, size_t nbSlice, size_t nbEvent) {

    // We iterate on the events
    std::for_each(std::execution::par_unseq, tabEventIdx, tabEventIdx + nbEvent,
        [&](int evtIdx) {   

            // We create a new global mask for each event
            std::vector<float> globalMask(nbPixel, 0.0f);


            // We iterate on the slices
            std::for_each(std::execution::unseq, tabSliceIdx, tabSliceIdx + nbSlice,
                [&](int sliceIdx) {

                    long firstPixelIdx = evtIdx * (nbPixel * nbSlice) + sliceIdx * nbPixel;


                    // Update the mask
                    std::transform(std::execution::unseq, tabCalibSignal + firstPixelIdx, tabCalibSignal + firstPixelIdx + nbPixel, globalMask.begin(), globalMask.begin(),
                        [=](float signal, float globalVal) {
                            return (signal >= threshold) ? globalVal : 1.0f; 
                        });
                });
            
            
            

            
            // On applique le mask a tabPedLow, ressort un tabPedLow avec 0 si le signal ne dépasse jamais le seuil et la valeur initial sinon
            std::transform(std::execution::unseq, tabPedLow, tabPedLow + nbPixel, globalMask.begin(), tabtemppedlow,
                [=](float pedlow, float mask) -> float {
                    return mask * pedlow;
                }
            );
            // tabPed va contenir les valeurs de tabPedLow sauf pour les éléments égaux à zéro, où il va contenir les valeurs correspondantes de tabPedHigh
            std::transform(std::execution::unseq, tabtemppedlow, tabtemppedlow + nbPixel, tabPedHigh, tabPed,
                [=](float pedlow, float pedhigh) -> float {
                    return (pedlow == 0.0f) ? pedhigh : pedlow ;
                }
            );

            
            // On applique le mask a tabGainLow, ressort un tabGainLow avec 0 si le signal ne dépasse jamais le seuil et la valeur initial sinon
            std::transform(std::execution::unseq, tabGainLow, tabGainLow + nbPixel, globalMask.begin(), tabtempgainlow,
                [=](float gainlow, float mask) -> float {
                    return mask * gainlow;
                }
            );
            // tabGain va contenir les valeurs de tabGainLow sauf pour les éléments égaux à zéro, où il va contenir les valeurs correspondantes de tabGainHigh
            std::transform(std::execution::unseq, tabtempgainlow, tabtempgainlow + nbPixel, tabGainHigh, tabGain,
                [=](float gainlow, float gainhigh) -> float {
                    return (gainlow == 0.0f) ? gainhigh : gainlow ;
                }
            );
            
            

            
            // We iterate on the slices to apply calibration
            std::for_each(std::execution::unseq, tabSliceIdx, tabSliceIdx + nbSlice,
                [&](int sliceIdx) {
                    long firstPixelIdx = evtIdx * (nbPixel * nbSlice) + sliceIdx * nbPixel;

                    std::transform(std::execution::unseq, tabCalibSignal + firstPixelIdx, tabCalibSignal + firstPixelIdx + nbPixel, tabPed, tabCalibSignal + firstPixelIdx,
                        [=](float signal, float ped){
                            return signal - ped;
                        }
                    );
                    // And then gain
                    std::transform(std::execution::unseq, tabCalibSignal + firstPixelIdx, tabCalibSignal + firstPixelIdx + nbPixel, tabGain, tabCalibSignal + firstPixelIdx,
                        [=](float signal, float gain){
                            return signal * gain;
                        }
                    );
                }); 
        }
    );

}










///Get the number of nanoseconds per elements of the Calibration
/** @param nbEvent : number of event of the tables
*/
void evaluateCalibration(size_t nbEvent){
    // Let's define size of data:
    size_t eventSize(NB_PIXEL * NB_SLICE);
    size_t nbElement(nbEvent * eventSize);
    VecData vecGainHigh(NB_PIXEL), vecGainLow(NB_PIXEL), vecPedestalHigh(NB_PIXEL), vecPedestalLow(NB_PIXEL), vecGain(NB_PIXEL), vecPedestal(NB_PIXEL), vectemppedlow(NB_PIXEL), vectempgainlow(NB_PIXEL);



    std::fill(vecGainHigh.begin(), vecGainHigh.end(), 0.02f);
    std::fill(vecGain.begin(), vecGain.end(), 0.0f);
    std::fill(vecGainLow.begin(), vecGainLow.end(), 0.01f);
    std::fill(vecPedestalHigh.begin(), vecPedestalHigh.end(), 40.0f);
    std::fill(vecPedestal.begin(), vecPedestal.end(), 0.0f);
    std::fill(vecPedestalLow.begin(), vecPedestalLow.end(), 38.0f);

    
    VecData vecADCSignal(nbElement);
    std::fill(vecADCSignal.begin(), vecADCSignal.end(), 42.0f);
    

    VecIndex vecEventIdx(nbEvent);   
    std::iota(vecEventIdx.begin(), vecEventIdx.end(), 0);

    VecIndex vecSliceIdx(NB_SLICE);   
    std::iota(vecSliceIdx.begin(), vecSliceIdx.end(), 0);


    // We have to create pointer to be able to catch them by copy without losing any time
    float *tabADC = vecADCSignal.data(), *tabGainHigh = vecGainHigh.data(), *tabGainLow = vecGainLow.data(), *tabPedHigh = vecPedestalHigh.data(), *tabPedLow = vecPedestalLow.data(), *tabPed = vecPedestal.data(), *tabGain = vecGain.data(), *tabtemppedlow = vectemppedlow.data(), *tabtempgainlow = vectempgainlow.data();


    int *tabSliceIdx = vecSliceIdx.data(), *tabEventIdx = vecEventIdx.data();

    size_t fullNbElement(nbElement);
    float threshold = 45.0f; 

    
    micro_benchmarkAutoNsPrint("evaluateCalibration 2d transform", fullNbElement, compute_calibration, tabADC, tabPedHigh, tabPedLow, tabGainHigh, tabGainLow, tabGain, tabPed, tabtemppedlow, tabtempgainlow, threshold, tabEventIdx, tabSliceIdx, NB_PIXEL, NB_SLICE, nbEvent);    
}



int main(int argc, char** argv){
    return micro_benchmarkParseArg(argc, argv, evaluateCalibration);
}
