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


///Compute the local peak integration
/** @param signal - Pointer to the table of calibrated signal.
 *  @param vecIndex - Vector of index
 *  @param nbSlice - number of slices per event
 *  @param nbPixel - number of pixels per slice
 *  @param nbSliceWindow - width of integration window
 *  @param nbSliceBeforePeak - number of slices before reach max value of signal 
 *  @param[out] integ - pointer to the table of integrated signal
 */
void compute_integration(float* signal,
						const VecIndex& vecIndex, const VecIndex& vecPixelIdx,
						size_t nbSlice, size_t nbPixel,
						size_t nbSliceWindow, size_t nbSliceBeforePeak,
						float* integ)
{	
	// itération sur le vecteur d'indices de taille nbEvent
	std::for_each(std::execution::par_unseq, std::begin(vecIndex), std::end(vecIndex),
        [=](int evtidx) {
            
            // on initialise la valeur max et son indice
			float fullmax = 0.0f;
			size_t fullmaxpos = 0;

			// Itération sur le vecteur d'indice de taille nbPixel
            std::for_each(std::execution::unseq, std::begin(vecPixelIdx), std::end(vecPixelIdx),
            [&](int pixelidx) {

                // itération sur toutes les slices d'un pixel d'un event
                for(size_t j(0lu); j < nbSlice; ++j){
                    float value = signal[evtidx * (nbSlice * nbPixel) + j*nbPixel + pixelidx];
                    float condition = value > fullmax;
                    fullmax = value * condition + (1.0f - condition) * fullmax;
                    fullmaxpos = j * condition + (1.0f - condition) * fullmaxpos;
                }

            });

			// Define integration window
			size_t begin = std::min(std::max(0lu, fullmaxpos - nbSliceBeforePeak), nbSlice - nbSliceWindow);
			size_t end = begin + nbSliceWindow;
			
			
			std::for_each(std::execution::unseq, std::begin(vecPixelIdx), std::end(vecPixelIdx),
            [=](int pixelidx) {
                // Compute local peak integration
                for(size_t j=begin; j < end; j++){
                    integ[evtidx*nbPixel + pixelidx] += signal[evtidx * (nbSlice * nbPixel) + j*nbPixel + pixelidx];
                }
            });

		});		
}





///Get the number of nanoseconds per elements of the Calibration
/** @param nbEvent : number of event of the tables
*/
void evaluateIntegration(size_t nbEvent){
    // Let's define size of data:
    size_t eventSize(NB_PIXEL * NB_SLICE);
    size_t nbElement(nbEvent * eventSize);

	// Les paramètre d'intégration
	size_t nbSliceWindow = 2;
	size_t nbSliceBeforePeak = 1;


	// Le signal à intégrer
	VecData vecSignal(nbElement);
	std::fill(vecSignal.begin(), vecSignal.end(), 42.0f);

	// Le vecteur d'indice pour parcourir tout les pixels et tout les events
	VecIndex vecIdx(nbEvent), vecPixelIdx(NB_PIXEL);
	std::iota(vecIdx.begin(), vecIdx.end(), 0);
    std::iota(vecPixelIdx.begin(), vecPixelIdx.end(), 0);

	// Le vecteur qui va stocker l'intégration
	VecData integ(NB_PIXEL*nbEvent);
	

	float *tabSignal = vecSignal.data(), *tabinteg = integ.data();

	// Appel de microBenchMark pour tester la performance de la fonction d'intégration
	size_t fullNbElement(nbElement);
	micro_benchmarkAutoNsPrint("evaluateCalibration 2d tranform", fullNbElement, compute_integration, tabSignal, vecIdx, vecPixelIdx, NB_SLICE, NB_PIXEL, nbSliceWindow, nbSliceBeforePeak, tabinteg);
	
}

int main(int argc, char** argv){
	return micro_benchmarkParseArg(argc, argv, evaluateIntegration);
}