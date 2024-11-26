#include <iostream>
#include <vector>
#include <numeric>

#include <execution>
#include <algorithm>

#include "micro_benchmark.h"

#include <iostream>
#include <vector>
#include <numeric>

#include <execution>
#include <algorithm>



#include "AlignedAllocator.h"


///Defines a vector of index
typedef std::vector<int> VecIndex;

///Defines a vector of data
typedef std::vector<float, AlignedAllocator<float> > VecData;



///Compute the calibration
/** @param[out] CalibSignal - Pointer to the table of calibrated signals (input tabADC).
 *  @param signalhi - Pointer to the input high signal data.
 *  @param signallow - Pointer to the input low signal data. 
 *  @param mask - Pointer to the mask data.
 *  @param vecIndex - Vector of indices.
 *  @param pedHi - Pointer to the table of high value of pedestal.
 *  @param pedLow -Pointer to the table of low value of pedestal.
 *  @param gainHi - Pointer to the table of high value of gain.
 *  @param gainLow - Pointer to the table of low value of gain.
 *  @param threshold - Threshold value.
 *  @param nbSlice - Number of slices.
 *  @param nbPixel - Number of pixels.
 *  @param nbEvent - Number of events.
 */

void compute_calibration(float* calibSignal,
						float* signallow, float* signalhi,
						float* mask,
						const VecIndex& vecIndex,
						float* pedLow, float* pedHi,
						float* gainLow, float* gainHi,
						float threshold,
						size_t nbSlice, size_t nbPixel, size_t nbEvent)
{	

	// itération sur le vecteur d'indices de taille nbEvent * nbPixel
	std::for_each(std::execution::par_unseq, std::begin(vecIndex), std::end(vecIndex),
        [=](int idx) {
			// définition de l'indice de l'évenement et du pixel courant
			int evtidx = idx / nbPixel;
			int pixelidx = idx % nbPixel;

			// itération sur toutes les slices d'un pixel d'un event pour calculer le mask
			for(size_t j(0lu); j < nbSlice; ++j){
				mask[idx] *= (signalhi[evtidx * (nbSlice * nbPixel) + j*nbPixel + pixelidx] < threshold);
			}
			
			// itération sur toutes les slices d'un pixel d'un event pour appliquer le mask et faire la calibration
			for(size_t j(0lu); j < nbSlice; ++j){
				float gaincalib = (1.0 - mask[idx]) * gainLow[pixelidx] + mask[idx] * gainHi[pixelidx];
				float pedcalib = (1.0 - mask[idx]) * pedLow[pixelidx] + mask[idx] * pedHi[pixelidx];
				float sigcalib = (1.0f-mask[idx]) * signallow[evtidx * (nbSlice * nbPixel) + j*nbPixel + pixelidx] + mask[idx] * signalhi[evtidx * (nbSlice * nbPixel) + j*nbPixel + pixelidx];

				calibSignal[evtidx * (nbSlice * nbPixel) + j*nbPixel + pixelidx] = (sigcalib - pedcalib) * gaincalib;
				
			}
				
				
		});		
}



///Get the number of nanoseconds per elements of the Calibration
/** @param nbEvent : number of event of the tables
*/
void evaluateCalibration(size_t nbEvent){
    // Let's define size of data:
    size_t eventSize(NB_PIXEL * NB_SLICE);
    size_t nbElement(nbEvent * eventSize);

	// définition des paramètres de la calibration
	VecData vecGainHigh(NB_PIXEL), vecGainLow(NB_PIXEL), vecPedestalHigh(NB_PIXEL), vecPedestalLow(NB_PIXEL);
    std::fill(vecGainHigh.begin(), vecGainHigh.end(), 0.02f);
    std::fill(vecGainLow.begin(), vecGainLow.end(), 0.01f);
    std::fill(vecPedestalHigh.begin(), vecPedestalHigh.end(), 40.0f);
    std::fill(vecPedestalLow.begin(), vecPedestalLow.end(), 38.0f);
	
	
	float threshold = 45.0;
	

	// Le signal à calibrer ( en version high et low )
	VecData vecSignalhi(nbElement), vecSignallow(nbElement);
	std::fill(vecSignalhi.begin(), vecSignalhi.end(), 42.0f);
	std::fill(vecSignallow.begin(), vecSignallow.end(), 41.0f);

	// Le vecteur pour stocker les données calibrées
	VecData vecCalibSignal(nbElement);
	std::fill(vecCalibSignal.begin(), vecCalibSignal.end(), 0.0f);

	// Le vecteur d'indices pour parcourir tout les pixels et les events
	VecIndex vecIdx(NB_PIXEL*nbEvent);
	std::iota(vecIdx.begin(), vecIdx.end(), 0);

	// Le mask pour choisir le channel high ou low
	VecData vecMask(NB_PIXEL*nbEvent);
	std::fill(vecMask.begin(), vecMask.end(), 1.0f);

	
	// We have to create pointer to be able to catch them by copy without losing any time
	float *tabCalibSignal = vecCalibSignal.data(), *tabSignalhi = vecSignalhi.data(), *tabSignallow = vecSignallow.data(), *tabMask = vecMask.data();

	float *tabgainhi = vecGainHigh.data(), *tabgainlow = vecGainLow.data(), *tabpedhi = vecPedestalHigh.data(), *tabpedlow = vecPedestalLow.data();
	
	// Appel de microBenchMark pour tester la performance de la fonction de calibration
	size_t fullNbElement(nbElement);

	micro_benchmarkAutoNsPrint("evaluateCalibration 2d tranform", fullNbElement, compute_calibration, tabCalibSignal, tabSignallow, tabSignalhi, tabMask, vecIdx, tabpedlow, tabpedhi, tabgainlow, tabgainhi, threshold, NB_SLICE, NB_PIXEL, nbEvent);
	
}

int main(int argc, char** argv){
	return micro_benchmarkParseArg(argc, argv, evaluateCalibration);
}
