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






/// Compute timing parameters
/**	@param vecEvtIndex : vecteur d'itérations de taille nb event
 *  @param vecPixelIdx : vecteur d'itérations de taille nb pixel
 *  @param signal : signal before cleaning
 *  @param cleanedsignal : signal cleaned (nbevent * nbpixel)
 *  @param fullmaxpos : position du maximum du signal dans les slices
 *  @param fullmax : maximum du signal
 *  @param nbpixel : nombre de pixel
 *  @param nbslice : nombre de slices
 *  @param weightY : somme du signal
 */
void compute_timing(const VecIndex& vecEvtIndex,
                    const VecIndex& vecPixelIdx,
                    float* signal, float* signalclean,
                    int* fullmaxpos,
                    float* fullmax,
					size_t nbPixel, size_t nbSlice,
					float* weightY, float* sumX, float* sumY, float* sumXY, float* sumX2, float* slope, float* intercept, float* longi,
					float* gx, float* gy, float* psi,
					float* tabPosX, float* tabPosY

) 
{	
	// itération sur le vecteur d'indices de taille nbEvent
	std::for_each(std::execution::par_unseq, std::begin(vecEvtIndex), std::end(vecEvtIndex),
        [=](int evtidx) {
			
			// Recherche du maximum du signal et sa position dans les slices
			// Itération sur le vecteur d'indice de taille nbPixel
            std::for_each(std::execution::unseq, std::begin(vecPixelIdx), std::end(vecPixelIdx),
            [&](int pixelidx) {
				fullmax[pixelidx] = 0.0f;
				fullmaxpos[pixelidx] = 0.0f;
				int fullmaxposvalue = fullmaxpos[pixelidx];
				float fullmaxvalue = fullmax[pixelidx];


                // itération sur toutes les slices d'un pixel d'un event
                for(size_t j(0lu); j < nbSlice; ++j){
                    float value = signal[evtidx * (nbSlice * nbPixel) + j*nbPixel + pixelidx];
                    float condition = value > fullmaxvalue;
                    fullmaxvalue = value * condition + (1.0f - condition) * fullmaxvalue;
                    fullmaxposvalue = j * condition + (1.0f - condition) * fullmaxposvalue;
                }
				
				fullmaxpos[pixelidx] /= fullmax[pixelidx];
				fullmaxpos[pixelidx] = fullmax[pixelidx]*(fullmax[pixelidx] >= 0.0f && fullmax[pixelidx] < nbSlice) + nbSlice*(fullmax[pixelidx] > nbSlice);
				
            });
			
			
			
			long firstPixelIdx = evtidx * nbPixel;

			float cosPsi(cosf(psi[evtidx])), sinPsi(sinf(psi[evtidx]));

			weightY[evtidx] = std::transform_reduce(std::execution::unseq, signalclean + firstPixelIdx, signalclean + firstPixelIdx + nbPixel, 0.0f, std::plus{},
                    [](auto val) {return val;});

			
			sumY[evtidx] = std::transform_reduce(std::execution::unseq, signalclean + firstPixelIdx, signalclean + firstPixelIdx + nbPixel, fullmaxpos, 0.0f, std::plus{},
                    [](auto val, auto y) {return val*y;});
			
			float g_x = gx[evtidx];
			float g_y = gy[evtidx];

			std::transform(std::execution::unseq, tabPosX, tabPosX + nbPixel, tabPosY, longi,
				[&](auto posx, auto posy) {
					return ((posx - g_x)*cosPsi + (posy - g_y)*sinPsi);
				});
			
			sumX[evtidx] = std::transform_reduce(std::execution::unseq, signalclean + firstPixelIdx, signalclean + firstPixelIdx + nbPixel, longi, 0.0f, std::plus{},
                    [](auto val, auto longi) {return val*longi;});
			
			sumXY[evtidx] = std::transform_reduce(std::execution::unseq, signalclean + firstPixelIdx, signalclean + firstPixelIdx + nbPixel, longi, 0.0f, std::plus{},
                    [](auto val, auto longi) {return val*longi;});

			sumX2[evtidx] = std::transform_reduce(std::execution::unseq, signalclean + firstPixelIdx, signalclean + firstPixelIdx + nbPixel, longi, 0.0f, std::plus{},
                    [](auto val, auto longi) {return val*longi*longi;});
			
			sumX[evtidx] /= (float)weightY[evtidx];
			sumY[evtidx] /= (float)weightY[evtidx];
			sumXY[evtidx] /= (float)weightY[evtidx];
			sumX2[evtidx] /= (float)weightY[evtidx];


			float sig2 = sumX2[evtidx] - sumX[evtidx]*sumX[evtidx];
			slope[evtidx] = (sumXY[evtidx]-sumX[evtidx]*sumY[evtidx])/sig2;
			intercept[evtidx] = sumY[evtidx] - sumX[evtidx]*slope[evtidx];
			
			
        });	
	
}





///Get the number of nanoseconds per elements of the Calibration
/** @param nbEvent : number of event of the tables
*/
void evaluateTiming(size_t nbEvent){
    // Let's define size of data:
    size_t eventSize(NB_PIXEL * NB_SLICE);
    size_t nbElement(nbEvent * eventSize);


    VecData Max(NB_PIXEL);

	VecData sumX(nbEvent), sumY(nbEvent), sumXY(nbEvent), sumX2(nbEvent), slope(nbEvent), intercept(nbEvent), weightY(nbEvent);

	VecData gx(nbEvent), gy(nbEvent), psi(nbEvent);

	VecData posX(NB_PIXEL), posY(NB_PIXEL), longi(NB_PIXEL);

	VecIndex PosMax(NB_PIXEL);

   
	// Le signal
	VecData vecSignal(nbElement), vecSignalClean(nbEvent*NB_PIXEL);
	std::fill(vecSignal.begin(), vecSignal.end(), 42.0f);
	

	// Le vecteur d'indices pour parcourir tout les pixels et les events
	VecIndex vecEvtIdx(nbEvent), vecPixelIdx(NB_PIXEL);
	std::iota(vecEvtIdx.begin(), vecEvtIdx.end(), 0);
	std::iota(vecPixelIdx.begin(), vecPixelIdx.end(), 0);


	
	// We have to create pointer to be able to catch them by copy without losing any time
	float *tabMax = Max.data(), *tabSignal = vecSignal.data(), *tabsumX = sumX.data(), *tabsumY = sumY.data(), *tabsumXY = sumXY.data(), *tabsumX2 = sumX2.data(), *tabslope = slope.data(), *tabintercept = intercept.data(), *tabweightY = weightY.data(), *tabsignalclean = vecSignalClean.data(), *tabgx = gx.data(), *tabgy = gy.data(), *tabpsi = psi.data(), *tabPosX = posX.data(), *tabPosY = posY.data(), *tablongi = longi.data(); 
	int *tabPosMax = PosMax.data();

	// Appel de microBenchMark pour tester la performance de la fonction de calibration
	size_t fullNbElement(nbEvent * NB_SLICE * NB_PIXEL);

    micro_benchmarkAutoNsPrint("evaluateCalibration 2d tranform", fullNbElement, compute_timing, vecEvtIdx, vecPixelIdx, tabSignal, tabsignalclean, tabPosMax, tabMax, NB_PIXEL, NB_SLICE, tabweightY, tabsumX, tabsumY, tabsumXY,tabsumX2, tabslope, tabintercept, tablongi, tabgx, tabgy, tabpsi, tabPosX, tabPosY);
}

int main(int argc, char** argv){
	return micro_benchmarkParseArg(argc, argv, evaluateTiming);
}