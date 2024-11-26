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
/** @param[out] tabSignalHex - Pointer to the table of calibrated and integrated signal.
 *  @param tabSignalmat - Pointer to the table of signal in matrix format
 *  @param vecIndex - Vector of index to reach all of the pixels and events
 *  @param threshold_center - Threshold for the center pixel (current pixel)
 *  @param threshold_neighbour -Threshold for the neighbour of the current pixel
 *  @param nbPixel - Number of pixel in the signal
 *  @param nbEvent - Number of event in the signal
 *  @param nbCol - Number of columns in the matrix format
 *  @param nbRow - Number of row in the matrix format
 */
void compute_cleaning(float* tabSignalHex, float* tabSignalMat, float* first_pass, float* second_pass, float* third_pass,
                        const VecIndex& vecevtIndex, const VecIndex& vecpixelIndex,
                        float* tabInjonction,
                        float threshold_center, float threshold_neighbour,
                        long int minNumberPixelNeighbor,
						size_t nbPixel, size_t nbEvent,
                        size_t nbCol, size_t nbRow)
{	
	// itération sur le vecteur d'indices de taille nbEvent, paraelle seulement pour un nombre nbevent
	std::for_each(std::execution::par_unseq, std::begin(vecevtIndex), std::end(vecevtIndex),
        [=](int evtidx) {
        
        for (size_t i = 0; i < nbPixel; i++) {
            long current_pix = evtidx * nbPixel + i;
            //size_t quadEvtSize(55 * 55);
            //size_t quadPixelIndex(evtidx * quadEvtSize + tabInjonction[pixelidx]);

            // Convertir le signal hexagonal en signal matriciel
            tabSignalMat[current_pix] = tabSignalHex[current_pix];
            }
        
        // const size_t offsets[] = {
        //     -nbCol - 1lu, -nbCol, -1lu, 1lu, nbCol, nbCol + 1lu
        // }; 
        for (size_t i = 0; i < nbPixel; i++) {
            long current_pix = evtidx * nbPixel + i;
            float first_pass = (tabSignalMat[current_pix] > threshold_center);

            // float second_pass = ( (tabSignalMat[current_pix + offsets[0]] > threshold_neighbour) + 
            //     (tabSignalMat[current_pix + offsets[1]] > threshold_neighbour) +
            //     (tabSignalMat[current_pix + offsets[2]] > threshold_neighbour) +
            //     (tabSignalMat[current_pix + offsets[3]] > threshold_neighbour) +
            //     (tabSignalMat[current_pix + offsets[4]] > threshold_neighbour) +
            //     (tabSignalMat[current_pix + offsets[5]] > threshold_neighbour) ) >= minNumberPixelNeighbor;

            // float third_pass = ( (tabSignalMat[current_pix] > threshold_neighbour) * (
            //     (tabSignalMat[current_pix + offsets[0]] > threshold_center) + 
            //     (tabSignalMat[current_pix + offsets[1]] > threshold_center) +
            //     (tabSignalMat[current_pix + offsets[2]] > threshold_center) +
            //     (tabSignalMat[current_pix + offsets[3]] > threshold_center) +
            //     (tabSignalMat[current_pix + offsets[4]] > threshold_center) +
            //     (tabSignalMat[current_pix + offsets[5]] > threshold_center) ) >= 0.5f
            // );
            float second_pass = 1.0f;

            float third_pass = 1.0f;
                        
            float cond = (first_pass && second_pass) || third_pass;
                            
            // Convert matrix cleaned signal to hexagonal cleaned signal
            tabSignalHex[current_pix] = tabSignalHex[current_pix] * cond ;	

        }
        
        
        /*
        long firstPixelIdx = evtidx * nbPixel;
        std::transform(std::execution::unseq, tabSignalMat + firstPixelIdx, tabSignalMat + firstPixelIdx + nbPixel, first_pass,
            [=](float signal) {
                float res = (signal > threshold_center) ? 1.0f : 0.0f;
                return res;
            });
        // Précalcul des offsets pour les pixels voisins
        const size_t offsets[] = {
            -nbCol - 1lu, -nbCol, -1lu, 1lu, nbCol, nbCol + 1lu
        };    
    
        for (size_t i = 0; i < nbPixel; i++) {
            size_t quadPixelIndex = evtidx * nbPixel + i;
            float sum_second_pass = 0;

            for (int offset : offsets) {
                sum_second_pass += (tabSignalMat[quadPixelIndex + offset] > threshold_neighbour);
            }
            second_pass[i] = (sum_second_pass >= minNumberPixelNeighbor);
        }
        
        for (size_t i = 0; i < nbPixel; i++) {
            size_t quadPixelIndex = evtidx * nbPixel + i;
            float sum_third_pass = 0;
            for (int offset : offsets) {
                sum_third_pass += (tabSignalMat[quadPixelIndex + offset] > threshold_center);
            }
            third_pass[i] = (sum_third_pass > 0.5f);
        }

        for (size_t i = 0; i < nbPixel; i++) {
            float cond1 = first_pass[i] * second_pass[i];
            float cond2 = third_pass[i];

            tabSignalHex[evtidx * nbPixel + i] = tabSignalHex[evtidx * nbPixel + i] * cond1 * cond2;	
        }
        */
        
        
        //long firstPixelIdx = evtidx * nbPixel;
        /*
        for (size_t i = 0; i < nbPixel; i++) {
            //size_t quadEvtSize(55 * 55);
            //size_t quadPixelIndex(evtidx * quadEvtSize + tabInjonction[pixelidx]);

            // Convertir le signal hexagonal en signal matriciel
            tabSignalMat[evtidx * nbPixel + i] = tabSignalHex[evtidx * nbPixel + i];
            }
        
        
        // Précalcul des offsets pour les pixels voisins
        const size_t offsets[] = {
            -nbCol - 1lu, -nbCol, -1lu, 1lu, nbCol, nbCol + 1lu
        };   
        

        for (size_t i = 0; i < nbPixel; i++) {
            size_t quadPixelIndex = evtidx * nbPixel + i;

            first_pass[i] = (tabSignalMat[quadPixelIndex] > threshold_center);
        }
        
        for (size_t i = 0; i < nbPixel; i++) {
            size_t quadPixelIndex = evtidx * nbPixel + i;
            float sum_second_pass = 0;
            for (int offset : offsets) {
                sum_second_pass += (tabSignalMat[quadPixelIndex + offset] > threshold_neighbour);
            }
            second_pass[i] = (sum_second_pass >= minNumberPixelNeighbor);
        }

        for (size_t i = 0; i < nbPixel; i++) {

            tabSignalHex[evtidx * nbPixel + i] = tabSignalHex[evtidx * nbPixel + i] * first_pass[i] * second_pass[i];	
        }
        */
          
				
		});	
        	
}




///Get the number of nanoseconds per elements of the Calibration
/** @param nbEvent : number of event of the tables
*/
void evaluateCleaning(size_t nbEvent){
    // Let's define size of data:
    size_t nbElement(nbEvent * NB_PIXEL);
    size_t nbCol(55), nbRow(55);
    size_t minNumberPixelNeighbor(2);

	// définition des paramètres de la cleaning
    float threshold_center = 45.0, threshold_neighbour = 40.0f;
    VecData vecTabInjonction(NB_PIXEL), first_pass(NB_PIXEL), second_pass(NB_PIXEL), third_pass(NB_PIXEL);
    std::fill(vecTabInjonction.begin(), vecTabInjonction.end(), 1.0f);
	

	// Le signal à clean
	VecData vecSignalHex(nbElement), vecSignalMat(nbElement);
	std::fill(vecSignalHex.begin(), vecSignalHex.end(), 42.0f);


	// Le vecteur d'indices pour parcourir tout les pixels et les events
	VecIndex vecevtIdx(nbEvent), vecpixelIndex(NB_PIXEL);
	std::iota(vecevtIdx.begin(), vecevtIdx.end(), 0);
    std::iota(vecpixelIndex.begin(), vecpixelIndex.end(), 0);



	
	// We have to create pointer to be able to catch them by copy without losing any time
	float *tabSignalHex = vecSignalHex.data(), *tabSignalMat = vecSignalMat.data(), *tabtabInjonction = vecTabInjonction.data(), *tabfirstpass = first_pass.data(), *tabsecondpass = second_pass.data(), *tabthirdpass = third_pass.data();

	
	// Appel de microBenchMark pour tester la performance de la fonction de calibration
	size_t fullNbElement(nbElement);

	micro_benchmarkAutoNsPrint("evaluateCalibration 2d tranform", fullNbElement, compute_cleaning, tabSignalHex, tabSignalMat, tabfirstpass, tabsecondpass, tabthirdpass, vecevtIdx, vecpixelIndex, tabtabInjonction, threshold_center, threshold_neighbour, minNumberPixelNeighbor, NB_PIXEL, nbEvent, nbCol, nbRow);
	
}

int main(int argc, char** argv){
	return micro_benchmarkParseArg(argc, argv, evaluateCleaning);
}