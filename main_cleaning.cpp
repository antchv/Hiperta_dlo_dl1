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
 *  @param tabInjonction - Table to convert hexagonal signal to/from matrix signal
 *  @param threshold_center - Threshold for the center pixel (current pixel)
 *  @param threshold_neighbour - Threshold for the neighbour of the current pixel
 *  @param nbPixel - Number of pixel in the signal
 *  @param nbCol - Number of columns in the matrix format
 *  @param nbRow - Number of row in the matrix format
 */
void compute_cleaning(float* tabSignalHex, float* tabSignalMat, int* offsets,
                        float* hightresh, float* keepsignal,
                        int* vecIndex,
                        float* tabInjonction,
                        float threshold_center, float threshold_neighbour,
                        long int minNumberPixelNeighbor,
						size_t nbPixel, size_t nbElement,
                        size_t nbCol, size_t nbRow)
{	// convert hexagonal signal to matrix
    std::for_each(std::execution::par_unseq, vecIndex, vecIndex + nbElement,
        [=](int idx) {
            // définition de l'indice de l'évenement et du pixel courant
			int evtidx = idx / nbPixel;
			int pixelidx = idx % nbPixel;
            int currentPixel = evtidx * nbPixel + pixelidx;

            // Convert hexagonal signal to matrix signal
            //size_t quadEvtSize(nbCol * nbRow);
            //size_t quadPixelIndex(idx / nbPixel * quadEvtSize + tabInjonction[idx % nbPixel]);
            //tabSignalMat[quadPixelIndex] = tabSignalHex[( idx / nbPixel) * nbPixel + (idx % nbPixel)];
            
            tabSignalMat[currentPixel] = tabSignalHex[currentPixel];
        });	
    
    // first pass, pixel > high threshold
    std::for_each(std::execution::par_unseq, vecIndex, vecIndex + nbElement,
        [=](int idx) {
            // définition de l'indice de l'évenement et du pixel courant
			int evtidx = idx / nbPixel;
			int pixelidx = idx % nbPixel;

            int currentPixel = evtidx * nbPixel + pixelidx;

            // First condition to check if the signal is above center threshold && si il a  nbvoisins > seuil bas suffisant
            float first_pass = (tabSignalMat[currentPixel] > threshold_center);

            float second_pass = ((tabSignalMat[currentPixel + offsets[0]]) + 
                                (tabSignalMat[currentPixel + offsets[1]]) +
                                (tabSignalMat[currentPixel + offsets[2]]) +
                                (tabSignalMat[currentPixel + offsets[3]]) +
                                (tabSignalMat[currentPixel + offsets[4]]) +
                                (tabSignalMat[currentPixel + offsets[5]]) ) >= minNumberPixelNeighbor;

            float cond = first_pass*second_pass;
            hightresh[currentPixel] = cond * tabSignalMat[currentPixel];
        });
    /*
    // second pass, enough neigbhours above high threshold
    std::for_each(std::execution::seq, vecIndex, vecIndex + nbElement,
    [=](int idx) {
        // définition de l'indice de l'évenement et du pixel courant
        int evtidx = idx / nbPixel;
        int pixelidx = idx % nbPixel;

        int currentPixel = evtidx * nbPixel + pixelidx;

        float temp_pass = ( (hightresh[currentPixel + offsets[0]]) + 
        (hightresh[currentPixel + offsets[1]]) +
        (hightresh[currentPixel + offsets[2]]) +
        (hightresh[currentPixel + offsets[3]]) +
        (hightresh[currentPixel + offsets[4]]) +
        (hightresh[currentPixel + offsets[5]]) ) >= minNumberPixelNeighbor;

        keepsignal[currentPixel] = temp_pass * tabSignalMat[currentPixel];
    });
    */

    // third pass, keep neighbours greater than thresholdNeighbours next to a selected pixel with a thresholdCenter
    std::for_each(std::execution::par_unseq, vecIndex, vecIndex + nbElement,
        [=](int idx) {
        // définition de l'indice de l'évenement et du pixel courant
        int evtidx = idx / nbPixel;
        int pixelidx = idx % nbPixel;

        int currentPixel = evtidx * nbPixel + pixelidx;

        // float cond = ((hightresh[currentPixel + offsets[0]]) + 
        //     (hightresh[currentPixel + offsets[1]]) +
        //     (hightresh[currentPixel + offsets[2]]) +
        //     (hightresh[currentPixel + offsets[3]]) +
        //     (hightresh[currentPixel + offsets[4]]) +
        //     (hightresh[currentPixel + offsets[5]]) ) >= 0.5f;

        float cond = 1.0f;

        float cond2 = tabSignalMat[currentPixel] > threshold_neighbour;
        
        keepsignal[currentPixel] = cond * cond2;


        // apply the condition
        tabSignalHex[currentPixel] = tabSignalHex[currentPixel] * keepsignal[currentPixel];	
		});
        	
}



/*
void compute_cleaningbis(float* tabSignalHex, float* tabSignalMat, int* offsets,
                        const VecIndex& vecIndex,
                        float* tabInjonction,
                        float threshold_center, float threshold_neighbour,
                        long int minNumberPixelNeighbor,
						size_t nbPixel,
                        size_t nbCol, size_t nbRow)
{	
    for(size_t idx(0lu) ; idx < vecIndex.size(); idx++) {
        // Convert hexagonal signal to matrix signal
        tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel)] = tabSignalHex[( idx / nbPixel) * nbPixel + (idx % nbPixel)];
    }
    

    for(size_t idx(0lu) ; idx < vecIndex.size(); idx++) {
        
        float first_pass = ((tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel)] > threshold_center)) * ( (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel) + offsets[0]] > threshold_neighbour) + 
        (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel) + offsets[1]] > threshold_neighbour) +
        (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel) + offsets[2]] > threshold_neighbour) +
        (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel) + offsets[3]] > threshold_neighbour) +
        (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel) + offsets[4]] > threshold_neighbour) +
        (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel) + offsets[5]] > threshold_neighbour) ) >= minNumberPixelNeighbor;


        float second_pass = ( (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel)] > threshold_neighbour) * (
            (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel) + offsets[0]] > threshold_center) + 
            (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel) + offsets[1]] > threshold_center) +
            (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel) + offsets[2]] > threshold_center) +
            (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel) + offsets[3]] > threshold_center) +
            (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel) + offsets[4]] > threshold_center) +
            (tabSignalMat[( idx / nbPixel) * nbPixel + (idx % nbPixel) + offsets[5]] > threshold_center) ) >= 0.5f
        );

        //float cond = (first_pass && second_pass) || third_pass; 
        float cond = first_pass || second_pass;  
                        
        // Convert matrix cleaned signal to hexagonal cleaned signal
        tabSignalHex[( idx / nbPixel) * nbPixel + (idx % nbPixel)] = tabSignalHex[( idx / nbPixel) * nbPixel + (idx % nbPixel)] * cond;	
    }
    
    
    
		
}
*/





///Get the number of nanoseconds per elements of the Calibration
/** @param nbEvent : number of event of the tables
*/
void evaluateCleaning(size_t nbEvent){
    // Let's define size of data:
    size_t nbElement(nbEvent * NB_PIXEL);
    int nbCol(55), nbRow(55);
    size_t minNumberPixelNeighbor(2);

    VecIndex offsets = {
        -nbCol - 1, -nbCol, -1, 1, nbCol, nbCol + 1
    };

	// définition des paramètres de la cleaning
    float threshold_center = 45.0, threshold_neighbour = 40.0f;
    VecData vecTabInjonction(NB_PIXEL);
    std::fill(vecTabInjonction.begin(), vecTabInjonction.end(), 1.0f);
	

	// Le signal à clean
	VecData vecSignalHex(nbElement), vecSignalMat(nbElement), vecHighTresh(nbElement), vecKeepSignal(nbElement);
	std::fill(vecSignalHex.begin(), vecSignalHex.end(), 42.0f);


	// Le vecteur d'indices pour parcourir tout les pixels et les events
	VecIndex vecIdx(NB_PIXEL*nbEvent);
	std::iota(vecIdx.begin(), vecIdx.end(), 0);



	
	// We have to create pointer to be able to catch them by copy without losing any time
	float *tabSignalHex = vecSignalHex.data(), *tabSignalMat = vecSignalMat.data(), *tabtabInjonction = vecTabInjonction.data(), *tabHighTresh = vecHighTresh.data(), *tabKeepSignal = vecKeepSignal.data();

	int *taboffsets = offsets.data(), *tabindex = vecIdx.data();
	// Appel de microBenchMark pour tester la performance de la fonction de calibration
	size_t fullNbElement(nbElement);

	micro_benchmarkAutoNsPrint("evaluateCalibration 2d tranform", fullNbElement, compute_cleaning, tabSignalHex, tabSignalMat, taboffsets, tabHighTresh, tabKeepSignal, tabindex, tabtabInjonction, threshold_center, threshold_neighbour, minNumberPixelNeighbor, NB_PIXEL, nbElement, nbCol, nbRow);
	
}

int main(int argc, char** argv){
	return micro_benchmarkParseArg(argc, argv, evaluateCleaning);
}