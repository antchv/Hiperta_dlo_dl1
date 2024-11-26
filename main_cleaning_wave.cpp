#include <iostream>
#include <vector>
#include <numeric>

#include <execution>
#include <algorithm>

#include "micro_benchmark.h"

#include "AlignedAllocator.h"
// #include <cblas.h>


///Defines a vector of index
typedef std::vector<int> VecIndex;

///Defines a vector of data
typedef std::vector<float, AlignedAllocator<float> > VecData;


// Fonction pour le produit matrice-vecteur en utilisant BLAS avec des VecData
// VecData produitMatriceVecteur(float* matrice, float* vecteur, int m) {
//     VecData resultat(m);
//     // Effectuer le produit matrice-vecteur en utilisant la fonction cblas_sgemv de BLAS (pour des données VecData)
//     cblas_sgemv(CblasRowMajor, CblasNoTrans, m, 1, 1.0f, matrice, 1, vecteur, 1, 1.0f, resultat.data(), 1);

//     return resultat;
// }




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

void compute_cleaning(float* tabSignalHex, float* tabSignalMat, int* offsetsNeighours,
                        const VecIndex& vecIndex, const VecIndex& vecEvtIndex,
                        float* tabInjonction,
                        float* threshold_wavelet, float* mu, float* somme, float* sigm,
						size_t nbPixel, size_t alpha,
                        size_t nbCol, size_t nbRow, float coef_center, float coef_neighbours,
                        float* tabone)
{	
    // peut etre faire un produit matrice vecteur pour faire la réduction ?
    //VecData matproj_sig = produitMatriceVecteur(tabSignalHex, tabone, vecEvtIndex.size());
    
    std::for_each(std::execution::par_unseq, std::begin(vecIndex), std::end(vecIndex),
        [=](int idx) {
            int evtidx = idx / nbPixel;
			int pixelidx = idx % nbPixel;
            int currentPixel = evtidx * nbPixel + pixelidx;

            // Faire un reduce?
            somme[evtidx] += tabSignalHex[currentPixel];
            //somme[evtidx] = matproj_sig[evtidx];
        });

    std::for_each(std::execution::par_unseq, std::begin(vecEvtIndex), std::end(vecEvtIndex),
        [=](int evtidx) {
            mu[evtidx] = somme[evtidx] / nbPixel;
        });

    std::for_each(std::execution::par_unseq, std::begin(vecIndex), std::end(vecIndex),
        [=](int idx) {
            int evtidx = idx / nbPixel;
			int pixelidx = idx % nbPixel;
            int currentPixel = evtidx * nbPixel + pixelidx;

            float c = tabSignalHex[currentPixel] - (tabSignalHex[currentPixel] / nbPixel);
            // faire un transform pour calculer un vecteur C de taille nbpixel ?
            sigm[evtidx] += (c*c);
            // faire sur vecevtidx et faire un reduce?
        });

    std::for_each(std::execution::par_unseq, std::begin(vecEvtIndex), std::end(vecEvtIndex),
        [=](int evtidx) {
            sigm[evtidx] = sqrt((sigm[evtidx]/((float)nbPixel) - mu[evtidx]*mu[evtidx]));
            threshold_wavelet[evtidx] = sqrt(alpha * log2(somme[evtidx])) * sigm[evtidx] + mu[evtidx];
        });
    
    // convert hexagonal signal to matrix
    // faire comme la version tailcut cpu
    std::for_each(std::execution::par_unseq, std::begin(vecIndex), std::end(vecIndex),
        [=](int idx) {
            int evtidx = idx / nbPixel;
			int pixelidx = idx % nbPixel;
            int currentPixel = evtidx * nbPixel + pixelidx;


            // Convert hexagonal signal to matrix signal
            //size_t quadEvtSize(nbCol * nbRow);
            //size_t quadPixelIndex(idx / nbPixel * quadEvtSize + tabInjonction[idx % nbPixel]);
            //tabSignalMat[quadPixelIndex] = tabSignalHex[( idx / nbPixel) * nbPixel + (idx % nbPixel)];
            
            tabSignalMat[currentPixel] = tabSignalHex[currentPixel];
        });	

    
    
    // compute condition to keep (or not) the signal in a pixel based on approximation coef
    // faire comme version tailcut cpu, sur vecevtindex puis boucle for sur les pixel
    std::for_each(std::execution::par_unseq, std::begin(vecIndex), std::end(vecIndex),
        [=](int idx) {
            int evtidx = idx / nbPixel;
			int pixelidx = idx % nbPixel;
            int currentPixel = evtidx * nbPixel + pixelidx;

            // Second condtion to check if there are enough neighbours which are above neighbours threshold
            // float cond = ( (tabSignalMat[currentPixel] * coef_center) +
            // (tabSignalMat[currentPixel + offsetsNeighours[0]] * coef_neighbours) + 
            // (tabSignalMat[currentPixel + offsetsNeighours[1]] * coef_neighbours) +
            // (tabSignalMat[currentPixel + offsetsNeighours[2]] * coef_neighbours) +
            // (tabSignalMat[currentPixel + offsetsNeighours[3]] * coef_neighbours) +
            // (tabSignalMat[currentPixel + offsetsNeighours[4]] * coef_neighbours) +
            // (tabSignalMat[currentPixel + offsetsNeighours[5]] * coef_neighbours) ) >= threshold_wavelet[idx];
            float cond = 1.0f;
            
            tabSignalHex[currentPixel] = cond*tabSignalHex[currentPixel];
		});	
    

}











/*
void compute_cleaningbis(float* tabSignalHex, float* tabSignalMat, int* offsetsNeighours,
                        const VecIndex& vecIndex, const VecIndex& vecEvtIndex,
                        float* tabInjonction,
                        float* threshold, float* mu, float* somme, float* sigm,
						size_t nbPixel, size_t alpha,
                        size_t nbCol, size_t nbRow, float coef_center, float coef_neighbours)
{	

    for(size_t idx(0lu) ; idx < vecIndex.size(); idx++) {
        // Convert hexagonal signal to matrix signal
        somme[( idx / nbPixel)] += tabSignalHex[( idx / nbPixel) * nbPixel + (idx % nbPixel)];
    }

    for(size_t idx(0lu) ; idx < vecEvtIndex.size(); idx++) {
        // Convert hexagonal signal to matrix signal
        mu[idx] = somme[idx] / nbPixel;
    }

    for(size_t idx(0lu) ; idx < vecIndex.size(); idx++) {
        // Convert hexagonal signal to matrix signal
        float c = tabSignalHex[( idx / nbPixel) * nbPixel + (idx % nbPixel)] - (tabSignalHex[( idx / nbPixel) * nbPixel + (idx % nbPixel)] / nbPixel);
        sigm[idx / nbPixel] += (c*c);
    }

    for(size_t idx(0lu) ; idx < vecEvtIndex.size(); idx++) {
        // Convert hexagonal signal to matrix signal
        sigm[idx] = sqrt((sigm[idx]/((float)nbPixel) - mu[idx]*mu[idx]));
        threshold[idx] = sqrt(alpha * log2(somme[idx])) * sigm[idx] + mu[idx];
    }

    

    for(size_t idx(0lu) ; idx < vecIndex.size(); idx++) {
        int currentPixel = ( idx / nbPixel) * nbPixel + (idx % nbPixel);

        float cond = ( (tabSignalMat[currentPixel] * coef_center) +
            (tabSignalMat[currentPixel + offsetsNeighours[0]] * coef_neighbours) + 
            (tabSignalMat[currentPixel + offsetsNeighours[1]] * coef_neighbours) +
            (tabSignalMat[currentPixel + offsetsNeighours[2]] * coef_neighbours) +
            (tabSignalMat[currentPixel + offsetsNeighours[3]] * coef_neighbours) +
            (tabSignalMat[currentPixel + offsetsNeighours[4]] * coef_neighbours) +
            (tabSignalMat[currentPixel + offsetsNeighours[5]] * coef_neighbours) ) >= threshold[idx];


                        
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
    int nbCol(55), nbRow(55), alpha(1);

    float coef_center = 20.0f / 55.0f;
    float coef_neighbours = 15.0f / 55.0f;

    VecIndex offsetsNeighours = {
        -nbCol - 1, -nbCol, -1, 1, nbCol, nbCol + 1
    };

	// définition des paramètres de la cleaning
    VecData threshold(nbEvent), somme(nbEvent), sigm(nbEvent), mu(nbEvent);
    VecData vecTabInjonction(NB_PIXEL);
    std::fill(vecTabInjonction.begin(), vecTabInjonction.end(), 1.0f);
	

	// Le signal à clean
	VecData vecSignalHex(nbElement), vecSignalMat(nbElement);
	std::fill(vecSignalHex.begin(), vecSignalHex.end(), 42.0f);


	// Le vecteur d'indices pour parcourir tout les pixels et les events
	VecIndex vecIdx(NB_PIXEL*nbEvent);
	std::iota(vecIdx.begin(), vecIdx.end(), 0);

    VecIndex vecEvtIndex(nbEvent);
	std::iota(vecEvtIndex.begin(), vecEvtIndex.end(), 0);


    VecData ones(NB_PIXEL);
	
	// We have to create pointer to be able to catch them by copy without losing any time
	float *tabSignalHex = vecSignalHex.data(), *tabSignalMat = vecSignalMat.data(), *tabtabInjonction = vecTabInjonction.data(), *tabsomme = somme.data(), *tabsigm = sigm.data(), *tabmu = mu.data(), *tabthreshold = threshold.data(), *tabone = ones.data();

	int *taboffsetsNeighours = offsetsNeighours.data();
	// Appel de microBenchMark pour tester la performance de la fonction de calibration
	size_t fullNbElement(nbElement);

	micro_benchmarkAutoNsPrint("evaluateCalibration 2d tranform", fullNbElement, compute_cleaning, tabSignalHex, tabSignalMat, taboffsetsNeighours, vecIdx, vecEvtIndex, tabtabInjonction, tabthreshold, tabmu, tabsomme, tabsigm, NB_PIXEL, alpha, nbCol, nbRow, coef_center, coef_neighbours, tabone);
	
}

int main(int argc, char** argv){
	return micro_benchmarkParseArg(argc, argv, evaluateCleaning);
}