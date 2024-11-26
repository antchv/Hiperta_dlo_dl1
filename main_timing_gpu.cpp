#include <iostream>
#include <vector>
#include <numeric>

#include <execution>
#include <algorithm>

#include "micro_benchmark.h"

#include "AlignedAllocator.h"


#include <cblas.h> 


///Defines a vector of index
typedef std::vector<int> VecIndex;

///Defines a vector of data
typedef std::vector<float, AlignedAllocator<float> > VecData;




// Fonction pour le produit matrice-vecteur en utilisant BLAS avec des VecData
VecData produitMatriceVecteur(float* matrice, float* vecteur, int m) {
    VecData resultat(m);
    // Effectuer le produit matrice-vecteur en utilisant la fonction cblas_sgemv de BLAS (pour des données VecData)
    cblas_sgemv(CblasRowMajor, CblasNoTrans, m, 1, 1.0f, matrice, 1, vecteur, 1, 1.0f, resultat.data(), 1);

    return resultat;
}







/// Compute timing parameters
/**	@param vecIndex : vecteur d'indices de taille nbEvent * nbPixel
 *  @param vecEvtIdx : vecteur d'indices de taille nbEvent
 *  @param[out] timingSlope : Slope of the time fit (nbEvent)
 * 	@param[out] timingIntercept : Intercept of the time fit (nbEvent)
 * 	@param[out] matTimingSumX : Matrix of the longitudinal projection (nbEvent * nbPixel)
 * 	@param[out] matTimingSumY : Matrix of the tabPosMax projection (nbEvent * nbPixel)
 * 	@param[out] matTimingSumXY : Matrix of the cross projection (nbEvent * nbPixel)
 * 	@param[out] matTimingSumX2 : Matrix of the 2nd moment projection (nbEvent * nbPixel)
 * 	@param[out] tabTimingWeight : Table of the reduced weight (nbEvent)
 * 	@param[out] tabTimingSumX : Table of the longitudinal projection (nbEvent)
 * 	@param[out] tabTimingSumY : Table of the tabPosMax projection (nbEvent)
 * 	@param[out] tabTimingSumXY : Table of the cross projection (nbEvent)
 * 	@param[out] tabTimingSumX2 : Table of the 2nd moment projection (nbEvent)
 * 	@param tabKeepSignalHex : table of the cleaned signal of all event to be analysed (nbEvent, nbPixel)
 * 	@param hillasX : Barycenter on x axis (m) (size nbPixel)
 * 	@param hillasY : Barycenter on y axis (m) (size nbPixel)
 * 	@param tabCosPsi : Table of cos psi (nbEvent)
 * 	@param tabSinPsi : Table of sin psi (nbEvent)
 * 	@param tabPosPixelX : Positions of the pixels on the x axis (nbPixel)
 * 	@param tabPosPixelY : Positions of the pixels on the y axis (nbPixel)
 * 	@param tabPosMax : table of the index of the maximum value (nbEvent * nbPixel)
 * 	@param tabone : table of ones use to reduce events (nbPixel)
 * 	@param nbEvent : number of events to be analysed
 * 	@param nbPixel : number of pixels per event
 */

void compute_timing(const VecIndex& vecIndex, const VecIndex& vecEvtIdx,
                    float* tabone,
                    float* tabKeepSignalHex,
                    float* tabPosMax,
                    float* hillasX,
                    float* hillasY,
                    float* tabCosPsi,
                    float* tabSinPsi,
                    float* tabPosPixelX,
                    float* tabPosPixelY,
                    float* tabTimingWeight,
                    float* matSqrtKeepSignalHex,
                    float* matTimingSumY,
                    float* matTimingSumXY,
                    float* matTimingSumX,
                    float* matTimingSumX2,
                    size_t nbPixel,
                    size_t nbEvent,
                    float* timingSlope,
                    float* timingIntercept)
{	
	
    
    std::for_each(EXECUTION_POLICY, std::begin(vecIndex), std::end(vecIndex),
        [=](int idx) {
            int eventIndex = idx / nbPixel;
			int pixelIndex = idx % nbPixel;
            
            
            // Compute sum
            size_t indexMat(eventIndex*nbPixel + pixelIndex);
            float scal(fabs(tabKeepSignalHex[indexMat]));
            
            float y(tabPosMax[indexMat]);
            matTimingSumY[indexMat] = y * scal;
           
            float longi = ((tabPosPixelX[pixelIndex] - hillasX[eventIndex])*tabCosPsi[eventIndex] + (tabPosPixelY[pixelIndex] - hillasY[eventIndex])*tabSinPsi[eventIndex]);
            
            matTimingSumXY[indexMat] = longi * y * scal;
            matTimingSumX[indexMat] = longi * scal;
            matTimingSumX2[indexMat] = longi * longi * scal;

        });
    
    // matrix vector multiplication to reduce
    VecData matTimingWeight_red = produitMatriceVecteur(matSqrtKeepSignalHex, tabone, vecEvtIdx.size());
    VecData matTimingSumX_red = produitMatriceVecteur(matTimingSumX, tabone, vecEvtIdx.size());
    VecData matTimingSumY_red = produitMatriceVecteur(matTimingSumY, tabone, vecEvtIdx.size());
    VecData matTimingSumX2_red = produitMatriceVecteur(matTimingSumX2, tabone, vecEvtIdx.size());
    VecData matTimingSumXY_red = produitMatriceVecteur(matTimingSumXY, tabone, vecEvtIdx.size());

    
    // macro pour execution policy
    std::for_each(std::execution::par_unseq, std::begin(vecEvtIdx), std::end(vecEvtIdx),
        [=](int eventIndex) {
        // Compute slope and intercept
        float scale(1.0f/matTimingWeight_red[eventIndex]);
        float sumX(matTimingSumX_red[eventIndex]*scale);
        float sumY(matTimingSumY_red[eventIndex]*scale);
        float sumX2(matTimingSumX2_red[eventIndex]*scale);
        float sig2 = sumX2 - sumX*sumX;
        float sumXY(matTimingSumXY_red[eventIndex]*scale);
        float slope((sumXY - sumX*sumY)/sig2);
        timingSlope[eventIndex] = slope;
        timingIntercept[eventIndex] = sumY - sumX*slope;
        
    
    });	
    
    
}




///Get the number of nanoseconds per elements of the Calibration
/** @param nbEvent : number of event of the tables
*/
void evaluateTiming(size_t nbEvent){
    // Let's define size of data:
    size_t nbPixel = NB_PIXEL;
    size_t nbElement = nbEvent * nbPixel;



    // Declare and initialize vectors for data
    VecData tabKeepSignalHex(nbElement, 42.0f);
    VecData tabPosMax(nbElement, 42.0f);
    VecData hillasX(nbEvent, 42.0f);
    VecData hillasY(nbEvent, 42.0f);
    VecData tabCosPsi(nbEvent, 42.0f);
    VecData tabSinPsi(nbEvent, 42.0f);
    VecData tabPosPixelX(nbElement, 42.0f);
    VecData tabPosPixelY(nbElement, 42.0f);
    VecData tabTimingWeight(nbEvent, 42.0f);
    VecData matSqrtKeepSignalHex(nbElement);
    VecData matTimingSumY(nbElement);
    VecData matTimingSumXY(nbElement);
    VecData matTimingSumX(nbElement);
    VecData matTimingSumX2(nbElement);
    VecData tabOnePixel(nbPixel, 42.0f);
    VecData timingSlope(nbEvent);
    VecData timingIntercept(nbEvent);

    // Initialize vecIdx
    VecIndex vecIdx(nbEvent*nbPixel), vecEvtIdx(nbEvent);
    std::iota(vecIdx.begin(), vecIdx.end(), 0);

    // Create pointers to data
    float* ptrTabKeepSignalHex = tabKeepSignalHex.data();
    float* ptrTabPosMax = tabPosMax.data();
    float* ptrHillasX = hillasX.data();
    float* ptrHillasY = hillasY.data();
    float* ptrTabCosPsi = tabCosPsi.data();
    float* ptrTabSinPsi = tabSinPsi.data();
    float* ptrTabPosPixelX = tabPosPixelX.data();
    float* ptrTabPosPixelY = tabPosPixelY.data();
    float* ptrTabTimingWeight = tabTimingWeight.data();
    float* ptrMatSqrtKeepSignalHex = matSqrtKeepSignalHex.data();
    float* ptrMatTimingSumY = matTimingSumY.data();
    float* ptrMatTimingSumXY = matTimingSumXY.data();
    float* ptrMatTimingSumX = matTimingSumX.data();
    float* ptrMatTimingSumX2 = matTimingSumX2.data();
    float* ptrTimingSlope = timingSlope.data();
    float* ptrTimingIntercept = timingIntercept.data();
    float* ptrtabone = tabOnePixel.data();

    // Call micro_benchmark
    size_t fullNbElement = nbElement; // Assuming you need fullNbElement elsewhere
    micro_benchmarkAutoNsPrint("evaluateCalibration 2d transform", fullNbElement, compute_timing,
                               vecIdx, vecEvtIdx, ptrtabone, ptrTabKeepSignalHex, ptrTabPosMax, ptrHillasX, ptrHillasY,
                               ptrTabCosPsi, ptrTabSinPsi, ptrTabPosPixelX, ptrTabPosPixelY,
                               ptrTabTimingWeight, ptrMatSqrtKeepSignalHex, ptrMatTimingSumY,
                               ptrMatTimingSumXY, ptrMatTimingSumX, ptrMatTimingSumX2,
                               nbPixel, nbEvent, ptrTimingSlope,
                               ptrTimingIntercept);
}




int main(int argc, char** argv){
	return micro_benchmarkParseArg(argc, argv, evaluateTiming);
}