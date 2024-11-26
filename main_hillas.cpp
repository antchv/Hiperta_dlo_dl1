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

// Fonction pour le produit matrice-matrice en utilisant BLAS avec des VecData
VecData produitMatriceMatrice(float* matriceA, float* matriceB, int m, int n, int p) {
    VecData resultat(m * p);

    // Effectuer le produit matrice-matrice en utilisant la fonction cblas_sgemm de BLAS (pour des données VecData)
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, p, n, 1.0f, matriceA, n, matriceB, p, 1.0f, resultat.data(), p);

    return resultat;
}







void compute_hillas(const VecIndex& vecEvtIndex, const VecIndex& vecIndex,
                    float* signal, float* projection,
                    float* hillasIntensity, float* hillasX, float* hillasY, float* hillasPhi,
                    float* hillasWidth, float* hillasLength,
                    float* hillasR, float* hillasPsi,
                    float* tabCosPsi, float* tabSinPsi,
                    size_t nbcol, size_t nbPixel)
{	
    //VecData matProjectedSignal = produitMatriceMatrice(signal, projection, vecEvtIndex.size(), nbPixel, nbcol);

    VecData matproj_sig = produitMatriceVecteur(signal, projection, vecEvtIndex.size());
    VecData matproj_x = produitMatriceVecteur(signal, projection, vecEvtIndex.size());
    VecData matproj_y = produitMatriceVecteur(signal, projection, vecEvtIndex.size());
    VecData matproj_x2 = produitMatriceVecteur(signal, projection, vecEvtIndex.size());
    VecData matproj_y2 = produitMatriceVecteur(signal, projection, vecEvtIndex.size());
    VecData matproj_xy = produitMatriceVecteur(signal, projection, vecEvtIndex.size());


    std::for_each(std::execution::par_unseq, std::begin(vecEvtIndex), std::end(vecEvtIndex),
        [=](int eventIndex) {
            /*
            float intensity = matProjectedSignal[eventIndex*nbcol];
            float invsumsig = (1.0f/intensity);

            float xm = matProjectedSignal[eventIndex*nbcol + 1lu] * invsumsig;
            hillasX[eventIndex] = xm;
            float ym = matProjectedSignal[eventIndex*nbcol + 2lu] * invsumsig;
            hillasY[eventIndex] = ym;
            float x2m = matProjectedSignal[eventIndex*nbcol + 3lu] * invsumsig;
            float y2m = matProjectedSignal[eventIndex*nbcol + 4lu] * invsumsig;
            float xym = matProjectedSignal[eventIndex*nbcol + 5lu] * invsumsig;
            */
            
            float intensity = matproj_sig[eventIndex];
            float invsumsig(1./intensity);
            
            float xm = matproj_x[eventIndex] * invsumsig;
            hillasX[eventIndex] = xm;
            float ym = matproj_y[eventIndex] * invsumsig;
            hillasY[eventIndex] = ym;
            float x2m = matproj_x2[eventIndex] * invsumsig;
            float y2m = matproj_y2[eventIndex] * invsumsig;
            float xym = matproj_xy[eventIndex] * invsumsig;


            float xm2(xm * xm);
            float ym2(ym * ym);
            float xmym(xm * ym);

            float vx2(x2m - xm2);
            float vy2(y2m - ym2);
            float vxy(xym - xmym);

            float d(vy2 + vx2);	
            hillasPhi[eventIndex] = atan2f(ym, xm);
            
            float z(sqrtf(d*d + 4.0f*(vxy*vxy - vy2*vx2)));	
            float eigenValueHigh = (vx2 + vy2 + z) / 2.0f;
            float eigenValueLow = (vy2 + vx2 - z) / 2.0f;
            
            hillasLength[eventIndex] = sqrtf(eigenValueHigh);	
            hillasWidth[eventIndex] = 2.0f*sqrtf(eigenValueLow);
            hillasR[eventIndex] = sqrt(xm2 + ym2);

            if(fabsf(vxy) < 1e-8){
                hillasPsi[eventIndex] = 0.0f;
            }else{
                float yOverX((eigenValueHigh - vx2)/vxy);
                hillasPsi[eventIndex] = atanf(yOverX);
            }
            tabCosPsi[eventIndex] = cosf(hillasPsi[eventIndex]);
            tabSinPsi[eventIndex] = sinf(hillasPsi[eventIndex]);
        }); 
    /*
    std::for_each(std::execution::seq, std::begin(vecIndex), std::end(vecIndex),
        [=](int idx) {
            int eventIndex = idx / nbPixel;
            int pixelIndex = idx % nbPixel;
            int currentpixel = eventIndex*nbPixel + pixelIndex;
            
            float longitudinal((tabPosPixelX[pixelIndex] - hillasX[eventIndex])*tabCosPsi[eventIndex] + (tabPosPixelY[pixelIndex] - hillasX[eventIndex])*tabSinPsi[eventIndex]);
	
            float l3(longitudinal*longitudinal*longitudinal);
            float a = tabSignalKept[currentpixel];
            matM3LongSignal[currentpixel] = l3*a;
            matM4LongSignal[currentpixel] = l3*longitudinal*a;
        });

    // Matrix vector multiplication to reduce on the pixel
    tabM3 = produitMatriceVecteur(matM3LongSignal, tabone, vecEvtIndex.size());
	tabM4 = produitMatriceVecteur(matM4LongSignal, tabone, vecEvtIndex.size());

    std::for_each(std::execution::seq, std::begin(vecEvtIndex), std::end(vecEvtIndex),
        [=](int eventIndex) {
            float length(hillasLength[eventIndex]), intensity(hillasIntensity[eventIndex]);
            float length3sumSig(length*length*length*intensity);
            hillasSkewness[eventIndex] = tabM3[eventIndex]/(length3sumSig);
            hillasKurtosis[eventIndex] = tabM4[eventIndex]/(length3sumSig*length);
            hillasLength[eventIndex] *= 2.0f;

        });

    */
    
        	
}




///Get the number of nanoseconds per elements of the Calibration
/** @param nbEvent : number of event of the tables
*/
void evaluateHillas(size_t nbEvent){
    // Let's define size of data:
    size_t nbElement(nbEvent * NB_PIXEL);
    size_t nbcol(6);

	// Le signal
	VecData vecSignal(nbElement);
	std::fill(vecSignal.begin(), vecSignal.end(), 42.0f);

    VecData xm(nbEvent), ym(nbEvent), intensity(nbEvent), length(nbEvent), width(nbEvent), r(nbEvent), phi(nbEvent), psi(nbEvent), sinpsi(nbEvent), cospsi(nbEvent);


	// Le vecteur d'indices pour parcourir tout les pixels et les events
	VecIndex vecIdx(NB_PIXEL*nbEvent), vecEvtIndex(nbEvent);
	std::iota(vecIdx.begin(), vecIdx.end(), 0);
	std::iota(vecEvtIndex.begin(), vecEvtIndex.end(), 0);


    VecData posPixelX(nbElement), posPixelY(nbElement), projection(nbEvent*NB_PIXEL);
	
	// We have to create pointer to be able to catch them by copy without losing any time
	float *tabSignal = vecSignal.data(), *tabxm = xm.data(), *tabym = ym.data(), *tabintensity = intensity.data(), *tablength = length.data(), *tabwidth = width.data(), *tabr = r.data(), *tabphi = phi.data(), *tabpsi = psi.data(), *tabCosPsi = cospsi.data(), *tabSinPsi = sinpsi.data(), *tabprojection = projection.data();

	// Appel de microBenchMark pour tester la performance de la fonction de calibration
	size_t fullNbElement(nbElement);

	micro_benchmarkAutoNsPrint("evaluateCalibration 2d tranform", fullNbElement, compute_hillas, vecEvtIndex, vecIdx, tabSignal, tabprojection, tabintensity, tabxm, tabym, tabphi, tabwidth, tablength, tabr, tabpsi, tabCosPsi, tabSinPsi, nbcol, NB_PIXEL);
	
}

int main(int argc, char** argv){
	return micro_benchmarkParseArg(argc, argv, evaluateHillas);
}