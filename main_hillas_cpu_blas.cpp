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
VecData produitMatriceVecteur(float* matrice, float* vecteur, int m, int n) {
    VecData resultat(m);
    cblas_sgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0f, matrice, n, vecteur, 1, 1.0f, resultat.data(), 1);

    return resultat;
}

// Fonction pour le produit matrice-matrice en utilisant BLAS avec des VecData
VecData produitMatriceMatrice(float* matriceA, float* matriceB, int m, int n, int p) {
    VecData resultat(m * p);

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, p, n, 1.0f, matriceA, n, matriceB, p, 1.0f, resultat.data(), p);

    return resultat;
}


///Make the Hillas reconstruction with a telescope
/**	@param signal : signal in the pixels telescope
 * 	@param[out] sumsig : Total aplitude of the signal in all pixels of the telescope
 *  @param[out] xm : X coordinate of the center of gravity of the shower in the camera
 * 	@param[out] ym : Y coordinate of the center of gravity of the shower in the camera
 * 	@param[out] phi : Orientation angle of Major Axis w.r.t. x-axis
 * 	@param[out] width : telescope image width
 * 	@param[out] length : telescope image length
 * 	@param[out] skewness : Assymetry of the image
 * 	@param[out] r : distance from the camera center
 * 	@param[out] kurtosis : Kurtosis excess, define the peakedness of the distribution
 * 	@param[out] psi : direction, angle of the shower
 * 	@param posPixelX : position of the pixels
 * 	@param posPixelY : position of the pixels
 *  @param tempXY : temporary tab to compute xym
 * 	@param nbPixel : number of pixels in the tables
 */
void compute_hillas(float* vecsignal, float* projection, float* signal, float* sumsig,
                    float* xm, float* ym, float* x2m, float* y2m, float* xym,
                    const VecIndex& vecEvtIndex,
                    float* posPixelX, float* posPixelY, float* tempXY,
                    float* phi, float* length, float* width, float* r, float* psi, float* longitudinal, float* longitudinal3,
                    float* m3_long, float* m4_long, float* skewness, float* kurtosis,
                    size_t nbPixel, size_t nbcol)
{   
    
    VecData matproj_sig = produitMatriceVecteur(vecsignal, projection, vecEvtIndex.size(), nbcol);
    VecData matproj_x = produitMatriceVecteur(vecsignal, projection, vecEvtIndex.size(), nbcol);
    VecData matproj_y = produitMatriceVecteur(vecsignal, projection, vecEvtIndex.size(), nbcol);
    VecData matproj_x2 = produitMatriceVecteur(vecsignal, projection, vecEvtIndex.size(), nbcol);
    VecData matproj_y2 = produitMatriceVecteur(vecsignal, projection, vecEvtIndex.size(), nbcol);
    VecData matproj_xy = produitMatriceVecteur(vecsignal, projection, vecEvtIndex.size(), nbcol);
    
    
    //VecData matproj = produitMatriceMatrice(vecsignal, projection, vecEvtIndex.size(), nbPixel, nbcol);

    // Itération sur tout les events
    std::for_each(std::execution::unseq, std::begin(vecEvtIndex), std::end(vecEvtIndex),
        [=](int evtidx) {    
            // Définition de l'indice du premier pixel de l'event courant
            long firstPixelIdx = evtidx * nbPixel;
            
            sumsig[evtidx] = matproj_sig[evtidx];
            
            float invsumsig(1./sumsig[evtidx]);
            
            // Calcul des moments xm, ym, x²m, y²m, xym
            xm[evtidx] = matproj_x[evtidx] * invsumsig;
            
            ym[evtidx] = matproj_y[evtidx] * invsumsig;

            x2m[evtidx] = matproj_x2[evtidx] * invsumsig;
            
            y2m[evtidx] = matproj_y2[evtidx] * invsumsig;
            
            xym[evtidx] = matproj_xy[evtidx] * invsumsig;
            
            /*
            sumsig[evtidx] = matproj[evtidx * nbcol];
            
            float invsumsig(1./sumsig[evtidx]);
            
            // Calcul des moments xm, ym, x²m, y²m, xym
            xm[evtidx] = matproj[evtidx * nbcol + 1] * invsumsig;
            
            ym[evtidx] = matproj[evtidx * nbcol + 2] * invsumsig;

            x2m[evtidx] = matproj[evtidx * nbcol + 3] * invsumsig;
            
            y2m[evtidx] = matproj[evtidx * nbcol + 4] * invsumsig;
            
            xym[evtidx] = matproj[evtidx * nbcol + 5] * invsumsig;
            
            */
            
            // Définition de grandeurs utiles pour le calcul des autres paramètres
            float xm2(xm[evtidx] * xm[evtidx]);
            float ym2(ym[evtidx] * ym[evtidx]);
            float xmym(xm[evtidx] * ym[evtidx]);

            float vx2(x2m[evtidx] - xm2);
            float vy2(y2m[evtidx] - ym2);
            float vxy(xym[evtidx] - xmym);

            float d(vy2 + vx2);	
            
            // Orientation de l'ellispe, angle entre la position de l'ellispe et l'axe des abscisses
            phi[evtidx] = atan2f(ym[evtidx], xm[evtidx]);
            
            float z(sqrtf(d*d + 4.0f*(vxy*vxy - vy2*vx2)));
            float eigenValueHigh = (vx2 + vy2 + z) / 2.0f;
            float eigenValueLow = (vy2 + vx2 - z) / 2.0f;

            // Longueur et largeur de l'ellispe
            length[evtidx] = sqrtf(eigenValueHigh);

            width[evtidx] = sqrtf(eigenValueLow);
            
            // Distance au centre
            r[evtidx] = sqrt(xm2 + ym2);


            // Direction de l'ellipse
            if(fabsf(vxy) < 1e-8){
                psi[evtidx] = 0.0f;
            }else{
                float yOverX((eigenValueHigh - vx2)/vxy);
                psi[evtidx] = atanf(yOverX);
            }

            
            // Calcul de skewness et kurtosis
            if(fabs(length[evtidx]) > 0.0){
                float cosPsi(cosf(psi[evtidx])), sinPsi(sinf(psi[evtidx]));
                float xm_value = xm[evtidx];
                float ym_value = ym[evtidx];
                
                std::transform(std::execution::unseq, posPixelX, posPixelX + nbPixel, posPixelY, longitudinal,
                    [&](float x, float y) {
                        return (x - xm_value)*cosPsi + (y - ym_value)*sinPsi;
                    });
                
                    


                m3_long[evtidx] = std::transform_reduce(std::execution::unseq, signal + firstPixelIdx, signal + firstPixelIdx + nbPixel, longitudinal, 0.0f, std::plus{},
                [](auto a, float l) {return l*l*l*a;});



                m4_long[evtidx] = std::transform_reduce(std::execution::unseq, signal + firstPixelIdx, signal + firstPixelIdx + nbPixel, longitudinal, 0.0f, std::plus{},
                [](auto a, float l) {return l*l*l*l*a;});
                
                float length3sumSig(length[evtidx]*length[evtidx]*length[evtidx]*sumsig[evtidx]);
                skewness[evtidx] = m3_long[evtidx]/(length3sumSig);
                kurtosis[evtidx] = m4_long[evtidx]/(length3sumSig*length[evtidx]);
            }

            length[evtidx] *= 2.0f;
            width[evtidx] *= 2.0f;

            

        });	
}







///Get the number of nanoseconds per elements of the Calibration
/** @param nbEvent : number of event of the tables
*/
void evaluateHillas(size_t nbEvent){
    // Let's define size of data:
    size_t nbElement(nbEvent * NB_PIXEL);
    size_t nbcol(1);

    VecData longitudinal(NB_PIXEL), longitudinal3(NB_PIXEL);

    VecData xm(nbEvent), ym(nbEvent), x2m(nbEvent), y2m(nbEvent), xym(nbEvent), intensity(nbEvent), length(nbEvent), width(nbEvent), phi(nbEvent), r(nbEvent), psi(nbEvent), m3_long(nbEvent), m4_long(nbEvent), skewness(nbEvent), kurtosis(nbEvent);

	// Le signal
	VecData vecSignal(nbElement), projection(nbcol * NB_PIXEL);
	std::fill(vecSignal.begin(), vecSignal.end(), 42.0f);


	// Le vecteur d'indices pour parcourir tout les pixels et les events
	VecIndex vecEvtIdx(nbEvent);
	std::iota(vecEvtIdx.begin(), vecEvtIdx.end(), 0);

    VecData posPixelX(NB_PIXEL), posPixelY(NB_PIXEL), tempXY(NB_PIXEL);


	
	// We have to create pointer to be able to catch them by copy without losing any time
	float *tabSignal = vecSignal.data(), *tabxm = xm.data(), *tabym = ym.data(), *tabx2m = x2m.data(), *taby2m = y2m.data(), *tabpospixelx = posPixelX.data(), *tabpospixely = posPixelY.data(), *tabintensity = intensity.data(), *tabtempXY = tempXY.data(), *tabxym = xym.data(), *tabphi = phi.data(), *tablength = length.data(), *tabwidth = width.data(), *tabr = r.data(), *tabpsi = psi.data(), *tablongitudinal = longitudinal.data(), *tablongitudinal3 = longitudinal3.data(), *tabproj = projection.data(), *tabm3 = m3_long.data(), *tabm4 = m4_long.data(), *tabskewness = skewness.data(), *tabkurtosis = kurtosis.data();

	size_t fullNbElement(nbElement);

    micro_benchmarkAutoNsPrint("evaluateCalibration 2d tranform", fullNbElement, compute_hillas, tabSignal, tabproj, tabSignal, tabintensity, tabxm, tabym, tabx2m, taby2m, tabxym, vecEvtIdx, tabpospixelx, tabpospixely, tabtempXY, tabphi, tablength, tabwidth, tabr, tabpsi, tablongitudinal, tablongitudinal3, tabm3, tabm4, tabskewness, tabkurtosis, NB_PIXEL, nbcol);
}

int main(int argc, char** argv){
	return micro_benchmarkParseArg(argc, argv, evaluateHillas);
}