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
void compute_hillas(float* signal, float* sumsig,
                    float* xm, float* ym, float* x2m, float* y2m, float* xym,
                    const VecIndex& vecEvtIndex,
                    float* posPixelX, float* posPixelY, float* tempXY,
                    float* phi, float* length, float* width, float* r, float* psi, float* longitudinal, float* longitudinal3,
                    float* m3_long, float* m4_long, float* skewness, float* kurtosis,
                    size_t nbPixel)
{	// Itération sur tout les events
    std::for_each(std::execution::par_unseq, std::begin(vecEvtIndex), std::end(vecEvtIndex),
        [=](int evtidx) {

            // Définition de l'indice du premier pixel de l'event courant
            long firstPixelIdx = evtidx * nbPixel;

            // Somme du signal sur tout les pixels de l'event
            sumsig[evtidx] = std::transform_reduce(std::execution::unseq, signal + firstPixelIdx, signal + firstPixelIdx +   nbPixel, 0.0f, std::plus{},
                            [](auto val) {return val;});

            // Calcul des moments xm, ym, x²m, y²m, xym
            xm[evtidx] = std::transform_reduce(std::execution::unseq, signal + firstPixelIdx, signal + firstPixelIdx + nbPixel, posPixelX + firstPixelIdx, 0.0f, std::plus{},
                    [](auto val, float x) {return val*x;});
            
            ym[evtidx] = std::transform_reduce(std::execution::unseq, signal + firstPixelIdx, signal + firstPixelIdx + nbPixel, posPixelY + firstPixelIdx, 0.0f, std::plus{},
                    [](auto val, float y) {return val*y;});

            x2m[evtidx] = std::transform_reduce(std::execution::unseq, signal + firstPixelIdx, signal + firstPixelIdx + nbPixel, posPixelX + firstPixelIdx, 0.0f, std::plus{},
                    [](auto val, float x) {return val*x*x;});
            
            y2m[evtidx] = std::transform_reduce(std::execution::unseq, signal + firstPixelIdx, signal + firstPixelIdx + nbPixel, posPixelY + firstPixelIdx, 0.0f, std::plus{},
                    [](auto val, float y) {return val*y*y;});

            // Utilisation d'un temporaire xy pour le calcul de xym
            std::transform(std::execution::unseq, posPixelX + firstPixelIdx, posPixelX + firstPixelIdx + nbPixel, posPixelY + firstPixelIdx, tempXY + firstPixelIdx,
                [=](float x, float y) {
                    return x*y;
                });
            
            xym[evtidx] = std::transform_reduce(std::execution::unseq, signal + firstPixelIdx, signal + firstPixelIdx + nbPixel, tempXY + firstPixelIdx, 0.0f, std::plus{},
                    [](auto val, float xy) {return val*xy;});


            
            float invsumsig(1./sumsig[evtidx]);

            xm[evtidx] *= invsumsig;
            ym[evtidx]*= invsumsig;

            x2m[evtidx] *= invsumsig;
            y2m[evtidx] *= invsumsig;
            
            xym[evtidx] *= invsumsig;


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
                
                std::transform(std::execution::unseq, posPixelX + firstPixelIdx, posPixelX + firstPixelIdx + nbPixel, posPixelY + firstPixelIdx, longitudinal,
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

    VecData longitudinal(NB_PIXEL), longitudinal3(NB_PIXEL);

    VecData xm(nbEvent), ym(nbEvent), x2m(nbEvent), y2m(nbEvent), xym(nbEvent), intensity(nbEvent), length(nbEvent), width(nbEvent), phi(nbEvent), r(nbEvent), psi(nbEvent), m3_long(nbEvent), m4_long(nbEvent), skewness(nbEvent), kurtosis(nbEvent);

	// Le signal
	VecData vecSignal(nbElement);
	std::fill(vecSignal.begin(), vecSignal.end(), 42.0f);


	// Le vecteur d'indices pour parcourir tout les pixels et les events
	VecIndex vecEvtIdx(nbEvent);
	std::iota(vecEvtIdx.begin(), vecEvtIdx.end(), 0);

    VecData posPixelX(nbElement), posPixelY(nbElement), tempXY(nbElement);


	
	// We have to create pointer to be able to catch them by copy without losing any time
	float *tabSignal = vecSignal.data(), *tabxm = xm.data(), *tabym = ym.data(), *tabx2m = x2m.data(), *taby2m = y2m.data(), *tabpospixelx = posPixelX.data(), *tabpospixely = posPixelY.data(), *tabintensity = intensity.data(), *tabtempXY = tempXY.data(), *tabxym = xym.data(), *tabphi = phi.data(), *tablength = length.data(), *tabwidth = width.data(), *tabr = r.data(), *tabpsi = psi.data(), *tablongitudinal = longitudinal.data(), *tablongitudinal3 = longitudinal3.data(), *tabm3 = m3_long.data(), *tabm4 = m4_long.data(), *tabskewness = skewness.data(), *tabkurtosis = kurtosis.data();

	size_t fullNbElement(nbElement);

    micro_benchmarkAutoNsPrint("evaluateCalibration 2d tranform", fullNbElement, compute_hillas, tabSignal, tabintensity, tabxm, tabym, tabx2m, taby2m, tabxym, vecEvtIdx, tabpospixelx, tabpospixely, tabtempXY, tabphi, tablength, tabwidth, tabr, tabpsi, tablongitudinal, tablongitudinal3, tabm3, tabm4, tabskewness, tabkurtosis, NB_PIXEL);
}

int main(int argc, char** argv){
	return micro_benchmarkParseArg(argc, argv, evaluateHillas);
}