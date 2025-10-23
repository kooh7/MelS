//
//  MelSpecD.m
//  MelS
//
//  Created by Ken on 10/10/2025.
//

#import "MelSpecD.h"
#import <math.h>
#import <stdlib.h>
#import <float.h>
#import <Accelerate/Accelerate.h>
#import <CoreML/CoreML.h>
#import <Foundation/Foundation.h>



@implementation MelSpecD

- (instancetype)init{
    [self initWindow];
    [self initDFT];
    [self initMelWindows];
    
    return self;
}

- (void)initWindow{  // hanning window
    double ds = (double)winSize - 1.0;
    for( double i = 0.0; i < winSize; i+=1.0){
        sigWin[(int)i] = pow(sin( fmod(M_PI * i / ds, M_PI)), 2.0);
    }
}

- (void)initDFT{
    for ( int k = 0; k < dftSize; k++){
        dftmatrix[k] = 1.0f;
        if (k%2 == 0){
            dftmatrix[((dftSize-1) *dftSize) + k] = 1.0;
        }
        else{
            dftmatrix[((dftSize-1) *dftSize) + k] = -1.0;
        }
    }
    
    for (int i = 1; i < (dftSize/2); i++){
        for ( int k = 0; k < dftSize; k++){
            double angle = (2 * i * k) * M_PI / dftSize;
            dftmatrix[(((i*2)-1) *dftSize)  + k] = cos(angle);
            dftmatrix[(((i*2)) *dftSize)  + k] = sin(angle);
        }
    }
}


- (double)melToHz:(double)melFreq {
    return 700 * (exp(melFreq / 1127.0) - 1.0);
}

- (double)HzToMel:(double)hzf{
    return 1127.0 * log(1.0 + hzf / 700.0);
}

- (float)HzToMelF:(float)hzf{
    return 1127.0f * logf(1.0f + hzf / 700.0f);
}


-(void)initMelWindows{
    
  
    double melMin = [self HzToMel:minFreq];  //1127.0 * log(1.0 + minFreq / 700.0);
    double melMax = [self HzToMel:maxFreq];  //1127.0 * log(1.0 + maxFreq / 700.0);
    double melStep = (melMax - melMin) / (nMelWins + 1.0);
    double fftbinWidth = (fs / dftSize);
    
    /*
    double cmelFreqs[nMelWins], lmelFreqs[nMelWins], rmelFreqs[nMelWins];
    double melFs[dftCoeffSize];
    
   
    float melMin = [self HzToMelF:(float)minFreq];
    float melMax = [self HzToMelF:(float)maxFreq];
    float melStep = (melMax - melMin) / (nMelWins + 1.0f);
    float fftbinWidth = fs / dftSize;
     */
     
    //float cmelFreqs[nMelWins], lmelFreqs[nMelWins], rmelFreqs[nMelWins];
    float melFreqs[nMelWins+2];
    float melFs[dftCoeffSize];
     
    
    
    
    
    for (int i = 0; i < nMelWins + 2; i++) {
        //melFreqs[i] = (float)(melMin) + (i * (float)(melStep));        // 6.3 -9    000996     // melscale in increments
        //melFreqs[i] = (float)(melMin) + (float)(i * melStep);          // 8.5 -9      00128
        melFreqs[i] = (float)(melMin + (i * melStep));                 // 3.3 -9      00054
        //melFreqs[i] = (melMin + (i * melStep));
        //printf("mel %.20f   : ix  %d \n", melFreqs[i], i  );
    }
    
    /*
    for (int i = 0; i < nMelWins; i++) {
        lmelFreqs[i] = (float)(melMin + (i * melStep));                 // 3.3 -9      00054
        cmelFreqs[i] = (float)(melMin + ((i + 1.0) * melStep));
        rmelFreqs[i] = (float)(melMin + ((i + 2.0) * melStep));
        //printf("mel %.20f   : ix  %d \n", melFreqs[i], i  );
    }
     */
    
    /*
    for (int i = 0; i < nMelWins; i++) {
        lmelFreqs[i] = (float)(melMin + (i * melStep));                 // 3.3 -9      00054
        cmelFreqs[i] = (float)(melMin + ((i + 1.0) * melStep));
        rmelFreqs[i] = (float)(melMin + ((i + 2.0) * melStep));
        //printf("mel %.20f   : ix  %d \n", melFreqs[i], i  );
    }
     */
    
    
    for(int j = 0; j < dftCoeffSize-1; j++){
        float fr = (float)(j * fftbinWidth);
        melFs[j] = (float)[self HzToMelF:(float)fr];
        
        //double fr = j * fftbinWidth;
        //melFs[j] = [self HzToMel:fr];
        //printf("MelF %d    %.30lf \n", j, melFs[j]);
    }
    melFs[dftCoeffSize-1] = 0.0f;
    
    for( int i=0; i < nMelWins; i++){
        for(int j = 0; j < dftCoeffSize; j++){
            float melFr = melFs[j];
            float fval = 0.0f;
            
            if (melFr > melFreqs[i] && melFr < melFreqs[i+1]){
                fval = (melFr - melFreqs[i]) / (melFreqs[i+1] - melFreqs[i]);
            } else if (melFr >= melFreqs[i+1]  && melFr < melFreqs[i+2] ) {
                
                fval = (melFreqs[i+2] - melFr) / (melFreqs[i+2] - melFreqs[i+1]);
            }
            
            melWindows[i*dftCoeffSize + j] = fval >= 0.0f  ? fval : 0.0f;
        }
    }
}
    



-(void)setSignal:(float*)isig :(int)sigl{
    sig = isig;
    siglen = sigl;
    printf("Signal set %d \n", siglen);
    [self windowSignal];
    printf("Window signalled \n");
    [self calculateSpectro];
    [self calculateMelSpec];
}




- (void)windowSignal{
    
    nframes = 2 + (int)floor( (siglen - dftSize)/(hopSize) );
    windowedSignal = (float**)malloc(sizeof(float*) * nframes);
    windowedSig = (float**)malloc(sizeof(float*) * nframes);
    
    
    float coltot, colmean;
    int olap = (dftSize - winSize) / 2;
    int type = 1;  // 1 is the skewed window    // 2 is the centred window
    
    // Fill the WindowedSignal and take away the mean
    for(int i = 0; i < nframes; i++){
        
        windowedSignal[i] = (float*)calloc(dftSize, sizeof(float));
        windowedSig[i] = (float*)calloc(dftSize, sizeof(float));
        
        int wsize;
        int instix = hopSize * i;
        coltot = 0.0f;
        
        if (siglen - instix < winSize){
            wsize = siglen - instix;
        }else{
            wsize = winSize;
        }
        
        if ( type == 1){
            for ( int j = 0; j < wsize; j++){
                windowedSignal[i][j] = sig[instix + j];
                coltot += sig[j + instix];
                windowedSig[i][j] = 0.0f;
            }
            for ( int j = wsize; j < dftSize; j++){
                windowedSignal[i][j] = 0.0f;
                windowedSig[i][j] = 0.0f;
            }
            colmean = coltot / winSize;
            for (int j = 0; j < winSize; j++){
                windowedSignal[i][j] -= colmean;
            }
            
            //for (int col = 0; col < nframes; col++) {
                //windowedSig[col][olap] = (1.0-preCoef) * windowedSignal[col][olap];
                //for (int row = olap+1; row < winSize+olap; row++) {
            windowedSig[i][0] = (1.0-preCoef) * windowedSignal[i][0];
            for (int row = 1; row < winSize; row++){
                windowedSig[i][row] = windowedSignal[i][row] - (preCoef * windowedSignal[i][row-1]);
            }
            
        }else{
            for (int j=0; j < olap; j++){
                //windowedSignal[i][j] = 0.0f;
                //windowedSig[i][j] = 0.0f;
            }
            for ( int j = 0; j < wsize; j++){
                windowedSignal[i][j+olap] = sig[instix + j];
                coltot += sig[j + instix];
                windowedSig[i][j+olap] = 0.0f;
            }
            for (int p = wsize+olap; p < dftSize; p++){
                windowedSignal[i][p] = 0.0f;
                windowedSig[i][p] = 0.0f;
            }
            
            colmean = coltot / winSize;
            for (int j = olap; j < dftSize-olap; j++){
                windowedSignal[i][j] -= colmean;
            }
            
            windowedSig[i][olap] = (float)((1.0-preCoef) * windowedSignal[i][olap]);
            for (int row = olap+1; row < winSize+olap; row++) {
                windowedSig[i][row] = (float)(windowedSignal[i][row] - (preCoef * windowedSignal[i][row-1]));
            }
        }
    }

    for( int i = 0; i < nframes; i++){
        free(windowedSignal[i]);
        windowedSignal[i] = NULL;
    }
    free(windowedSignal);
    windowedSignal = NULL;
        
    
    for (int i = 0; i < winSize; i++){
        for (int j = 0; j < nframes; j++){
            if (type == 1){
                windowedSig[j][i] *= (float)sigWin[i];
            }else{
                windowedSig[j][i + olap] *= (float)sigWin[i];
            }
        }
    }
}



- (void)calculateSpectro{
    
    
    spectro = malloc(nframes * sizeof(double*));
    for( int i = 0; i < nframes; i++){
        spectro[i] = calloc(dftCoeffSize, sizeof(double) );
    }
    
    int type = 1;   // 1 : DFT   2 : FFT
    
    //if (type == 2){
    int l2dft = (int)log2(dftSize);
    float *realpart = (float *)malloc((dftSize/2) * sizeof(float));
    float *imagpart = (float *)malloc((dftSize/2) * sizeof(float));
    DSPSplitComplex outputComplex = {realpart, imagpart};
    FFTSetup fftSetup = vDSP_create_fftsetup(l2dft, kFFTRadix2);
    //}
    
    double* specptr;
    float* sigptr;
    
    double oddcoeff, evencoeff;
    int oddix, evenix;
   
    
    for( int i = 0; i < nframes; i++){
        
        specptr = spectro[i];
        sigptr = windowedSig[i];
        oddcoeff = 0.0;
        evencoeff = 0.0;
        //int oddix, evenix;
        
        if (type == 1){
            
            // DFT version
            for (int j = 0; j < dftSize; j++){
                oddcoeff += (dftmatrix[j] * sigptr[j]);
                evencoeff += (dftmatrix[((dftSize-1) * dftSize) + j] * sigptr[j]);
            }
            //if (oddcoeff < 0.0f){
            //    oddcoeff = -oddcoeff;
            //}
            //if (evencoeff < 0.0f){
            //    evencoeff = - evencoeff;
            //}
            
            specptr[0] = oddcoeff*oddcoeff;
            specptr[dftCoeffSize-1] = evencoeff*evencoeff;
            
            for(int i = 1; i < dftCoeffSize-1; i++){
                oddix = (i*2) - 1;
                evenix = i*2;
                oddcoeff = 0.0;
                evencoeff = 0.0;
                
                for (int j = 0; j < dftSize; j++){
                    oddcoeff += (dftmatrix[(oddix*dftSize) + j] * (double)sigptr[j]);
                    evencoeff += (dftmatrix[(evenix*dftSize)+ j] * (double)sigptr[j]);
                }
                specptr[i] = ((oddcoeff*oddcoeff) + (evencoeff * evencoeff ));
            }
        }else{
            
            //printf("CAlling FFT %d \n", i);
            
            // Perform the FFT - Use vDSP_fft_r2c for real-to-complex FFT
            for (int j = 0; j < dftSize/2; j++) {
                outputComplex.realp[j] = sigptr[j*2] ;
                outputComplex.imagp[j] = sigptr[(j*2)+1]; // Imaginary part is zero for real input
            }
            //printf("Assigned ");
            
            vDSP_fft_zrip(fftSetup, &outputComplex, 1, l2dft, FFT_FORWARD);
                        
            //printf("FFT  ");
            float dfts2 = dftSize * dftSize;
            float dftsq = sqrt(dftSize);
            for (int j = 0; j < dftCoeffSize; j++) { // N/2 because we're only interested in the positive frequencies
                spectro[i][j] = 0.25 * ( ( outputComplex.realp[j] * outputComplex.realp[j]) + (outputComplex.imagp[j] * outputComplex.imagp[j]) );
            }
            //printf("Mag  ");
        }
    }
    free(realpart);
    free(imagpart);
    
    vDSP_destroy_fftsetup(fftSetup);
}





-(void)calculateMelSpec{
    
    //melSpectro = malloc(sizeof(double*) * maxLength);
    double stdsq = 4.5689974 * 2.0;
    double smean = -4.2677393;
    
    //melSpectro[i] = calloc(nMelWins, sizeof(double) );
    memset(melSpectro, 0, sizeof(melSpectro));
    
    double mspectrum[128];
    
    for (int i = 0; i < maxLength; i++){      // sig cols
        
        memset(mspectrum, 0, sizeof(mspectrum));
        
        if (i < nframes){
            //double totval = 0.0;
            for( int j = 0; j < nMelWins; j++){
                
                for (int k = 0; k < dftCoeffSize; k++){
                    //melSpectro[j][i] += spectro[k][i] * melWindows[k][j];
                    mspectrum[j] += (spectro[i][k] * (double)melWindows[(j*dftCoeffSize) + k]);
                }
                if (mspectrum[j] < (double)FLT_EPSILON){
                    mspectrum[j] = log((double)(FLT_EPSILON));
                }
                else{
                    mspectrum[j] = log(mspectrum[j]);
                }
            }
        }
        for( int j = 0; j < nMelWins; j++){
            mspectrum[j] -= smean;
            mspectrum[j] /= stdsq;
            melSpectro[(i*nMelWins) + j] = (float)mspectrum[j];
        }
    }
}




    


-(void) dealloc{
    
    if(windowedSignal != NULL) {
        for( int i = 0; i < nframes; i++){
            if (windowedSignal[i] != NULL) {
                free(windowedSignal[i]);
                windowedSignal[i] = NULL;
            }
        }
        free(windowedSignal);
        windowedSignal = NULL;
    }
    
    if(windowedSig != NULL) {
        for( int i = 0; i < nframes; i++){
            if (windowedSig[i] != NULL) {
                //printf("Freeing windowedSig[%d]: %p\n", i, windowedSig[i]);
                free(windowedSig[i]);
                windowedSig[i] = NULL;
            }
        }
        free(windowedSig);
        windowedSig = NULL;
    }
    
  
    
    if(spectro != NULL) {
        for( int i = 0; i < nframes; i++){
            if (spectro[i] != NULL) {
                free(spectro[i]);
                spectro[i] = NULL;
            }
        }
        free(spectro);
    }
    
    /*
    if(melSpectro != NULL) {
        for( int i = 0; i < nframes; i++){
            if (melSpectro[i] != NULL) {
                free(melSpectro[i]);
                melSpectro[i] = NULL;
            }
        }
        free(melSpectro);
        melSpectro = NULL;
    }
    */
}


-(void)convertToMLMultiArray{
        
    NSArray<NSNumber *> *shape = @[@1, @maxLength, @nMelWins];  // Shape: (1, 1024, 128)
        
    NSError *error = nil;
    MLMultiArray *mlArray = [[MLMultiArray alloc] initWithShape:shape dataType:MLMultiArrayDataTypeFloat32 error:&error];
        
    if (error) {
        NSLog(@"Error creating MLMultiArray: %@", error);
        return;
    }
        
    for (NSInteger i = 0; i < 131072; i++) {
        mlArray[i] = @(melSpectro[i]);
    }
        
    // Now `mlArray` has the shape (1, 1024, 128)
    NSLog(@"MLMultiArray: %@", mlArray);
}



// FUNCTIONS USED TO RETREIVE VARIABLES - some only used for testing


-(double*)sigwin{
    return sigWin;
}

-(int)getMaxLength{
    return maxLength;
}

- (MLMultiArray*)getMLarray{
    return mlArray;
}

-(int)winsize{
    return winSize;
}

-(int) getdftsize{
    return dftSize;
}

- (float*)melwin{
    return &melWindows[0];
}

- (int)melwinsize{
    return nMelWins;
}

- (int)coefsize{
    return dftCoeffSize;
}

-(float**)getWindowedSig{
    printf("Returning WindowedSig");
    return windowedSig;
}

- (int)getNumFrames{
    return nframes;
}


- (double**)getSpectro{
    return spectro;
}

- (float*)getMelSpectro{
    return melSpectro;
}

-(double*) getdftMat{
    return &dftmatrix[0];
}


@end








/*
float upslope = ((melFr - lmelFreqs[i]) / (cmelFreqs[i] - lmelFreqs[i]));
float downslope = ((rmelFreqs[i] - melFr) / (rmelFreqs[i] - cmelFreqs[i]));

//float upslope = ((melFr - melFreqs[i]) / (melFreqs[i+1] - melFreqs[i]));
//float downslope = ((melFreqs[i+2] - melFr) / (melFreqs[i+2] - melFreqs[i+1]));
//float upslope = (float)(melFr - melFreqs[i]) / melStep;
//float downslope = (float)(melFreqs[i+2] - melFr) / melStep;
fval = fminf(upslope, downslope);
fval = fmaxf(fval, 0.0f);
melWindows[i*dftCoeffSize + j] = fval;
*/



/*

- (void)initMelWindows{
    double melMin = 1127.0 * log(1.0 + minFreq / 700.0);
    double melMax = 1127.0 * log(1.0 + maxFreq / 700.0);
    double melStep = (melMax - melMin) / (nMelWins + 1.0);
    
    float melFreqs[nMelWins + 2];  //, melHz[nMelWins + 2];
    float melFs[dftCoeffSize];
    
    double fftbinWidth = fs / dftSize;
    double freq, val, melF;
    
    printf("MelMin %.30lf \n", melMin);
    printf("MelMax %.30lf \n", melMax);
    printf("MelStep %.30lf \n", melStep);
    printf("FFTBW %.30lf \n", fftbinWidth);
    
    for (int i = 0; i < nMelWins + 2; i++) {
        melFreqs[i] = melMin + (float)(i * melStep);         // melscale in increments of melstep
        //melHz[i] = [self melToHz: melFreqs[i]];     // the melscale in Hz ( with diff index )
        printf("mel %.20f   : ix  %d \n", melFreqs[i], i  );
    }
    
    for(int j = 0; j < dftCoeffSize; j++){
        float fr = j * (float)fftbinWidth;
        melFs[j] = [self HzToMelF:fr];
        printf("MelF %d    %.30lf \n", j, melFs[j]);
    }

    //melWindows = malloc(sizeof(double*) * nMelWins);
    for( int i=0; i < nMelWins; i++){
        //melWindows[i] = malloc(sizeof(double) * dftCoeffSize);
        
        for(int j = 0; j < dftCoeffSize; j++){
            //freq = j * fftbinWidth;
            //melF = 1127.0 * log(1.0 + freq / 700.0);           // melF is the mel at each point
            //melF = [self HzToMel:freq];
            
            //float fr = (float)(j * fftbinWidth);
            //float melFr = [self HzToMelF:fr];
            float melFr = melFs[j];
            //printf("MelF %d    %.30lf \n", j, melF);
            
            float fval = 0.0f;
            //float invMelStep = 1.0f / (float)melStep;
            
            if (melFr > melFreqs[i] && melFr < melFreqs[i+1]){
                //val = (melF - melFreqs[i]) / (melFreqs[i+1] - melFreqs[i]);
                fval = (melFr - melFreqs[i]) / melStep;
                //fval = (melFr - melFreqs[i]) / (melFreqs[i+1] - melFreqs[i]);
            } else if (melFr >= melFreqs[i+1]  && melFr < melFreqs[i+2] ) {
                //val = (melFreqs[i+2] - melF) / (melFreqs[i+2] - melFreqs[i+1]);
                fval = (melFreqs[i+2] - melFr) / melStep;
                //fval = (melFreqs[i+2] - melFr) / (melFreqs[i+2] - melFreqs[i+1]);
            }
                
            //melWindows[i][j] = val;
            melWindows[i*dftCoeffSize + j] = fval >= 0.0f  ? fval : 0.0f;
        }
    }
}
*/

//NS_ASSUME_NONNULL_END





/*
- (void)initWindow{  // hann
    sigWin = malloc(sizeof(float) * winSize);
    double ds = (double)winSize - 1.0;
    double A = M_PI / ds;
    for( int i = 0; i < winSize; i++){
        sigWin[i] = (float)pow(sin(A * i), 2.0);
    }
}
*/

/*
for i in range(1, int(dftSize/2)-1):
for k in range(dftSize):
   angle = 2 * np.pi * i * k / dftSize;
   dftmat[(i*2)-1,k] = np.cos(angle)
   dftmat[i*2,k] = np.sin(angle)
*/







/*
for( int i = 0; i < nframes; i++){
    if (windowedSig[i] != NULL) {free(windowedSig[i]);}
    if (windowedSignal != NULL) {free(windowedSignal[i]);}
    if (spectro[i] != NULL){ free(spectro[i]); }
    if (melSpectro[i] != NULL){ free(melSpectro[i]); }
}

if(spectro != NULL){ free(spectro);}
if(melSpectro != NULL){ free(melSpectro); }
if(windowedSig != NULL) { free(windowedSig); }
if(windowedSignal != NULL) { free(windowedSignal); }
*/

//for( int i=0; i < nMelWins; i++){
//    if(melWindows[i] != NULL) {free(melWindows[i]);}
//}
//if(melWindows != NULL){free(melWindows);}


//double melFr = melFs[j];
//double fval = 0.0;


/*
if (melFr > lmelFreqs[i] && melFr < cmelFreqs[i]){
    //fval = (melFr - melFreqs[i]) / melStep;
    fval = (melFr - lmelFreqs[i]) / (cmelFreqs[i] - lmelFreqs[i]);
} else if (melFr >= cmelFreqs[i]  && melFr < rmelFreqs[i] ) {
    //fval = (melFreqs[i+2] - melFr) / melStep;
    fval = (rmelFreqs[i] - melFr) / (rmelFreqs[i] - cmelFreqs[i]);
}
*/
