//
//  MelSpecD.m
//  MelS
//
//  Created by Ken on 10/10/2025.
//

#import "melspecD.h"
#import <math.h>
#import <stdlib.h>
#import <float.h>
//NS_ASSUME_NONNULL_BEGIN


@implementation MelSpecD

- (instancetype)init{
    [self initWindow];
    [self initDFT];
    [self initMelWindows];
    
    return self;
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

- (void)initWindow{  // hann
    //sigWin = malloc(sizeof(double) * winSize);
    double ds = (double)winSize - 1.0;
    //double A = M_PI / ds;
    for( double i = 0.0; i < winSize; i+=1.0){
        //double di = (double)i;
        sigWin[(int)i] = pow(sin( fmod(M_PI * i / ds, M_PI)), 2.0);
        //printf("%f ", i);
    }
}




- (void)initDFT{
    //dftmatrix = malloc(sizeof(float*) * dftSize);
    //for (int i = 0; i < dftSize; ++i){
    //    dftmatrix[i] = malloc(sizeof(float) * dftSize);
    //}
    for ( int k = 0; k < dftSize; k++){
        //dftmatrix[0][k] = 1.0f;
        dftmatrix[k] = 1.0f;
        if (k%2 == 0){
            //dftmatrix[dftSize-1][k] = 1.0;
            dftmatrix[((dftSize-1) *dftSize) + k] = 1.0;
        }
        else{
            //dftmatrix[dftSize-1][k] = -1.0f;
            dftmatrix[((dftSize-1) *dftSize) + k] = -1.0;
        }
    }
    
    double ds = (double)dftSize;
    
    for (int i = 1; i < (dftSize/2); i++){
        for ( int k = 0; k < dftSize; k++){
            double angle = (2 * i * k) * M_PI / ds;
            //dftmatrix[(i*2)-1][k] = cos(angle);
            //dftmatrix[(i*2)][k] = sin(angle);
            dftmatrix[(((i*2)-1) *dftSize)  + k] = cos(angle);
            dftmatrix[(((i*2)) *dftSize)  + k] = sin(angle);
        }
    }
}

- (void)windowSignal{
    
    nframes = 1 + (int)floor( (siglen - dftSize)/(hopSize) );
    
    winstart = (int*)malloc(sizeof(int) * nframes);
    windowedSignal = (float**)malloc(sizeof(float*) * nframes);
    windowedSig = (float**)malloc(sizeof(float*) * nframes);
    
    for(int i = 0; i < nframes; i++){
        winstart[i] = i * hopSize;
        windowedSignal[i] = (float*)calloc(dftSize, sizeof(float));
        windowedSig[i] = (float*)calloc(dftSize, sizeof(float));
    }
    
    float coltot, colmean;
    //int wstart;
    int olap = (dftSize - winSize) / 2;
    //printf("olap %d", olap);
    
    // Fill the WindowedSignal and take away the mean
    for(int i = 0; i < nframes; i++){
        
        int wsize;
        int instix = hopSize * i;
        
        coltot = 0.0f;
        //wstart = winstart[i];
        
        if (siglen - instix < winSize){
            wsize = siglen - instix;
            //psize = winSize - wsize;
        }else{
            wsize = winSize;
            //psize = 0;
        }
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
    }
    //row mean removal
    /*
    for (int j = olap; j <  winSize - olap; j++){
        double ctr = 0.0;
        for (int i = 0; i < nframes; i++){
            ctr += windowedSignal[i][j];
        }
        ctr /= nframes;
        for (int i = 0; i < nframes; i++){
            windowedSignal[i][j] -= ctr;
        }
    }*/
    
    // column based preemph
    /*
    for (int row = olap; row < dftSize-olap; row++) {
        windowedSig[0][row] = 0.03 * windowedSignal[0][row];
        for (int col = 1; col < nframes; col++) {
            windowedSig[col][row] = windowedSignal[col][row] - (preCoef * windowedSignal[col-1][row]);
        }
    }
     */
    
    for (int col = 0; col < nframes; col++) {
        windowedSig[col][olap] = 0.03 * windowedSignal[col][olap];
        for (int row = olap+1; row < winSize+olap; row++) {
            windowedSig[col][row] = windowedSignal[col][row] - (preCoef * windowedSignal[col][row-1]);
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
            windowedSig[j][i + olap] *= (float)sigWin[i];
        }
    }
}



- (void)calculateSpectro{
    
    // allocate memory for specs
    //spectro = malloc(sizeof(double*) * nframes);
    //melSpectro = malloc(sizeof(double*) * nframes);
    
    
    spectro = malloc(nframes * sizeof(double*));
    
    for( int i = 0; i < nframes; i++){
        spectro[i] = calloc(dftCoeffSize, sizeof(double) );
        
        double* specptr = spectro[i];
        float* sigptr = windowedSig[i];
        
        double oddcoeff = 0.0;
        double evencoeff = 0.0;
        int oddix, evenix;
             
       
             
        for (int j = 0; j < dftSize; j++){
            oddcoeff += (dftmatrix[j] * sigptr[j]);
            evencoeff += (dftmatrix[((dftSize-1) * dftSize) + j] * sigptr[j]);
        }
        if (oddcoeff < 0.0f){
            oddcoeff = -oddcoeff;
        }
        if (evencoeff < 0.0f){
            evencoeff = - evencoeff;
        }
        
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
    }
}

-(void)calculateMelSpec{
    melSpectro = malloc(sizeof(double*) * nframes);
    double stdsq = 4.5689974 * 2.0;
    double smean = -4.2677393;
    
    for (int i = 0; i < nframes; i++){      // sig cols
        
        melSpectro[i] = calloc(melWins, sizeof(double) );
        double totval = 0.0;
        for( int j = 0; j < melWins; j++){
            
            for (int k = 0; k < dftCoeffSize; k++){
                //melSpectro[j][i] += spectro[k][i] * melWindows[k][j];
                melSpectro[i][j] += (spectro[i][k] * (double)melWindows[(j*dftCoeffSize) + k]);
            }
            if (melSpectro[i][j] < (double)FLT_EPSILON){
                melSpectro[i][j] = log((double)(FLT_EPSILON));
            }
            else{
                melSpectro[i][j] = log(melSpectro[i][j]);
            }
            totval += melSpectro[i][j];
        }
        double totM = totval / melWins;
        for( int j = 0; j < melWins; j++){
            //melSpectro[i][j] -= totM;
            melSpectro[i][j] -= smean;
            melSpectro[i][j] /= stdsq;
        }
    }
}


/*
Eigen::MatrixXf build_mel_filterbank() {
    const int n_fft_bins = N_FFT / 2 + 1;
    Eigen::MatrixXf fb = Eigen::MatrixXf::Zero(N_MELS, n_fft_bins);

    const float mel_min = hz_to_mel_kaldi(F_MIN);
    const float mel_max = hz_to_mel_kaldi(F_MAX);
    const float mel_step = (mel_max - mel_min) / (N_MELS + 1);

    // FFT bin center freqs (Hz) → mel
    const float fft_bin_width = float(SAMPLE_RATE) / float(N_FFT);
    std::vector<float> mel_of_bin(n_fft_bins);
    for (int k = 0; k < n_fft_bins; ++k) {
        float hz = k * fft_bin_width;
        mel_of_bin[k] = hz_to_mel_kaldi(hz);
    }

    // For each mel filter, define left/center/right in mel,
    // then compute triangular weights against mel_of_bin
    for (int m = 0; m < N_MELS; ++m) {
        float left_m   = mel_min +  m      * mel_step;
        float center_m = mel_min + (m + 1) * mel_step;
        float right_m  = mel_min + (m + 2) * mel_step;

        for (int k = 0; k < n_fft_bins; ++k) {
            float mel_k = mel_of_bin[k];
            float w = 0.0f;

            if (mel_k > left_m && mel_k < center_m) {
                w = (mel_k - left_m) / std::max(center_m - left_m, 1e-12f);
            } else if (mel_k >= center_m && mel_k < right_m) {
                w = (right_m - mel_k) / std::max(right_m - center_m, 1e-12f);
            } else {
                w = 0.0f;
            }

            if (w < 0.0f) w = 0.0f;
            fb(m, k) = w;
        }
    }
    return fb;
}
*/

/*
// FFT bin center freqs (Hz) → mel
const float fft_bin_width = float(SAMPLE_RATE) / float(N_FFT);
std::vector<float> mel_of_bin(n_fft_bins);
for (int k = 0; k < n_fft_bins; ++k) {
    float hz = k * fft_bin_width;
    mel_of_bin[k] = hz_to_mel_kaldi(hz);
}
*/


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
    double melStep = (melMax - melMin) / (melWins + 1.0);
    double fftbinWidth = (fs / dftSize);
    
    //float melMin = [self HzToMelF:(float)minFreq];  //1127.0 * log(1.0 + minFreq / 700.0);
    //float melMax = [self HzToMelF:(float)maxFreq];  //1127.0 * log(1.0 + maxFreq / 700.0);
    //float melStep = (melMax - melMin) / (melWins + 1.0f);
    //float fftbinWidth = fs / dftSize;
    
    
    printf("MelMin %.30lf \n", melMin);
    printf("MelMax %.30lf \n", melMax);
    printf("MelStep %.30lf \n", melStep);
    printf("FFTBW %.30lf \n", fftbinWidth);
    
    //float melFreqs[melWins + 2];
    float cmelFreqs[melWins], lmelFreqs[melWins], rmelFreqs[melWins];
    float melFs[dftCoeffSize];
    //double melFreqs[melWins + 2];
    //double melFs[dftCoeffSize];
    
    /*
    for (int i = 0; i < melWins + 2; i++) {
        //melFreqs[i] = (float)(melMin) + (i * (float)(melStep));        // 6.3 -9    000996     // melscale in increments
        //melFreqs[i] = (float)(melMin) + (float)(i * melStep);          // 8.5 -9      00128
        melFreqs[i] = (float)(melMin + (i * melStep));                 // 3.3 -9      00054
        //melFreqs[i] = (melMin + (i * melStep));
        printf("mel %.20f   : ix  %d \n", melFreqs[i], i  );
    }
    
    
    for (int i = 0; i < melWins; i++) {
        lmelFreqs[i] = (float)(melMin + (i * melStep));                 // 3.3 -9      00054
        cmelFreqs[i] = (float)(melMin + ((i + 1.0) * melStep));
        rmelFreqs[i] = (float)(melMin + ((i + 2.0) * melStep));
        //printf("mel %.20f   : ix  %d \n", melFreqs[i], i  );
    }
     */
    
    for (int i = 0; i < melWins; i++) {
        lmelFreqs[i] = (melMin + (i * melStep));                 // 3.3 -9      00054
        cmelFreqs[i] = (melMin + ((i + 1.0f) * melStep));
        rmelFreqs[i] = (melMin + ((i + 2.0f) * melStep));
        //printf("mel %.20f   : ix  %d \n", melFreqs[i], i  );
    }
    
    for(int j = 0; j < dftCoeffSize-1; j++){
        float fr = (float)(j * fftbinWidth);
        melFs[j] = (float)[self HzToMelF:(float)fr];
        
        //double fr = j * fftbinWidth;
        //melFs[j] = (float)[self HzToMel:fr];
        //printf("MelF %d    %.30lf \n", j, melFs[j]);
    }
    melFs[dftCoeffSize-1] = 0.0;
    
    for( int i=0; i < melWins; i++){
        for(int j = 0; j < dftCoeffSize; j++){
            float melFr = melFs[j];
            float fval = 0.0;
            //float mstep = (float)melStep;
            
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
            
            if (melFr > lmelFreqs[i] && melFr < cmelFreqs[i]){
                //fval = (melFr - melFreqs[i]) / melStep;
                fval = (melFr - lmelFreqs[i]) / (cmelFreqs[i] - lmelFreqs[i]);
            } else if (melFr >= cmelFreqs[i]  && melFr < rmelFreqs[i] ) {
                //fval = (melFreqs[i+2] - melFr) / melStep;
                fval = (rmelFreqs[i] - melFr) / (rmelFreqs[i] - cmelFreqs[i]);
            }
            /*
            if (melFr > melFreqs[i] && melFr < melFreqs[i+1]){
                //fval = (melFr - melFreqs[i]) / melStep;
                fval = (melFr - melFreqs[i]) / (melFreqs[i+1] - melFreqs[i]);
            } else if (melFr >= melFreqs[i+1]  && melFr < melFreqs[i+2] ) {
                //fval = (melFreqs[i+2] - melFr) / melStep;
                fval = (melFreqs[i+2] - melFr) / (melFreqs[i+2] - melFreqs[i+1]);
            }
            */
            melWindows[i*dftCoeffSize + j] = fval >= 0.0f  ? (float)fval : 0.0f;
        }
    }
}
    
    


/*

- (void)initMelWindows{
    double melMin = 1127.0 * log(1.0 + minFreq / 700.0);
    double melMax = 1127.0 * log(1.0 + maxFreq / 700.0);
    double melStep = (melMax - melMin) / (melWins + 1.0);
    
    float melFreqs[melWins + 2];  //, melHz[melWins + 2];
    float melFs[dftCoeffSize];
    
    double fftbinWidth = fs / dftSize;
    double freq, val, melF;
    
    printf("MelMin %.30lf \n", melMin);
    printf("MelMax %.30lf \n", melMax);
    printf("MelStep %.30lf \n", melStep);
    printf("FFTBW %.30lf \n", fftbinWidth);
    
    for (int i = 0; i < melWins + 2; i++) {
        melFreqs[i] = melMin + (float)(i * melStep);         // melscale in increments of melstep
        //melHz[i] = [self melToHz: melFreqs[i]];     // the melscale in Hz ( with diff index )
        printf("mel %.20f   : ix  %d \n", melFreqs[i], i  );
    }
    
    for(int j = 0; j < dftCoeffSize; j++){
        float fr = j * (float)fftbinWidth;
        melFs[j] = [self HzToMelF:fr];
        printf("MelF %d    %.30lf \n", j, melFs[j]);
    }

    //melWindows = malloc(sizeof(double*) * melWins);
    for( int i=0; i < melWins; i++){
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

-(void) dealloc{
    //free(sigWin);
    free(winstart);
    winstart = NULL;
    //free(winend);
   
    //for (int i = 0; i < dftSize; ++i){
    //    free(dftmatrix[i]);
    //}
    //free(dftmatrix);
    
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
    
    //for( int i=0; i < melWins; i++){
    //    if(melWindows[i] != NULL) {free(melWindows[i]);}
    //}
    //if(melWindows != NULL){free(melWindows);}
    

}


-(double*)sigwin{
    return sigWin;
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
    return melWins;
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

- (double**)getMelSpectro{
    return melSpectro;
}

-(double*) getdftMat{
    return &dftmatrix[0];
}


@end

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
