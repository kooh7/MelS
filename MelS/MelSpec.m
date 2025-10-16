//
//  melspec.m
//  MelS
//
//  Created by Ken on 08/10/2025.
//

//#import <Foundation/Foundation.h>
#import<stdlib.h>
#import <math.h>
#import "melspec.h"


@implementation Melspec




- (instancetype)init{
    [self initWindow];
    [self initDFT];
    [self initMelWindows];
    return self;
}


-(void)setSignal:(float*)isig :(int)sigl{
    sig = isig;
    siglen = sigl;
    [self windowSignal];                       // done
    [self getSpectro];                         // add DFT function
    [self calculateMelSpec];
}


/*
- (void)initWindow{  // hann
    sigWin = malloc(sizeof(double) * winSize);
    double ds = (double)winSize - 1.0;
    double A = M_PI / ds;
    for( double i = 0.0; i < winSize; i+=1.0){
        //double di = (double)i;
        sigWin[(int)i] = pow(sin( fmod(M_PI * i / ds, M_PI)), 2.0);
        //printf("%f ", i);
    }
}
*/


- (void)initWindow{  // hann
    sigWin = malloc(sizeof(float) * winSize);
    double ds = (double)winSize - 1.0;
    //double A = M_PI / ds;
    for( int i = 0; i < winSize; i++){
        //sigWin[i] = (float)pow(sin(A * i), 2.0);
        sigWin[i] = pow(sin( fmod(M_PI * i / ds, M_PI)), 2.0);
    }
}


/*
for i in range(1, int(dftSize/2)-1):
for k in range(dftSize):
   angle = 2 * np.pi * i * k / dftSize;
   dftmat[(i*2)-1,k] = np.cos(angle)
   dftmat[i*2,k] = np.sin(angle)
*/


- (void)initDFT{
    dftmatrix = malloc(sizeof(float*) * dftSize);
    for (int i = 0; i < dftSize; ++i){
        dftmatrix[i] = malloc(sizeof(float) * dftSize);
     }
    for ( int k = 0; k <  dftSize; k++){
        dftmatrix[0][k] = 1.0f;
        if (k%2 == 0){
            dftmatrix[dftSize-1][k] = 1.0f;
        }
        else{
            dftmatrix[dftSize-1][k] = -1.0f;
        }
    }
    
    double ds = (double)dftSize;
    
    for (int i = 1; i < (dftSize/2) - 1; i++){
        for ( int k = 0; k < dftSize; k++){
            double angle = (double)(2 * i * k) * M_PI / ds;
            dftmatrix[(i*2)-1][k] = (float)cos(angle);
            dftmatrix[(i*2)][k] = (float)sin(angle);
        }
    }
}

- (void)windowSignal{
    nwindows = siglen / hopSize;   // TODO: get the correct value here
    winstart = malloc(sizeof(int) * nwindows);
    winend = malloc(sizeof(int) * nwindows);
    for (int i = 0; i < nwindows; ++i){
        winstart[i] = i * hopSize;
        winend[i] = (i * hopSize) + dftSize;
    }
    
    spectro = malloc(sizeof(float*) * nwindows);
    melSpectro = malloc(sizeof(float*) * nwindows);
    for( int i = 0; i < nwindows; i++){
        spectro[i] = malloc(sizeof(float) * dftCoeffSize);
        melSpectro[i] = malloc(sizeof(float) * melWins);
    }
}



- (void)getSpectro{
    float winSig[dftSize];
    for(int i = 0; i < nwindows; i++){
        int js = i * hopSize;
        for(int j = 0; j < dftSize; j++){
            winSig[j] = sig[js + j] * sigWin[j];
        }
        
        // TODO do the dft here to WinSIg
        
    }
}

-(void)calculateMelSpec{
    for (int i = 0; i < nwindows; i++){
        for( int j = 0; j < melWins; j++){
            for (int k = 0; k < dftCoeffSize; k++){
                melSpectro[j][i] += spectro[k][i] * melWindows[k][j];
            }
        }
    }
}


- (float)melToHz:(float)melFreq {
    return 700 * (exp(melFreq / 1127.0) - 1);
}


- (void)initMelWindows{
    float melMin = 1127.0f * logf(1.0f + minFreq / 700.0f);
    float melMax = 1127.0f * logf(1.0f + maxFreq / 700.0f);
    float melStep = (melMax - melMin) / (melWins + 1);
    float melFreqs[melWins + 2], melHz[melWins + 2];
        
    float freq, val, melF;
    
    for (int i = 0; i < melWins + 2; i++) {
        melFreqs[i] = melMin + i * melStep;
        melHz[i] = [self melToHz: melFreqs[i]];
    }
    

    melWindows = malloc(sizeof(float*) * melWins);
    for( int i=0; i < melWins; i++){
        melWindows[i] = malloc(sizeof(float) * dftCoeffSize);
        
        for(int j = 0; j < dftCoeffSize; j++){
            freq = (float)(j * fs/2.0f)  / ((float)dftCoeffSize-1.0f);   // or fs/dftCOeffSize   # is it -1
            melF = 1127.0f * logf(1.0f + freq / 700.0f);                // melF is the mel at each point oo
            val = 0.0f;
            if (melF > melFreqs[i] & melF < melFreqs[i+1]){
                val = (melF - melFreqs[i]) / (melFreqs[i+1] - melFreqs[i]);
            } else if (melF >= melFreqs[i+1]  && melF < melFreqs[i+2] ) {
                val = (melFreqs[i+2] - melF) / (melFreqs[i+2] - melFreqs[i+1]);
            }
                
            /*
            if( freq > melHz[i] && freq < melHz[i+1] ){
                val = (freq - melHz[i]) / (melHz[i+1] - melHz[i]);
            } else if (freq >= melHz[i+1]  && freq < melHz[i+2] ) {
                val = (melHz[i+2] - freq) / (melHz[i+2] - melHz[i+1]);
            }
            */
            melWindows[i][j] = val >= 0.0  ? val : 0.0f;
            //melWindows[i*dftCoeffSize + j] = val >= 0.0  ? val : 0.0f;
        }
    }
}


/*
- (void)initMelWindows{
    float melMin = 1127.0f * log(1.0f + minFreq / 700.0f);
    float melMax = 1127.0f * log(10.f + maxFreq / 700.0f);
    float melStep = (melMax - melMin) / (melWins + 1);
    float melFreqs[melWins + 2], melHz[melWins + 2];

    
    float freq, val;
    
    for (int i = 0; i < melWins + 2; i++) {
        melFreqs[i] = melMin + i * melStep;
        melHz[i] = [self melToHz: melFreqs[i]];
    }

    melWindows = malloc(sizeof(float*) * melWins);
    for( int i=0; i < melWins; i++){
        melWindows[i] = malloc(sizeof(float) * dftCoeffSize);
        for(int j = 0; j < dftCoeffSize; j++){
            freq = (float)(j * fs)  / dftCoeffSize;
            val = 0.0f;
            if( freq > melHz[i] && freq < melHz[i+1] ){
                val = (freq - melHz[i]) / (melHz[i+1] - melHz[i]);
            } else if (freq >= melHz[i+1]  && freq < melHz[i+2] ) {
                val = (melHz[i+2] - freq) / (melHz[i+2] - melHz[i+1]);
            }
            melWindows[i][j] = val;
        }
    }
}
*/

-(void) dealloc{
    free(sigWin);
    free(winstart);
    free(winend);
   
    for (int i = 0; i < dftSize; ++i){
        free(dftmatrix[i]);
    }
    free(dftmatrix);
    
    for( int i = 0; i < nwindows; i++){
        free(spectro[i]);
        free(melSpectro[i]);
    }
    
    free(spectro);
    free(melSpectro);
    
    for( int i=0; i < melWins; i++){
        free(melWindows[i]);
    }
    free(melWindows);
}


-(float*)sigwin{
    return sigWin;
}

-(int)winsize{
    return winSize;
}

- (float**)melwin{
    return melWindows;
}

- (int)melwinsize{
    return melWins;
}

- (int)coefsize{
    return dftCoeffSize;
}

@end

