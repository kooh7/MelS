//
//  MelSpecD.h
//  MelS
//
//  Created by Ken on 10/10/2025.
//

#import <Foundation/Foundation.h>
#import <CoreML/CoreML.h>

#ifndef melspecD_h
#define melspecD_h

#define dftSize 512
#define winSize 400
#define hopSize 160
#define fs 16000.0
#define dftCoeffSize 257
#define nMelWins 128
#define minFreq 20.0
#define maxFreq 8000.0
#define preCoef 0.97
#define maxLength 1024



@interface MelSpecD : NSObject{

    //double* sig;
    float* sig;
    long siglen;

    int* winstart;
    float** windowedSig;
    float** windowedSignal;

    double sigWin[winSize];
    double dftmatrix[dftSize*dftSize];
    float melWindows[nMelWins*dftSize];
    
    int nframes;
    double** spectro;
    float melSpectro[nMelWins*maxLength];
    
    MLMultiArray* mlArray;
}


- (void)initWindow;
- (void)initDFT;
- (void)initMelWindows;

- (void)setSignal:(float*)isig :(int)sigl;

- (void)windowSignal;
- (void)calculateSpectro;
- (void)calculateMelSpec;
- (void)convertToMLMultiArray;

- (double)melToHz:(double)melFreq;
- (double)HzToMel:(double)hzf;
- (float)HzToMelF:(float)hzf;


//the functions are for testing - allow values to be called and saved from main.m
- (int)winsize;
- (int)getdftsize;
- (double*)sigwin;
- (float*)melwin;

- (int)melwinsize;
- (int)coefsize;
- (int)getNumFrames;
- (int)getMaxLength;

- (double**)getSpectro;
//- (double**)getMelSpectro;
//- (double*)getMelSpectro;
- (float*)getMelSpectro;
- (float**)getWindowedSig;
- (double*)getdftMat;
- (MLMultiArray*)getMLarray;

@end


#endif /* MelSpecD_h */
