//
//  MelSpecD.h
//  MelS
//
//  Created by Ken on 10/10/2025.
//

#import <Foundation/Foundation.h>

#ifndef melspecD_h
#define melspecD_h

#define dftSize 512
#define winSize 400
#define hopSize 160
#define fs 16000.0
#define dftCoeffSize 257
#define melWins 128
#define minFreq 20.0
#define maxFreq 8000.0
#define preCoef 0.97



@interface MelSpecD : NSObject{

    //double* sig;
    float* sig;
    long siglen;

    int* winstart;
    float** windowedSig;
    float** windowedSignal;

    
    double sigWin[winSize];
    //double psigWin[dftSize];
    double dftmatrix[dftSize*dftSize];
    float melWindows[melWins*dftSize];
    

    //double** melStarts;
    //int** melLens;

    int nframes;
    double** spectro;
    double** melSpectro;
}



// Method to process the command-line arguments
//- (void)processArguments:(int)argc argv:(const char **)argv;

//- (void)setFileName:(NSString*)fileName;    // done
//- (void)loadFromFile;                       // done
//- (void)saveToFile;

- (void)initWindow;                         // done
- (void)initDFT;                            // done
- (void)initMelWindows;                     // done // add normalisation

- (void)setSignal:(float*)isig :(int)sigl;
- (void)windowSignal;                       // done


- (void)calculateSpectro;                         // add DFT function

- (void)calculateMelSpec;                   // done


- (double)melToHz:(double)melFreq;
- (double)HzToMel:(double)hzf;
- (float)HzToMelF:(float)hzf;
//- (void)dealloc;

//testing functions
- (int)winsize;
- (int)getdftsize;
- (double*)sigwin;
- (float*)melwin;

- (int)melwinsize;
- (int)coefsize;
- (int)getNumFrames;

- (double**)getSpectro;
- (double**)getMelSpectro;
- (float**)getWindowedSig;
- (double*)getdftMat;

@end




#endif /* MelSpecD_h */
