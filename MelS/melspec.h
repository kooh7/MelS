//
//  melspec.h
//  MelS
//
//  Created by Ken on 08/10/2025.
//

#import <Foundation/Foundation.h>

#ifndef melspec_h
#define melspec_h

#define dftSize 512
#define winSize 400
#define hopSize 160
#define fs 16000.0
#define dftCoeffSize 257
#define melWins 128
#define minFreq 20.0
#define maxFreq 8000.0



@interface Melspec : NSObject{

    //NSString* filename;
    float* sig;
    long siglen;

    float* winstart;
    float* winend;

    //float sigWin[dftSize];
    
    float* sigWin;
    float** dftmatrix;
    //const int melWins; //= 128;
    float** melWindows;
    
    //float** melStarts;
    //int** melLens;

    float** spectro;

    int nwindows;
    float** melSpectro;
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
- (void)getSpectro;                         // add DFT function
- (void)calculateMelSpec;                   // done


- (float)melToHz:(float)melFreq;
//- (void)dealloc;

//testing functions
- (int)winsize;
- (float*)sigwin;
- (float**)melwin;
- (int)melwinsize;
- (int)coefsize;

@end





#endif /* melspec_h */
