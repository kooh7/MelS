//
//  main.m
//  MelS
//
//  Created by Ken on 08/10/2025.
//

#import <Foundation/Foundation.h>
#import "audioloader.h"
#import "melspec.h"
#import "melspecD.h"
#import "binwriter.h"
//#import <math.h>


NSArray* convertC2NSarray(float*, int, int);

int main(int argc, const char * argv[]) {
    @autoreleasepool {
        // insert code here...
        //NSLog(@"Hello, World!");
        //NSString *inputStr = [NSString stringWithUTF8String:argv[1]];
        
        //printf("%f ", M_PI);
        
        //NString *inputStr = '/Users/ken/Documents/'
        //AudioLoader* aloader = [[AudioLoader alloc] initWithFilePath:inputStr];
        //[aloader loadAudioFileAndConvertToFloatArray];
        //float* sig = [aloader audioFloatArray];
        //long siglen = [aloader numberOfSamples];
        
        //Melspec* mspec = [[Melspec alloc] init];
        MelSpecD* mspecd = [[MelSpecD alloc] init];

        BinWriter* bwriter = [[BinWriter alloc] init];
    
        NSString* pathD = @"/Users/ken/Documents/sigwinD400.bin";
        int entries = [mspecd winsize];
        [bwriter save1DCdblArray:[mspecd sigwin] entries:entries toFile:pathD];
        
        NSString* pathMD = @"/Users/ken/Documents/sigwinMD400.bin";
        int rows = [mspecd coefsize];
        int cols = [mspecd melwinsize];
        [bwriter save1DCArray:[mspecd melwin] rows:rows cols:cols toFile:pathMD];
        //[bwriter save1DCdblArray:[mspecd melwin] rows:rows cols:cols toFile:pathMD];
        // melwindows[dftCoeffSize][mel_wins]
        
        
        
        NSString *inputStr = @"/Users/ken/Desktop/2_06.wav";
        AudioLoader* aloader = [[AudioLoader alloc] initWithFilePath:inputStr];
        [aloader loadAudioFileAndConvertToFloatArray];
        float* sig = [aloader audioFloatArray];
        NSString* pathSF = @"/Users/ken/Documents/sigF.bin";
        long siglen = [aloader numberOfSamples];
        [bwriter save1DCArray:sig entries: siglen toFile:pathSF];
        
        
        [mspecd setSignal:sig :siglen];
        float** winput = [mspecd getWindowedSig];
        int nframes = [mspecd getNumFrames];
        int dftwinsize = [mspecd getdftsize];
        printf("WinSig  %d   %d ", nframes, dftwinsize);
        
        NSString* pathWS = @"/Users/ken/Documents/winSig.bin";
        [bwriter save2DCArray:(float**)winput
                              rows:(NSInteger)nframes
                             cols:(NSInteger)dftwinsize
                        toFile:(NSString *)pathWS];
        
        
        double* dftmat = [mspecd getdftMat];
        NSString* dftpath = @"/Users/ken/Documents/dft.bin";
        [bwriter save1DCdblArray:(double*)dftmat
                              rows:(NSInteger)dftwinsize
                             cols:(NSInteger)dftwinsize
                        toFile:(NSString *)dftpath];
        

        
        double** spectro = [mspecd getSpectro];
        double** melSpectro = [mspecd getMelSpectro];
        //int nframes = [mspecd getNumFrames];
        //int dftwinsize = [mspecd getdftsize];
        int melwins = [mspecd melwinsize];
        
        NSString* pathSpec = @"/Users/ken/Documents/Spec.bin";
        [bwriter save2DCdblArray:(double**)spectro
                              rows:(NSInteger)nframes
                             cols:(NSInteger)dftCoeffSize
                        toFile:(NSString *)pathSpec];
     
        
        NSString* pathMel = @"/Users/ken/Desktop/2_06.bin";
        [bwriter save2DCdblArray:(float**)melSpectro
                              rows:(NSInteger)nframes
                             cols:(NSInteger)melwins
                        toFile:(NSString *)pathMel];
        
        
    }
    return 0;
}

/*
NSArray* convertC2NSarray(float*, int rows, int cols){
    NSMutableArray array =   ;
    for (int i = 0; i < rows; i++){
        NSArray rowarr =
        array[i] = arr;
    }
}
*/


// call file loader with the file path
// call the melspec class with the pointer to the raw data and the size of the raw data
// call the save data file with a pointer to the 
// dealloc / exit
