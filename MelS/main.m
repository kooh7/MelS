//
//  main.m
//  MelS
//
//  Created by Ken on 08/10/2025.
//

#import <Foundation/Foundation.h>
#import <CoreML/CoreML.h>
#import "AudioLoader.h"
#import "MelSpecD.h"
#import "BinWriter.h"


MLMultiArray *createMLMultiArrayFromData(double **data, int numRows, int numCols);


int main(int argc, const char * argv[]) {
    @autoreleasepool {
        
        /*
        if (argc < 2) {
            NSLog(@"Usage: %@ <filename>", [NSString stringWithUTF8String:argv[0]]);
            return 1;
        }
        */
        
        // Get the filename from command line arguments  and parse the name of the output file from that
        //NSString *inputStr = [NSString stringWithUTF8String:argv[1]];
        NSString *inputStr = @"/Users/ken/Desktop/OmarWavs/split/1_chunk0.wav";
        NSUInteger slen = [inputStr length];
        NSString *baseStr = [inputStr substringToIndex:slen-4];
        NSString *binStr = @".bin";
        NSString *savePath = [baseStr stringByAppendingString:binStr];
        
    
        
        NSLog(@"Base : %@ <filename>", baseStr  ) ;
        NSLog(@"Bin  : %@ <filename>", savePath  );
        
        MelSpecD* mspecd = [[MelSpecD alloc] init];
        BinWriter* bwriter = [[BinWriter alloc] init];                           // wries to a binary file
        
        AudioLoader* aloader = [[AudioLoader alloc] initWithFilePath:inputStr];  // loads audio into a variable
        [aloader loadAudioFileAndConvertToFloatArray];
        float* sig = [aloader audioFloatArray];
        long siglen = [aloader numberOfSamples];
        
        mspecd = [[MelSpecD alloc] init];
        [mspecd setSignal:sig :siglen];                                         // starts the processing of mel spec
        //int nframes = [mspecd getNumFrames];
        int nframes = [mspecd getMaxLength];
        int melwins = [mspecd melwinsize];
        float* melSpectro = [mspecd getMelSpectro];                           // gets reference to the MelSpec
        [bwriter saveFlatTo2DCArray:(float*)melSpectro                           // writes the MelSpec to a binary file
                            rows:(NSInteger)nframes
                            cols:(NSInteger)melwins
                          toFile:(NSString *)savePath];
        
        /*
        [bwriter save2DCdblArray:(double**)melSpectro                           // writes the MelSpec to a binary file
                            rows:(NSInteger)nframes
                            cols:(NSInteger)melwins
                          toFile:(NSString *)savePath];
        */
        
        //[mspecd convertToMLMultiArray];
        //MLMultiArray* mlma = [mspecd getMLarray];
        
        
        // *mlma = convertToMLMA(melSpectro, nframes, melwins);       // converts the MelSpec C-array to a MLMultiArray for CoreML
    }
    return 0;
}
  








/*
 
 #import "melspec.h"
 NSArray* convertC2NSarray(float**, int, int);
int main(int argc, const char * argv[]) {
    @autoreleasepool {
        //Melspec* mspec = [[Melspec alloc] init];
        
        
        NSString* pathD = @"/Users/ken/Documents/sigwinD400.bin";
        int entries = [mspecd winsize];
        [bwriter save1DCdblArray:[mspecd sigwin] entries:entries toFile:pathD];
        
        NSString* pathMD = @"/Users/ken/Documents/sigwinMD400.bin";
        int rows = [mspecd coefsize];
        int cols = [mspecd melwinsize];
        [bwriter save1DCArray:[mspecd melwin] rows:rows cols:cols toFile:pathMD];
        //[bwriter save1DCdblArray:[mspecd melwin] rows:rows cols:cols toFile:pathMD];
        // melwindows[dftCoeffSize][mel_wins]
        
        NSString *directoryPath = @"/Users/ken/Desktop/OmarWavs/split/";
        NSMutableArray *fileList = [NSMutableArray array];
        
        // Simulate adding filenames
        [fileList addObject:@"1_chunk0"];
        [fileList addObject:@"1_chunk1"];
        [fileList addObject:@"1_chunk2"];
        [fileList addObject:@"1_chunk3"];
        [fileList addObject:@"1_chunk4"];
        [fileList addObject:@"1_chunk5"];
        [fileList addObject:@"1_chunk6"];
        [fileList addObject:@"1_chunk7_rem"];
        [fileList addObject:@"2_chunk0"];
        [fileList addObject:@"2_chunk1"];
        [fileList addObject:@"2_chunk2"];
        [fileList addObject:@"2_chunk3"];
        [fileList addObject:@"2_chunk4"];
        [fileList addObject:@"2_chunk5"];
        [fileList addObject:@"2_chunk6_rem"];
        
        NSString *wavext = @".wav";
        NSString *binext = @".bin";
        
        // Process full paths
        for (NSString *filename in fileList) {
            NSString *ofilename = [filename stringByAppendingString:wavext];
            NSString *sfilename = [filename stringByAppendingString:binext];
            NSString *fullPath = [directoryPath stringByAppendingPathComponent:ofilename ];
            NSString *savePath = [directoryPath stringByAppendingPathComponent:sfilename ];
            NSLog(@"Full path: %@", fullPath);
     
            [aloader setAudioFilePath: fullPath];
            [aloader loadAudioFileAndConvertToFloatArray];
            float* sig = [aloader audioFloatArray];
            long siglen = [aloader numberOfSamples];
            
            mspecd = [[MelSpecD alloc] init];
            [mspecd setSignal:sig :siglen];
            //int nframes = [mspecd getNumFrames];
            int nframes = 1024;
            int dftwinsize = [mspecd getdftsize];
            int melwins = [mspecd melwinsize];
            double** melSpectro = [mspecd getMelSpectro];
            [bwriter save2DCdblArray:(double**)melSpectro
                                rows:(NSInteger)nframes
                                cols:(NSInteger)melwins
                              toFile:(NSString *)savePath];
        }
 
        */
    
        /*
         
         
         NSString *inputStr = @"/Users/ken/Desktop/1_03.wav";
         //NSString *inputStr = @"/Users/ken/Documents/Khel.wav";
         AudioLoader* aloader = [[AudioLoader alloc] initWithFilePath:inputStr];
         [aloader loadAudioFileAndConvertToFloatArray];
         float* sig = [aloader audioFloatArray];
         NSString* pathSF = @"/Users/ken/Documents/sigF.bin";
         long siglen = [aloader numberOfSamples];
         [bwriter save1DCArray:sig entries: siglen toFile:pathSF];
         
        float** winput = [mspecd getWindowedSig];
       
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
        
     
        
        NSString* pathMel = @"/Users/ken/Desktop/1_03.bin";
        //NSString* pathMel = @"/Users/ken/Documents/Khel.bin";
        [bwriter save2DCdblArray:(float**)melSpectro
                              rows:(NSInteger)nframes
                             cols:(NSInteger)melwins
                        toFile:(NSString *)pathMel];

        
        
    }
    return 0;
}
         */




/*
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
         MelSpecD* mspecd;
         BinWriter* bwriter = [[BinWriter alloc] init];
         AudioLoader* aloader = [[AudioLoader alloc] init];
         
         NSString* pathD = @"/Users/ken/Documents/sigwinD400.bin";
         int entries = [mspecd winsize];
         [bwriter save1DCdblArray:[mspecd sigwin] entries:entries toFile:pathD];
         
         NSString* pathMD = @"/Users/ken/Documents/sigwinMD400.bin";
         int rows = [mspecd coefsize];
         int cols = [mspecd melwinsize];
         [bwriter save1DCArray:[mspecd melwin] rows:rows cols:cols toFile:pathMD];
         //[bwriter save1DCdblArray:[mspecd melwin] rows:rows cols:cols toFile:pathMD];
         // melwindows[dftCoeffSize][mel_wins]
         
         NSString *directoryPath = @"/Users/ken/Desktop/OmarWavs/split/";
         NSMutableArray *fileList = [NSMutableArray array];
         
         // Simulate adding filenames
         [fileList addObject:@"1_chunk0"];
         [fileList addObject:@"1_chunk1"];
         [fileList addObject:@"1_chunk2"];
         [fileList addObject:@"1_chunk3"];
         [fileList addObject:@"1_chunk4"];
         [fileList addObject:@"1_chunk5"];
         [fileList addObject:@"1_chunk6"];
         [fileList addObject:@"1_chunk7_rem"];
         [fileList addObject:@"2_chunk0"];
         [fileList addObject:@"2_chunk1"];
         [fileList addObject:@"2_chunk2"];
         [fileList addObject:@"2_chunk3"];
         [fileList addObject:@"2_chunk4"];
         [fileList addObject:@"2_chunk5"];
         [fileList addObject:@"2_chunk6_rem"];
         
         NSString *wavext = @".wav";
         NSString *binext = @".bin";
         
         // Process full paths
         for (NSString *filename in fileList) {
             NSString *ofilename = [filename stringByAppendingString:wavext];
             NSString *sfilename = [filename stringByAppendingString:binext];
             NSString *fullPath = [directoryPath stringByAppendingPathComponent:ofilename ];
             NSString *savePath = [directoryPath stringByAppendingPathComponent:sfilename ];
             NSLog(@"Full path: %@", fullPath);
      
             [aloader setAudioFilePath: fullPath];
             [aloader loadAudioFileAndConvertToFloatArray];
             float* sig = [aloader audioFloatArray];
             long siglen = [aloader numberOfSamples];
             
             mspecd = [[MelSpecD alloc] init];
             [mspecd setSignal:sig :siglen];
             //int nframes = [mspecd getNumFrames];
             int nframes = 1024;
             int dftwinsize = [mspecd getdftsize];
             int melwins = [mspecd melwinsize];
             double** melSpectro = [mspecd getMelSpectro];
             [bwriter save2DCdblArray:(double**)melSpectro
                                 rows:(NSInteger)nframes
                                 cols:(NSInteger)melwins
                               toFile:(NSString *)savePath];
         }
     
         ### Old start here
          
          
          NSString *inputStr = @"/Users/ken/Desktop/1_03.wav";
          //NSString *inputStr = @"/Users/ken/Documents/Khel.wav";
          AudioLoader* aloader = [[AudioLoader alloc] initWithFilePath:inputStr];
          [aloader loadAudioFileAndConvertToFloatArray];
          float* sig = [aloader audioFloatArray];
          NSString* pathSF = @"/Users/ken/Documents/sigF.bin";
          long siglen = [aloader numberOfSamples];
          [bwriter save1DCArray:sig entries: siglen toFile:pathSF];
          
         float** winput = [mspecd getWindowedSig];
        
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
         
      
         
         NSString* pathMel = @"/Users/ken/Desktop/1_03.bin";
         //NSString* pathMel = @"/Users/ken/Documents/Khel.bin";
         [bwriter save2DCdblArray:(float**)melSpectro
                               rows:(NSInteger)nframes
                              cols:(NSInteger)melwins
                         toFile:(NSString *)pathMel];
        
         
         
     }
     return 0;
 }
 */

