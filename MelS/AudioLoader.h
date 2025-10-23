//
//  audioloader.h
//  MelS
//
//  Created by Ken on 08/10/2025.
//

#import <Foundation/Foundation.h>
#import <AVFoundation/AVFoundation.h>


#ifndef audioloader_h
#define audioloader_h


@interface AudioLoader : NSObject

@property (nonatomic, strong) AVAudioFile *audioFile;
@property (nonatomic, strong) NSData *audioData;  // Raw audio data
@property (nonatomic, assign) float *audioFloatArray; // Pointer to the float array of audio data
@property (nonatomic, assign) NSInteger numberOfSamples;  // Number of samples read
@property (nonatomic, assign) NSInteger numberOfChannels; // Number of channels (stereo = 2)

- (instancetype)initWithFilePath:(NSString *)path;
- (void)setAudioFilePath:(NSString *)path;
- (BOOL)loadAudioFileAndConvertToFloatArray;

@end


#endif /* audioloader_h */
