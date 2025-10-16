//
//  AudioLoader.m
//  MelS
//
//  Created by Ken on 08/10/2025.
//

#import "AudioLoader.h"

//NS_ASSUME_NONNULL_BEGIN

@implementation AudioLoader


- (instancetype)initWithFilePath:(NSString *)path {
    self = [super init];
    if (self) {
        // Convert the file path to NSURL
        NSURL *url = [NSURL fileURLWithPath:path];
        NSError *error = nil;

        // Initialize AVAudioFile for reading the audio file
        self.audioFile = [[AVAudioFile alloc] initForReading:url error:&error];
        if (error) {
            NSLog(@"Error loading audio file: %@", error.localizedDescription);
            return nil;
        }
    }
    return self;
}

- (BOOL)loadAudioFileAndConvertToFloatArray {
    NSError *error = nil;

    // Read audio data into an AVAudioPCMBuffer
    AVAudioFormat *audioFormat = self.audioFile.processingFormat;
    AVAudioFrameCount frameCount = (AVAudioFrameCount)self.audioFile.length;
    AVAudioPCMBuffer *pcmBuffer = [[AVAudioPCMBuffer alloc] initWithPCMFormat:audioFormat frameCapacity:frameCount];

    // Load the audio file into the PCM buffer
    [self.audioFile readIntoBuffer:pcmBuffer error:&error];
    if (error) {
        NSLog(@"Error reading audio file into buffer: %@", error.localizedDescription);
        return NO;
    }

    // Get the number of channels and samples
    self.numberOfChannels = pcmBuffer.format.channelCount;
    self.numberOfSamples = (NSInteger)pcmBuffer.frameLength;

    // Allocate memory for the float array (samples * channels)
    self.audioFloatArray = (float *)malloc(sizeof(float) * self.numberOfSamples * self.numberOfChannels);
    if (self.audioFloatArray == NULL) {
        NSLog(@"Error allocating memory for float array.");
        return NO;
    }

    // Convert the audio data into a float array
    for (AVAudioFrameCount frame = 0; frame < pcmBuffer.frameLength; frame++) {
        for (AVAudioChannelCount channel = 0; channel < self.numberOfChannels; channel++) {
            // Get the sample value for the current frame/channel, scale to float (-1.0 to 1.0 range)
            float sample = *((float *)pcmBuffer.floatChannelData[channel] + frame);
            self.audioFloatArray[frame * self.numberOfChannels + channel] = sample;
        }
    }

    NSLog(@"Successfully loaded audio file and converted to float array.");
    return YES;
}

- (void)dealloc {
    if (self.audioFloatArray) {
        free(self.audioFloatArray);  // Don't forget to free the allocated memory
    }
}

@end



//NS_ASSUME_NONNULL_END
