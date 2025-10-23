//
//  BinWriter.m
//  MelS
//
//  Created by Ken on 09/10/2025.
//

#import "BinWriter.h"

@implementation BinWriter


// Save a C 2D array of floats to a binary file
- (BOOL)save2DCdblArray:(double **)array
              rows:(NSInteger)rows
             cols:(NSInteger)cols
            toFile:(NSString *)filePath {
    
    // Create a mutable data object to hold the raw binary data
    NSMutableData *data = [NSMutableData data];
    
    // Write the dimensions of the array (rows and cols) to the file first
    [data appendBytes:&rows length:sizeof(NSInteger)];
    [data appendBytes:&cols length:sizeof(NSInteger)];
    
    // Write the 2D array data (floats) to the binary file
    for (NSInteger i = 0; i < rows; i++) {
        [data appendBytes:array[i] length:cols * sizeof(double)];
    }
    
    // Write the data to the file
    NSError *error = nil;
    BOOL success = [data writeToFile:filePath atomically:YES];
    
    if (!success) {
        NSLog(@"Failed to save file: %@", error.localizedDescription);
    }
    
    return success;
}


// Save a C 2D array of floats to a binary file
- (BOOL)save2DCArray:(float **)array
              rows:(NSInteger)rows
             cols:(NSInteger)cols
            toFile:(NSString *)filePath {
    
    // Create a mutable data object to hold the raw binary data
    NSMutableData *data = [NSMutableData data];
    
    // Write the dimensions of the array (rows and cols) to the file first
    [data appendBytes:&rows length:sizeof(NSInteger)];
    [data appendBytes:&cols length:sizeof(NSInteger)];
    
    // Write the 2D array data (floats) to the binary file
    for (NSInteger i = 0; i < rows; i++) {
        [data appendBytes:array[i] length:cols * sizeof(float)];
    }
    
    // Write the data to the file
    NSError *error = nil;
    BOOL success = [data writeToFile:filePath atomically:YES];
    
    if (!success) {
        NSLog(@"Failed to save file: %@", error.localizedDescription);
    }
    
    return success;
}

-(BOOL) saveFlatTo2DCdblArray:(double*)melSpectro                           // writes the MelSpec to a binary file
                    rows:(NSInteger)nframes
                    cols:(NSInteger)melwins
                  toFile:(NSString *)filePath{
    
    NSMutableData* data =[NSMutableData data];
    
    [data appendBytes:&nframes length:sizeof(NSInteger)];
    [data appendBytes:&melwins length:sizeof(NSInteger)];
    
    // Write the 2D array data (floats) to the binary file
    for (NSInteger i = 0; i < nframes; i++) {
        [data appendBytes:&melSpectro[i*melwins] length:melwins * sizeof(double)];
    }
    
    // Write the data to the file
    NSError *error = nil;
    BOOL success = [data writeToFile:filePath atomically:YES];
    
    if (!success) {
        NSLog(@"Failed to save file: %@", error.localizedDescription);
    }
    
    return success;
}
    
-(BOOL) saveFlatTo2DCArray:(float*)melSpectro                           // writes the MelSpec to a binary file
                    rows:(NSInteger)nframes
                    cols:(NSInteger)melwins
                  toFile:(NSString *)filePath{
    
    NSMutableData* data =[NSMutableData data];
    
    [data appendBytes:&nframes length:sizeof(NSInteger)];
    [data appendBytes:&melwins length:sizeof(NSInteger)];
    
    // Write the 2D array data (floats) to the binary file
    for (NSInteger i = 0; i < nframes; i++) {
        [data appendBytes:&melSpectro[i*melwins] length:melwins * sizeof(float)];
    }
    
    // Write the data to the file
    NSError *error = nil;
    BOOL success = [data writeToFile:filePath atomically:YES];
    
    if (!success) {
        NSLog(@"Failed to save file: %@", error.localizedDescription);
    }
    
    return success;
}





- (BOOL) save1DCArray:(float*)array
         entries:(NSInteger)entries
              toFile:(NSString *)filePath{
    
    NSMutableData *data = [NSMutableData data];
    
    // Write the dimensions of the array (rows and cols) to the file first
    //[data appendBytes:&entries length:sizeof(NSInteger)];
    [data appendBytes : array length:entries * sizeof(float)];
    
    NSError *error = nil;
    BOOL success = [data writeToFile:filePath atomically:YES];
    
    if (!success) {
        NSLog(@"Failed to save file: %@", error.localizedDescription);
    }
    
    return success;
}

- (BOOL) save1DCdblArray:(double*)array
         entries:(NSInteger)entries
              toFile:(NSString *)filePath{
    
    NSMutableData *data = [NSMutableData data];
    
    // Write the dimensions of the array (rows and cols) to the file first
    //[data appendBytes:&entries length:sizeof(NSInteger)];
    [data appendBytes : array length:entries * sizeof(double)];
    
    
    
    NSError *error = nil;
    BOOL success = [data writeToFile:filePath atomically:YES];
    
    if (!success) {
        NSLog(@"Failed to save file: %@", error.localizedDescription);
    }
    
    return success;
}

- (BOOL) save1DCdblArray:(double[])array
         rows:(NSInteger)rows
         cols:(NSInteger)cols
        toFile:(NSString *)filePath{
    
    NSMutableData *data = [NSMutableData data];
    
    // Write the dimensions of the array (rows and cols) to the file first
    [data appendBytes:&rows length:sizeof(NSInteger)];
    [data appendBytes:&cols length:sizeof(NSInteger)];
    
    [data appendBytes : array length:rows*cols * sizeof(double)];
    
    NSError *error = nil;
    BOOL success = [data writeToFile:filePath atomically:YES];
    if (!success) {
        NSLog(@"Failed to save file: %@", error.localizedDescription);
    }
    return success;
}

- (BOOL) save1DCArray:(float[])array
         rows:(NSInteger)rows
         cols:(NSInteger)cols
        toFile:(NSString *)filePath{
    
    NSMutableData *data = [NSMutableData data];
    
    // Write the dimensions of the array (rows and cols) to the file first
    [data appendBytes:&rows length:sizeof(NSInteger)];
    [data appendBytes:&cols length:sizeof(NSInteger)];
    
    [data appendBytes : array length:rows*cols * sizeof(float)];
    
    NSError *error = nil;
    BOOL success = [data writeToFile:filePath atomically:YES];
    if (!success) {
        NSLog(@"Failed to save file: %@", error.localizedDescription);
    }
    return success;
}

        
@end


