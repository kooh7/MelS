//
//  BinWriter.h
//  MelS
//
//  Created by Ken on 09/10/2025.
//

#import <Foundation/Foundation.h>

#ifndef binwriter_h
#define binwriter_h

@interface BinWriter : NSObject

- (BOOL)save2DCdblArray:(double **)array
              rows:(NSInteger)rows
             cols:(NSInteger)cols
                 toFile:(NSString *)filePath;


- (BOOL)save2DCArray:(float **)array
          rows:(NSInteger)rows
         cols:(NSInteger)cols
         toFile:(NSString *)filePath;

- (BOOL) save1DCArray:(float*)array
         entries:(NSInteger)entries
        toFile:(NSString *)filePath;

- (BOOL) save1DCdblArray:(double[])array
         rows:(NSInteger)rows
         cols:(NSInteger)cols
        toFile:(NSString *)filePath;
    
- (BOOL) save1DCdblArray:(double*)array
         entries:(NSInteger)entries
        toFile:(NSString *)filePath;

- (BOOL) save1DCArray:(float[])array
         rows:(NSInteger)rows
         cols:(NSInteger)cols
               toFile:(NSString *)filePath;

@end

#endif /* binwriter */
