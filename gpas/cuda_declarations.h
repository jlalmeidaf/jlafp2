#ifndef _CUDA_DECLARATIONS_H_
#define _CUDA_DECLARATIONS_H_


/*******************************************************************************
 * __CUDACC__ is defined only if we include this file into *.cu file           *
 * and compile it with nvcc                                                    *
 *                                                                             *
 * DO NOT INCLUDE THIS FILE TO *.CPP                                           *
 *******************************************************************************/

    #ifndef __CUDACC__

        //this will be compiled for example in cpp files (although it shouldn't)

        #include "texture_types.h"

        // define the keywords, so that the IDE does not complain about them
        #define __global__
        #define __device__
        #define __shared__
        #define __constant__

        template <typename T, int I = 1, enum cudaTextureReadMode E = cudaReadModeElementType> struct texture
        {
            int                          normalized;
            enum cudaTextureFilterMode   filterMode;
            enum cudaTextureAddressMode  addressMode[3];
            struct cudaChannelFormatDesc channelDesc;
        };

        char tex1Dfetch(texture<char, 1, cudaReadModeElementType> t, int x);

        int3 threadIdx, blockIdx, blockDim, gridDim;

        __device__ float signbit(float);
        
    #endif




    //INLINES
    inline void checkError(const char *errorMessage)
    {
        cudaError_t err = cudaGetLastError();
        if(cudaSuccess != err)
        {
            printf("CUDA error: %s [%s]\n", errorMessage, cudaGetErrorString(err));
        }
        err = cudaThreadSynchronize();
        if(cudaSuccess != err)
        {
            printf("CUDA error: %s [%s]\n", errorMessage, cudaGetErrorString(err));
        }
    }

    //PREPROCESSOR MACROS
    //PACK_TO_LEFT pack a signed short number to more significant two bytes of int
    #define PACK_TO_LEFT(x)         (((int)((unsigned short)(x))) << 16)
    #define PACK_TO_RIGHT(x)        ((int)((unsigned short)(x)))
    #define PACK(x, y)              (PACK_TO_LEFT(x)|PACK_TO_RIGHT(y))

    #define PACK_0_BYTE(x)          (((int)((unsigned char)(x))) << 24)
    #define PACK_1_BYTE(x)          (((int)((unsigned char)(x))) << 16)
    #define PACK_2_BYTE(x)          (((int)((unsigned char)(x))) << 8)
    #define PACK_3_BYTE(x)          (((int)((unsigned char)(x))))
    #define PACK_BYTES(a,b,c,d)     (PACK_0_BYTE(a)|PACK_1_BYTE(b)|PACK_2_BYTE(c)|PACK_3_BYTE(d))  



    #define UNPACK_FROM_LEFT(x)     ((short)((x) >> 16))
    #define UNPACK_FROM_RIGHT(x)    ((short)((x)&0xFFFF))

    #define UNPACK_0_BYTE(x)        ((char)( (x) >> 24))
    #define UNPACK_1_BYTE(x)        ((char)(((x) >> 16)&0xFF))
    #define UNPACK_2_BYTE(x)        ((char)(((x) >>  8)&0xFF))
    #define UNPACK_3_BYTE(x)        ((char)( (x)       &0xFF))


    #define EXCHANGE_UP_LEFT(x, y)  x = (((y)&0xFFFF0000) | (x >> 16));
    //#define CHANGE_UP_LEFT(x, y)    x = ((x&0xFFFF0000) | ((unsigned short)(y)));



        
    #define signum(x) (signbit(-(x)) - signbit(x))

    //CONSTANTS
    #define SHORT_MIN -32768

#endif
