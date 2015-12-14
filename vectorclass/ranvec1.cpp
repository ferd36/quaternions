/****************************  ranvec1.cpp   **********************************
* Author:        Agner Fog
* Date created:  2014-09-09
* Last modified: 2014-10-17
* Version:       1.15
* Project:       vector classes
* Description:
* Header file defining pseudo random number generators with vector output.
* Two pseudo random number generators are combined:
* 1. "Mersenne Twister for Graphic Processor" (MTGP).
*    (Saito & Matsumoto: Variants of Mersenne Twister Suitable for Graphic
*    Processors". ACM Transactions on Mathematical Software, v. 39, no. 2,
*    2013).
* 2. "Multiply-With-Carry Generator" (MWC).
*    (Press, et. al.: Numerical Recipes: The Art of Scientific Computing,
*    3rd. edition. Cambridge Univ. Press, 2007).
*
* Instructions:
* Make an object of the class Ranvec1. The constructor has a parameter "gtype"
* to indicate the desired random number generator:
* gtype = 1: MWC.  Use for smaller projects and where speed is critical
* gtype = 2: MTGP. Use for large projects with many streams
* gtype = 3: Both generators combined. Use for the most demanding projects
* 
* Multi-threaded programs must make one instance of Ranvec1 for each thread,
* with different seeds. It is not safe to access the same random number 
* generator instance in multiple threads. 
*
* The Ranvec1 object must be initialized with one or more seeds, by calling
* one of the init functions. The same seed will always generate the same 
* sequence of random numbers (with the same value of gtype). A different seed
* or set of seeds will generate a different sequence. You can use one of the
* following member functions for initialization:
*
* void init(int seed):  General initialization function
* void init(int seed1, int seed2):  Use seed1 for the MWC generator and
*        seed2 for the MTGP.
* void initByArray(int seeds[], int numSeeds):  Use an array of seeds.
*        The sequence will change if at least one of the seeds is changed.
*        If gtype = 3 then seeds[0] will be used for the MWC generator and
*        all the remaining seeds will be used for MTGP.
*
* The following member functions can be used for random number outputs:
* Scalars:
* uint32_t random32b():                 Returns an integer of 32 random bits
* int      random1i(int min, int max):  One integer in the interval min <= x <= max
* int      random1ix(int min, int max): Same, with extra precision
* uint64_t random64b():                 Returns an integer of 64 random bits
* float    random1f():                  One floating point number in the interval 0 <= x < 1
* double   random1d():                  One double in the interval 0 <= x < 1
*
* 128 bit vectors:
* Vec4ui   random128b():                Returns 128 random bits as a vector of 4 integers
* Vec4i    random4i(int min, int max):  4 integers in the interval min <= x <= max
* Vec4i    random4ix(int min, int max): Same, with extra precision
* Vec4f    random4f():                  4 floating point numbers in the interval 0 <= x < 1
* Vec2d    random2d():                  2 doubles in the interval 0 <= x < 1
*
* 256 bit vectors:
* Vec8ui   random256b():                Returns 256 random bits as a vector of 8 integers
* Vec8i    random8i(int min, int max):  8 integers in the interval min <= x <= max
* Vec8i    random8ix(int min, int max): Same, with extra precision
* Vec8f    random8f():                  8 floating point numbers in the interval 0 <= x < 1
* Vec4d    random4d():                  4 doubles in the interval 0 <= x < 1
*
* 512 bit vectors:
* Vec16ui  random512b():                Returns 512 random bits as a vector of 16 integers
* Vec16i   random16i(int min, int max): 16 integers in the interval min <= x <= max
* Vec16i   random16ix(int min, int max):Same, with extra precision
* Vec16f   random16f():                 16 floating point numbers in the interval 0 <= x < 1
* Vec8d    random8d():                  8 doubles in the interval 0 <= x < 1
*
* The 256 bit vector functions are available only if MAX_VECTOR_SIZE >= 256.
* The 512 bit vector functions are available only if MAX_VECTOR_SIZE >= 512.
*
* For detailed instructions, see VectorClass.pdf
* For theoretical explanation, see the article: "Pseudo-Random Number Generators
* for Vector Processors and Multicore processors". www.agner.org/random/theory
*
* (c) Copyright 2014 GNU General Public License www.gnu.org/licenses
******************************************************************************/

#include "ranvec1.h"


/******************************************************************************
                      Member functions for Ranvec1base
******************************************************************************/

// Constructor
Ranvec1base::Ranvec1base(int gtype) {
    // Set everything to zero
    // memset(this, 0, sizeof(*this));
    int i;
    for (i=0; i<sizeof(buffer1)/sizeof(*buffer1); i++) buffer1[i] = 0;
    for (i=0; i<sizeof(buffer2)/sizeof(*buffer2); i++) buffer2[i] = 0;
    iw = 0;
    idx = idm = 0;
    nextx = xj = 0;
    if (gtype > 3) gtype = 3;
    if (gtype < 1) gtype = 1;
    gentype = gtype;                               // Save generator type
}

// Initialize with one seed
void Ranvec1base::init(int seed) {
    switch (gentype) {
    case 1:                                      // MWC only
        initMWC(seed);
        break;
    case 2:                                      // MTGP only
        initMTGP(seed);
        break;
    case 3:                                      // Both
        initMWC(seed);
        initMTGP(seed);
        break;
    }
}

// Initialize with two seeds
void Ranvec1base::init(int seed1, int seed2) {
    int seeds[2] = {seed1, seed2};
    switch (gentype) {
    case 1: case 2:
        initByArray(seeds, 2);
        break;
    case 3:
        initMWC(seed1);
        initMTGP(seed2);
        break;
    }
}

// Initialize by array of seeds
void Ranvec1base::initByArray(int32_t const seeds[], int numSeeds) {
    if (numSeeds < 1) init(0);                   // Error if no seeds
    int i;
    switch (gentype) {
    case 1:
        initMWC(0);                              // Just for setting iw = 0 or whatever may be needed in future implementations
        // Initialize buffer2 from array
        initMTGPByArray(seeds, numSeeds);
        // Copy to buffer1
        for (i = 0; i < sizeof(buffer1)/sizeof(buffer1[0]); i++) {
            buffer1[i] = buffer2[i];             // Copy to buffer1
        } 
        // Randomize one round
        for (i = 0; i < 4 * 16 / sizeof(next1()); i++)  next1();
        // Check for illegal seeds
        for (i = 0; i < sizeof(buffer1)/sizeof(buffer1[0]); i+=2) {
            static const uint32_t factorsm1[8] = {
                mwcfac0-1, mwcfac1-1, mwcfac2-1, mwcfac3-1, mwcfac4-1, mwcfac5-1, mwcfac6-1, mwcfac7-1};
            if ((buffer1[i]==0 && buffer1[i+1]==0) || (buffer1[i]==0xFFFFFFFF && buffer1[i+1]==factorsm1[i/2])) {
                // seeds (0,0) and (-1,factor-1) are avoided because they will reproduce themselves
                buffer1[i]   = (buffer2[i+16] & -4) | 1;
                buffer1[i+1] = buffer2[i+17];
            }
        }
        break;
    case 2:
        initMTGPByArray(seeds, numSeeds);
        break;
    case 3:
        if (numSeeds > 1) {
            initMWC(seeds[0]);
            initMTGPByArray(seeds+1, numSeeds-1);
        }
        else {
            init(seeds[0]);
        }
        break;
    }
}

// Get next 512 bits into a buffer
void Ranvec1base::next(uint32_t * dest) {
#if INSTRSET < 8                                 // SSE2 - AVX: use 128 bit vectors
    int i;
    switch (gentype) {
    case 1:                                      // MWC
        for (i = 0; i < 4; i++) {
            next1().store(dest + 4*i);
        }
        break;
    case 2:                                      // MTGP
        for (i = 0; i < 4; i++) {
            next2().store(dest + 4*i);
        }
        break;
    case 3:                                      // Both
        for (i = 0; i < 4; i++) {
            (next1()+next2()).store(dest + 4*i);
        }
        break;
    }
#elif INSTRSET < 9                               // AVX2: use 256 bit vectors
    switch (gentype) {
    case 1:                                      // MWC
        next1().store(dest);
        next1().store(dest+8);
        break;
    case 2:                                      // MTGP
        next2().store(dest);
        next2().store(dest+8);
        break;
    case 3:                                      // Both
        (next1()+next2()).store(dest);
        (next1()+next2()).store(dest+8);
        break;
    }
#else                                            // AVX512: use 512 bit vectors
    switch (gentype) {
    case 1:                                      // MWC
        next1().store(dest);
        break;
    case 2:                                      // MTGP
        next2().store(dest);
        break;
    case 3:                                      // Both
        (next1()+next2()).store(dest);
        break;
    }
#endif
}


/******************************************************************************
                      Member functions for Ranvec1base, MWC generator
******************************************************************************/

// Initialize with seed
void Ranvec1base::initMWC(int seed) {
    static const int vectorsize = sizeof(next1());         // Vector size depends on instruction set
    int i;
    // Fill buffer with function of seed
    uint32_t tmp = seed;

    // Seeds (0,0) and (-1,factor-1) will degenerate into constant output.
    // This seeding method should avoid these forbidden seeds:

    for (i = 0; i < 16; i++) {
        tmp = buffer1[i] = 1566083941u * (tmp ^ (tmp >> 27)) + i;
    }
    // Randomize 4 rounds
    for (i = 0; i < 4 * 64 / vectorsize; i++)  next1();
    // Initialize index
    iw = 0;
}

// Next state. Generate next vector and output
#if INSTRSET < 8                                           // SSE2 - AVX: use 128 bit vectors
Vec4ui Ranvec1base::next1() {                              // Get 128 bits from MWC
    // Factors for multiply-with-carry
    static const uint32_t factors[16] = {
        mwcfac0, 0, mwcfac1, 0, mwcfac2, 0, mwcfac3, 0,
        mwcfac4, 0, mwcfac5, 0, mwcfac6, 0, mwcfac7, 0
    };
    Vec4ui x, f;
    Vec2uq y;
    x.load(buffer1 + iw);                                  // previous x and carry
    f.load(factors + iw);                                  // factors
    y = _mm_mul_epu32(x, f);                               // 32*32->64 bit unsigned multiply
    y += Vec2uq(x) >> 32u;                                 // add old carry
    y.store(buffer1 + iw);                                 // new x and carry
    // The uppermost bits of the carry are not sufficiently random. Randomize some more for output
    y ^= y << shw1;                                        // output function
    y ^= y >> shw2;
    y ^= y << shw3;
    iw = (iw + 4) & 15;
    return Vec4ui(y);
}

#elif INSTRSET < 9                                         // AVX2: use 256 bit vectors
Vec8ui Ranvec1base::next1() {                              // Get 128 bits from MWC
    // Factors for multiply-with-carry
    static const uint32_t factors[16] = {
        mwcfac0, 0, mwcfac1, 0, mwcfac2, 0, mwcfac3, 0,
        mwcfac4, 0, mwcfac5, 0, mwcfac6, 0, mwcfac7, 0
    };
    Vec8ui x, f;
    Vec4uq y;
    x.load(buffer1 + iw);                                  // previous x and carry
    f.load(factors + iw);                                  // factors
    y = _mm256_mul_epu32(x, f);                            // 32*32->64 bit unsigned multiply
    y += Vec4uq(x) >> 32u;                                 // add old carry
    y.store(buffer1 + iw);                                 // new x and carry
    // The uppermost bits of the carry are not sufficiently random. Randomize some more for output
    y ^= y << shw1;                                        // output function
    y ^= y >> shw2;
    y ^= y << shw3;
    iw = (iw + 8) & 15;
    return Vec8ui(y);
}
#else                                                      // AVX512: use 512 bit vectors
Vec16ui Ranvec1base::next1() {                             // Get 512 bits from MWC
    // Factors for multiply-with-carry
    static const uint32_t factors[16] = {
        mwcfac0, 0, mwcfac1, 0, mwcfac2, 0, mwcfac3, 0,
        mwcfac4, 0, mwcfac5, 0, mwcfac6, 0, mwcfac7, 0
    };
    Vec16ui x, f;
    Vec8uq y;
    x.load(buffer1);                                       // previous x and carry
    f.load(factors);                                       // factors
    y = _mm512_mul_epu32(x, f);                            // 32*32->64 bit unsigned multiply
    y += Vec8uq(x) >> 32u;                                 // add old carry
    y.store(buffer1);                                      // new x and carry
    // The uppermost bits of the carry are not sufficiently random. Randomize some more for output
    y ^= y << shw1;                                        // output function
    y ^= y >> shw2;
    y ^= y << shw3;
    return Vec16ui(y);
}
#endif



/******************************************************************************
                      Member functions for Ranvec1base, MTGP generator
******************************************************************************/

// Initialize MTGP with seed
void Ranvec1base::initMTGP(int seed) {
    int i;
    const uint32_t hidden_seed =  tbl2 ^ (tbl3 << 16);
    uint32_t tmp = hidden_seed;
    tmp += tmp >> 16;
    tmp += tmp >> 8;
    // Fill buffer
    tmp = (tmp & 0xFF) * 0x01010101;
    for (i=0; i<sizeof(buffer2)/sizeof(*buffer2); i++) buffer2[i] = tmp;
    // memset(buffer2, (uint8_t)tmp, sizeof(buffer2));
    // Seed buffer
    tmp = buffer2[bo+0] = seed;
    buffer2[bo+1] = hidden_seed;    
    for (i = 1; i < bsize; i++) {	
        tmp = buffer2[bo+i] ^= 1812433253u * (tmp ^ (tmp >> 30)) + i;
    }
    // Copy beginning of buffer to end for wrap-around
    xj.load(buffer2+bo);
    xj.store(buffer2+bo+bsize);

    // Initialize indexes
    idx = 0;
    idm = mpos;
    nextx.load(buffer2+bo + idx);
    xj.load(buffer2+bo + idm - sizeof(xj)/4);
}


// Function used in init_by_array
static uint32_t ini_func1(uint32_t x) {
    return (x ^ (x >> 27)) * 1664525u;
}

// Function used in init_by_array
static uint32_t ini_func2(uint32_t x) {
    return (x ^ (x >> 27)) * 1566083941u;
}

// Initialize by array of seeds
void Ranvec1base::initMTGPByArray(int32_t const seeds[], int numSeeds) {
    int i, j, count;
    uint32_t r;
    int lag;
    int mid;
    uint32_t hidden_seed;
    uint32_t tmp;
    uint32_t * buf0 = buffer2 + bo;

    if (bsize >= 623)      lag = 11;
    else if (bsize >= 68)  lag = 7;
    else if (bsize >= 39)  lag = 5;
    else                   lag = 3;

    mid = (bsize - lag) / 2;
    hidden_seed = tbl2 ^ (tbl3 << 16);
    tmp = hidden_seed;
    tmp += tmp >> 16;
    tmp += tmp >> 8;
    // memset(buffer2, (uint8_t)tmp, sizeof(buffer2));
    tmp = (tmp & 0xFF) * 0x01010101;
    for (i=0; i<sizeof(buffer2)/sizeof(*buffer2); i++) buffer2[i] = tmp;
    buf0[0] = hidden_seed;

    if (numSeeds > bsize - 1)  count = numSeeds;
    else                       count = bsize - 1;

    r = ini_func1(buf0[0] ^ buf0[mid] ^ buf0[bsize - 1]);
    buf0[mid] += r;
    r += numSeeds;
    buf0[(mid + lag) % bsize] += r;
    buf0[0] = r;
    i = 1;
    for (j = 0; j < count; j++) {
        r = ini_func1(buf0[i] ^ buf0[(i + mid) % bsize] ^ buf0[(i + bsize - 1) % bsize]);
        buf0[(i + mid) % bsize] += r;
        if (j < numSeeds) r += seeds[j];
        r += i; 
        buf0[(i + mid + lag) % bsize] += r;
        buf0[i] = r;
        i = (i + 1) % bsize;
    }
    for (j = 0; j < bsize; j++) {
        r = ini_func2(buf0[i] + buf0[(i + mid) % bsize] + buf0[(i + bsize - 1) % bsize]);
        buf0[(i + mid) % bsize] ^= r;
        r -= i;
        buf0[(i + mid + lag) % bsize] ^= r;
        buf0[i] = r;
        i = (i + 1) % bsize;
    }
    const uint32_t non_zero = 0x4d544750;
    if (buf0[bsize - 1] == 0) buf0[bsize - 1] = non_zero;

    // Copy beginning of buffer to end for wrap-around
    xj.load(buffer2+bo);
    xj.store(buffer2+bo+bsize);

    // Initialize indexes
    idx = 0;
    idm = mpos;
    nextx.load(buffer2+bo + idx);
    nextx.load(buffer2+bo + idx);
    xj.load(buffer2+bo + idm - sizeof(xj)/4);
}


// Next state. Return next vector from MTGP
#if INSTRSET < 8                                 // SSE2 - AVX: use 128 bit vectors
Vec4ui Ranvec1base::next2() {
    Vec4ui x, y, lastxj;
    int nextidx = idx + 4; 
    if (nextidx >= bsize) nextidx -= bsize;

    x = nextx;
    nextx.load(buffer2+bo + nextidx);
    y = blend4ui<1,2,3,4>(x, nextx);             // (This will use PALIGNR instruction if SSSE3 available. Inefficient for lower instruction sets)
    x = (x & mask) ^ y;
    lastxj = xj;
    y = xj.load(buffer2+bo + idm);
    x ^= x << sh1;
    y = x ^ (y >> sh2);

   // 4x32 matrix multiplication one column at a time
    Vec4i m0 = Vec4i(y) << 31 >> 31;             // copy bit 0 of each element
    Vec4i m1 = Vec4i(y) << 30 >> 31;             // copy bit 1 of each element
    Vec4i m2 = Vec4i(y) << 29 >> 31;             // copy bit 2 of each element
    Vec4i m3 = Vec4i(y) << 28 >> 31;             // copy bit 3 of each element
    Vec4i m = ((m0 & tbl0) ^ (m1 & tbl1)) ^ ((m2 & tbl2) ^ (m3 & tbl3));
    y ^= Vec4ui(m);

    y.store(buffer2+bo + idx);
    // Wrap-around: Store both at beginning and end of buffer.
    // Do this early to avoid store forwarding stall on subsequent read with different alignment
    if (idx < 4) {
        y.store(buffer2+bo + bsize + idx);       // Store also at end
    }
    else if (idx >= bsize - 3) {
        y.store(buffer2+bo - bsize + idx);       // Store also at beginning
    }

    // Tempering of output
    Vec4ui t1 = blend4ui<3,4,5,6>(lastxj, xj);
    t1 ^= t1 >> 16;
    t1 ^= t1 >> 8;

   // 4x32 matrix multiplication one column at a time
    Vec4i k0 = Vec4i(t1) << 31 >> 31;            // copy bit 0 of each element
    Vec4i k1 = Vec4i(t1) << 30 >> 31;            // copy bit 1 of each element
    Vec4i k2 = Vec4i(t1) << 29 >> 31;            // copy bit 2 of each element
    Vec4i k3 = Vec4i(t1) << 28 >> 31;            // copy bit 3 of each element
    Vec4i k = ((k0 & temper0) ^ (k1 & temper1)) ^ ((k2 & temper2) ^ (k3 & temper3));
    y ^= Vec4ui(k);

    // Advance pointers, with wrap-around
    idx = nextidx;
    idm += 4; 
    if (idm >= bsize) idm -= bsize;

    // Return vector
    return y;
}
#elif INSTRSET < 9                               // AVX2: use 256 bit vectors

// Get 256 bits from MTGP
Vec8ui Ranvec1base::next2() {
    Vec8ui x, y, lastxj;
    int nextidx = idx + 8; 
    if (nextidx >= bsize) nextidx -= bsize;

    x = nextx;
    nextx.load(buffer2+bo + nextidx);
    y = blend8ui<1,2,3,4,5,6,7,8>(x, nextx);
    x = (x & mask) ^ y;
    lastxj = xj;
    y = xj.load(buffer2+bo + idm);
    x ^= x << sh1;
    y = x ^ (y >> sh2);

    // 4x32 matrix multiplication by table lookup
    static const uint32_t tbl[16] = {
        0, tbl0, tbl1, tbl1 ^ tbl0,
        tbl2, tbl2 ^ tbl0, tbl2 ^ tbl1, tbl2 ^ tbl1 ^ tbl0,
        tbl3, tbl3 ^ tbl0, tbl3 ^ tbl1, tbl3 ^ tbl1 ^ tbl0,
        tbl3 ^ tbl2, tbl3 ^ tbl2 ^ tbl0, tbl3 ^ tbl2 ^ tbl1, tbl3 ^ tbl2 ^ tbl1 ^ tbl0};
    Vec8ui m = Vec8ui(lookup<16>(Vec8i(y & 0x0F), tbl));   // Table lookup. (will use permute which is faster than gather)
    y ^= m;

    y.store(buffer2+bo + idx);
    // Wrap-around: Store both at beginning and end of buffer.
    // Do this early to avoid store forwarding stall on subsequent read with different alignment
    if (idx < 8) {
        y.store(buffer2+bo + bsize + idx);       // Store also at end
    }
    else if (idx >= bsize - 7) {
        y.store(buffer2+bo - bsize + idx);       // Store also at beginning
    }

    // Tempering of output
    Vec8ui t1 = blend8ui<7,8,9,10,11,12,13,14>(lastxj, xj);
    t1 ^= t1 >> 16;
    t1 ^= t1 >> 8;

    // 4x32 matrix multiplication by table lookup
    static const uint32_t tempmat[16] = {
        0, temper0, temper1, temper1 ^ temper0,
        temper2, temper2 ^ temper0, temper2 ^ temper1, temper2 ^ temper1 ^ temper0,
        temper3, temper3 ^ temper0, temper3 ^ temper1, temper3 ^ temper1 ^ temper0,
        temper3 ^ temper2, temper3 ^ temper2 ^ temper0, temper3 ^ temper2 ^ temper1, temper3 ^ temper2 ^ temper1 ^ temper0};
    m = Vec8ui(lookup<16>(Vec8i(t1 & 0x0F), tempmat));
    y ^= m;

    // Advance pointers, with wrap-around
    idx = nextidx;
    idm += 8; 
    if (idm >= bsize) idm -= bsize;

    // Return vector
    return y;
}
#else                                            // AVX512: use 512 bit vectors
// Get 512 bits from MTGP
Vec16ui Ranvec1base::next2() {
    Vec16ui x, y, lastxj;
    int nextidx = idx + 16; 
    if (nextidx >= bsize) nextidx -= bsize;

    x = nextx;
    nextx.load(buffer2+bo + nextidx);
    y = blend16ui<1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16>(x, nextx);
    x = (x & mask) ^ y;
    lastxj = xj;
    y = xj.load(buffer2+bo + idm);
    x ^= x << sh1;
    y = x ^ (y >> sh2);

    // 4x32 matrix multiplication by table lookup
    static const uint32_t tbl[16] = {
        0, tbl0, tbl1, tbl1 ^ tbl0,
        tbl2, tbl2 ^ tbl0, tbl2 ^ tbl1, tbl2 ^ tbl1 ^ tbl0,
        tbl3, tbl3 ^ tbl0, tbl3 ^ tbl1, tbl3 ^ tbl1 ^ tbl0,
        tbl3 ^ tbl2, tbl3 ^ tbl2 ^ tbl0, tbl3 ^ tbl2 ^ tbl1, tbl3 ^ tbl2 ^ tbl1 ^ tbl0};
    Vec16ui m = Vec16ui(lookup<16>(Vec16i(y & 0x0F), tbl));
    y ^= m;

    y.store(buffer2+bo + idx);
    // Wrap-around: Store both at beginning and end of buffer.
    // Do this early to avoid store forwarding stall on subsequent read with different alignment
    if (idx < 16) {
        y.store(buffer2+bo + bsize + idx);       // Store also at end
    }
    else if (idx >= bsize - 15) {
        y.store(buffer2+bo - bsize + idx);       // Store also at beginning
    }

    // Tempering of output
    Vec16ui t1 = blend16ui<15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30>(lastxj, xj);
    t1 ^= t1 >> 16;
    t1 ^= t1 >> 8;

    // 4x32 matrix multiplication by table lookup
    static const uint32_t tempmat[16] = {
        0, temper0, temper1, temper1 ^ temper0,
        temper2, temper2 ^ temper0, temper2 ^ temper1, temper2 ^ temper1 ^ temper0,
        temper3, temper3 ^ temper0, temper3 ^ temper1, temper3 ^ temper1 ^ temper0,
        temper3 ^ temper2, temper3 ^ temper2 ^ temper0, temper3 ^ temper2 ^ temper1, temper3 ^ temper2 ^ temper1 ^ temper0};
    m = Vec16ui(lookup<16>(Vec16i(t1 & 0x0F), tempmat));
    y ^= m;

    // Advance pointers, with wrap-around
    idx = nextidx;
    idm += 16; 
    if (idm >= bsize) idm -= bsize;

    // Return vector
    return y;
}
#endif


/******************************************************************************
                      Member functions for Ranvec1
******************************************************************************/

// Reset all output buffers
void Ranvec1::resetBuffers() {
    buf32.reset();
    buf64.reset();
    buf128.reset();
#if MAX_VECTOR_SIZE >= 256
    buf256.reset();
#endif
#if MAX_VECTOR_SIZE >= 512
    buf512.reset();
#endif
}


/******************************************************************************
                      Scalar output functions for Ranvec1
******************************************************************************/

// One integer in the interval min <= x <= max
// Relative error on frequencies < 2^-32
int Ranvec1::random1i(int min, int max) {
   if (max <= min) {
      if (max == min) {
          return min;                            // Interval has only one value
      }
      else {
          return 0x80000000;                     // Error: interval length is negative
      }
   }
   uint32_t interval;                            // Length of interval
   uint64_t longran;                             // Random bits * interval
   uint32_t iran;                                // Longran / 2^32
   uint32_t ranbits = random32b();               // Random bits
   interval = (uint32_t)(max - min) + 1u;        // Length of interval
   if (interval == 0) return ranbits;            // interval overflows
   longran  = (uint64_t)ranbits * interval;
   iran = (uint32_t)(longran >> 32);
   // Convert back to signed and return result
   return (int32_t)iran + min;
}

// Same, with extra precision
int Ranvec1::random1ix(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Each output value has exactly the same probability.
   // This is obtained by rejecting certain bit values so that the number
   // of possible bit values is divisible by the interval length
   if (max <= min) {
      if (max == min) {
          return min;                            // Interval has only one value
      }
      else {
          return 0x80000000;                     // Error: interval length is negative
      }
   }
   // Assume that 64 bit integers are available. Use multiply and shift method
   uint32_t interval;                            // Length of interval
   uint64_t longran;                             // Random bits * interval
   uint32_t iran;                                // longran / 2^32
   uint32_t remainder;                           // longran % 2^32

   interval = uint32_t(max - min) + 1u;          // length of interval
   if (interval == 0) return random32b();        // avoid division by 0
   if (interval != randomixInterval) {
      // Interval length has changed. Must calculate rejection limit.
      // Reject when remainder >= 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then
      randomixLimit = uint32_t(((uint64_t)1 << 32) / interval) * interval - 1;
      randomixInterval = interval;
   }
   do { // Rejection loop
      longran  = (uint64_t)random32b() * interval;
      iran = (uint32_t)(longran >> 32);
      remainder = (uint32_t)longran;
   } while (remainder > randomixLimit);
   // Convert back to signed and return result
   return (int32_t)iran + min;
}

// One floating point number in the interval 0 <= x < 1
// Resolution is extended to 2^-24
float Ranvec1::random1f() {
    union Uif {
        uint32_t i;
        float f;
    };
    Uif u1, u2;
    uint32_t r = random32b();                    // get 32 random bits
    // Insert exponent and random mantissa to get random number in the interval 1 <= x < 2
    // Subtract 1.0 if next bit is 0, or 1.0 - 2^-24 = 0.99999994f if next bit is 1
    u1.i = 0x3F800000 - ((r >> 8) & 1);          // bit 8
    u2.i = (r >> 9) | 0x3F800000;                // bit 9 - 31
    return u2.f - u1.f;
    //return u2.f - (r & 0x100 ? 0.99999994f : 1.0f);
}

// One double in the interval 0 <= x < 1
// Resolution is 2^-52
double Ranvec1::random1d() {
    union Uqd {
        uint64_t q;
        double d;
    };
    Uqd u1;
    uint64_t r = random64b();                    // get 64 random bits
    // Insert exponent and random mantissa to get random number in the interval 1 <= x < 2,
    // then subtract 1.0 to get the interval 0 <= x < 1.
    // A resolution of 2^-52 is sufficient. No need to put in an extra bit
    u1.q = (r >> 12) | 0x3FF0000000000000;       // bit 12 - 63
    return u1.d - 1.0;
}

/******************************************************************************
          Helper functions for random vector output functions
******************************************************************************/
// 32x32->64 bit unsigned multiplication
static inline Vec4ui mulExtended(Vec4ui const & a, Vec4ui const & b) {
    return _mm_mul_epu32(a, b);
}

// Shift down 32 bits
static inline Vec4ui shift32Down(Vec4ui const & a) {
    return Vec4ui(Vec2uq(a) >> 32u);
}

// Shift up 32 bits
static inline Vec4ui shift32Up(Vec4ui const & a) {
    return Vec4ui(Vec2uq(a) << 32u);
}

#if MAX_VECTOR_SIZE >= 256
// 32x32->64 bit unsigned multiplication
static inline Vec8ui mulExtended(Vec8ui const & a, Vec8ui const & b) {
#if INSTRSET >= 8  // AVX2
    return _mm256_mul_epu32(a, b);
#else
    return Vec8ui(mulExtended(a.get_low(), b.get_low()), mulExtended(a.get_high(), b.get_high()));
#endif
}

// Shift down 32 bits
static inline Vec8ui shift32Down(Vec8ui const & a) {
    return Vec8ui(Vec4uq(a) >> 32u);
}

// Shift up 32 bits
static inline Vec8ui shift32Up(Vec8ui const & a) {
    return Vec8ui(Vec4uq(a) << 32u);
}
#endif

#if MAX_VECTOR_SIZE >= 512
// 32x32->64 bit unsigned multiplication
static inline Vec16ui mulExtended(Vec16ui const & a, Vec16ui const & b) {
#if INSTRSET >= 9  // AVX512f
    return _mm512_mul_epu32(a, b);
#else
    return Vec16ui(mulExtended(a.get_low(), b.get_low()), mulExtended(a.get_high(), b.get_high()));
#endif
}

// Shift down 32 bits
static inline Vec16ui shift32Down(Vec16ui const & a) {
    return Vec16ui(Vec8uq(a) >> 32u);
}

// Shift up 32 bits
static inline Vec16ui shift32Up(Vec16ui const & a) {
    return Vec16ui(Vec8uq(a) << 32u);
}
#endif


/******************************************************************************
                      128 bit vector output functions for Ranvec1
******************************************************************************/
// 4 integers in the interval min <= x <= max
// Relative error on frequencies < 2^-32   
Vec4i Ranvec1::random4i(int min, int max) {
   if (max <= min) {
      if (max == min) {
          return Vec4i(min);                                                   // Interval has only one value
      }
      else {
          return Vec4i(0x80000000);                                            // Error: interval length is negative
      }
   }
   Vec4ui ranbits   = random128b();                                            // Random bits
   uint32_t interval = (uint32_t)(max - min) + 1u;                             // Length of interval
   if (interval == 0) return ranbits;
   Vec4ui prod_even = mulExtended(ranbits, interval);                          // 64-bit product of even numbered elements
   Vec4ui prod_odd  = mulExtended(shift32Down(ranbits), interval);             // 64-bit product of odd numbered elements
   static const Vec4ib select_odd(false,true,false,true);                      // Odd elements true
   Vec4ui prods     = select(select_odd, prod_odd,shift32Down(prod_even));     // combine high parts of products
   return Vec4i(prods) + Vec4i(min);
}

// 4 integers in the interval min <= x <= max, with extra precision
Vec4i Ranvec1::random4ix(int min, int max) {
   if (max <= min) {
      if (max == min) {
          return Vec4i(min);                                                   // Interval has only one value
      }
      else {
          return Vec4i(0x80000000);                                            // Error: interval length is negative
      }
   }
   uint32_t interval = (uint32_t)(max - min) + 1u;
   if (interval == 0) return random128b();                                     // avoid division by 0
   if (interval != randomixInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder >= 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then
      randomixLimit = uint32_t(((uint64_t)1 << 32) / interval) * interval - 1;
      randomixInterval = interval;
   }
   Vec4ui prods1;                                                              // Random number * interval
   Vec4ib reject1(true);                                                       // Initial values in prods1 are invalid
   static const Vec4ib select_odd(false,true,false,true);                      // Odd elements true

   bool incomplete;                                                            // when not all elements are accepted
   int n = 0;                                                                  // Counter for rejection loop

   // Rejection loop
   do {
       Vec4ui ranbits   = random128b();                                        // Random bits
       Vec4ui prod_even = mulExtended(ranbits, Vec4ui(interval));              // 64-bit product of even numbered elements
       Vec4ui prod_odd  = mulExtended(shift32Down(ranbits), Vec4ui(interval)); // 64-bit product of odd  numbered elements
       Vec4ui prods2    = select(select_odd, prod_odd, shift32Down(prod_even));// combine high parts of products
       Vec4ui remainders= select(select_odd, shift32Up(prod_odd), prod_even);  // combine low  parts of products
       Vec4ib reject2   = remainders > randomixLimit;                          // Reject these elements
       Vec4ib transfer  = andnot(reject1, reject2);                            // Transfer elements that are missing in prods1 and accepted in prods2
       prods1  = select(transfer, prods2, prods1);
       reject1 = reject1 & reject2;                                            // Indicates elements still missing
       incomplete = horizontal_or(reject1);                                    // True if any elements missing
       n++;
       if (n > 1 && incomplete) {
           // Still some elements missing after 2 or more iterations. Look for accepted elements that can be moved
           reject2 |= transfer;                                                // These elements have already been used
           int i, j = 0;
           for (i = 0; i < 4; i++) {                                           // Loop through missing elements in prods1
               if (reject1[i]) {
                   while (j < 4) {                                             // Search for accepted unused elements in prods2
                       if (!reject2[j]) {
                           prods1.insert(i, prods2[j]);                        // Move element from prods2[j] to prods1[i]
                           reject1.insert(i, false);                           // Remember this element is now accepted
                           j++;
                           break;
                       }
                       j++;
                   }
               }
           }
           incomplete = horizontal_or(reject1);                                // Stop if all elements have been fixed now
       }
   }
   while (incomplete);

   return Vec4i(prods1) + Vec4i(min);
}

// 4 floating point numbers with uniform distribution in the interval 0 <= x < 1.
// Extra bit inserted for resolution 2^-24
Vec4f Ranvec1::random4f() {
    Vec4ui const one = 0x3F800000;                                             // Binary representation of 1.0f
    Vec4ui ranbits = random128b();                                             // Get random bits
    Vec4ui r1 = one - ((ranbits >> 8) & 1);                                    // 1.0 if bit8 is 0, or 1.0-2^-24 = 0.99999994f if bit8 is 1
    Vec4ui r2 = (ranbits >> 9) | one;                                          // bits 9 - 31 inserted as mantissa
    return Vec4f(reinterpret_f(r2)) - Vec4f(reinterpret_f(r1));                // Get into interval 0 <= x < 1
}

// 2 doubles with uniform distribution in the interval 0 <= x < 1
// Resolution 2^-52. (One more bit is possible, but this would be overkill)
Vec2d Ranvec1::random2d() {
    Vec2uq const one = 0x3FF0000000000000;                                     // Binary representation of 1.0
    Vec2uq ranbits = Vec2uq(random128b());                                     // Get random bits
    Vec2uq r = (ranbits >> 12) | one;                                          // bits 12 - 63 inserted as mantissa
    return Vec2d(reinterpret_d(r)) - Vec2d(reinterpret_d(one));                // Get into interval 0 <= x < 1
}


/******************************************************************************
                      256 bit vector output functions for Ranvec1
******************************************************************************/

#if MAX_VECTOR_SIZE >= 256
// 8 integers in the interval min <= x <= max
// Relative error on frequencies < 2^-32   
Vec8i Ranvec1::random8i(int min, int max) {
   if (max <= min) {
      if (max == min) {
          return Vec8i(min);                                                   // Interval has only one value
      }
      else {
          return Vec8i(0x80000000);                                            // Error: interval length is negative
      }
   }
   Vec8ui ranbits   = random256b();                                            // Random bits
   uint32_t interval = (uint32_t)(max - min) + 1u;                             // Length of interval
   if (interval == 0) return ranbits;                                          // interval overflows
   Vec8ui prod_even = mulExtended(ranbits, interval);                          // 64-bit product of even numbered elements
   Vec8ui prod_odd  = mulExtended(shift32Down(ranbits), interval);             // 64-bit product of odd numbered elements
   static const Vec8ib select_odd(false,true,false,true,false,true,false,true);// Odd elements true
   Vec8ui prods     = select(select_odd, prod_odd,shift32Down(prod_even));     // combine high parts of products
   return Vec8i(prods) + Vec8i(min);
}

// 8 integers in the interval min <= x <= max, with extra precision
// (Not guaranteed to give exactly the same sequence as random4ix called twice)
Vec8i Ranvec1::random8ix(int min, int max) {
   if (max <= min) {
      if (max == min) {
          return Vec8i(min);                                                   // Interval has only one value
      }
      else {
          return Vec8i(0x80000000);                                            // Error: interval length is negative
      }
   }
   uint32_t interval = (uint32_t)(max - min) + 1u;
   if (interval == 0) return random256b();                                     // avoid division by 0
   if (interval != randomixInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder >= 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then
      randomixLimit = uint32_t(((uint64_t)1 << 32) / interval) * interval - 1;
      randomixInterval = interval;
   }
   Vec8ui prods1;                                                              // Random number * interval
   Vec8ib reject1 = true;                                                      // Initial values in prods1 are invalid
   static const Vec8ib select_odd(false,true,false,true,false,true,false,true);// Odd elements true
   bool incomplete;                                                            // when not all elements are accepted
   int n = 0;                                                                  // Counter for rejection loop

   // Rejection loop
   do {
       Vec8ui ranbits   = random256b();                                        // Random bits
       Vec8ui prod_even = mulExtended(ranbits, Vec8ui(interval));              // 64-bit product of even numbered elements
       Vec8ui prod_odd  = mulExtended(shift32Down(ranbits), Vec8ui(interval)); // 64-bit product of odd  numbered elements
       Vec8ui prods2    = select(select_odd, prod_odd, shift32Down(prod_even));// combine high parts of products
       Vec8ui remainders= select(select_odd, shift32Up(prod_odd), prod_even);  // combine low  parts of products
       Vec8ib reject2   = remainders > randomixLimit;                          // Reject these elements
       Vec8ib transfer  = andnot(reject1, reject2);                            // Transfer elements that are missing in prods1 and accepted in prods2
       prods1  = select(transfer, prods2, prods1);
       reject1 = reject1 & reject2;                                            // Indicates elements still missing
       incomplete = horizontal_or(reject1);                                    // True if any elements missing
       n++;
       if (n > 1 && incomplete) {
           // Still some elements missing after 2 or more iterations. Look for accepted elements that can be moved
           reject2 |= transfer;                                                // These elements have already been used
           int i, j = 0;
           for (i = 0; i < 8; i++) {                                           // Loop through missing elements in prods1
               if (reject1[i]) {
                   while (j < 8) {                                             // Search for accepted unused elements in prods2
                       if (!reject2[j]) {
                           prods1.insert(i, prods2[j]);                        // Move element from prods2[j] to prods1[i]
                           reject1.insert(i, false);                           // Remember this element is now accepted
                           j++;
                           break;
                       }
                       j++;
                   }
               }
           }
           incomplete = horizontal_or(reject1);                                // Stop if all elements have been fixed now
       }
   }
   while (incomplete);

   return Vec8i(prods1) + Vec8i(min);
}

// 8 floating point numbers with uniform distribution in the interval 0 <= x < 1.
// Extra bit inserted for resolution 2^-24
Vec8f Ranvec1::random8f() {
    Vec8ui const one = 0x3F800000;                                             // Binary representation of 1.0f
    Vec8ui ranbits = random256b();                                             // Get random bits
    Vec8ui r1 = one - ((ranbits >> 8) & 1);                                    // 1.0 if bit8 is 0, or 1.0-2^-24 = 0.99999994f if bit8 is 1
    Vec8ui r2 = (ranbits >> 9) | one;                                          // bits 9 - 31 inserted as mantissa
    return Vec8f(reinterpret_f(r2)) - Vec8f(reinterpret_f(r1));                // Get into interval 0 <= x < 1
}

// 4 doubles with uniform distribution in the interval 0 <= x < 1
// Resolution 2^-52. (One more bit is possible, but this would be overkill)
Vec4d Ranvec1::random4d() {
    Vec4uq const one = 0x3FF0000000000000;                                     // Binary representation of 1.0
    Vec4uq ranbits = Vec4uq(random256b());                                     // Get random bits
    Vec4uq r = (ranbits >> 12) | one;                                          // bits 12 - 63 inserted as mantissa
    return Vec4d(reinterpret_d(r)) - Vec4d(reinterpret_d(one));                // Get into interval 0 <= x < 1
}

#endif

/******************************************************************************
                      512 bit vector output functions for Ranvec1
******************************************************************************/

#if MAX_VECTOR_SIZE >= 512
// 16 integers in the interval min <= x <= max
// Relative error on frequencies < 2^-32   
Vec16i Ranvec1::random16i(int min, int max) {
   if (max <= min) {
      if (max == min) {
          return Vec16i(min);                                                  // Interval has only one value
      }
      else {
          return Vec16i(0x80000000);                                           // Error: interval length is negative
      }
   }
   Vec16ui ranbits   = random512b();                                           // Random bits
   uint32_t interval = (uint32_t)(max - min) + 1u;
   if (interval == 0) return ranbits;
   Vec16ui prod_even = mulExtended(ranbits, interval);                         // 64-bit product of even numbered elements
   Vec16ui prod_odd  = mulExtended(shift32Down(ranbits), interval);            // 64-bit product of odd numbered elements
   static const Vec16ib select_odd(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1);           // Odd elements true
   Vec16ui prods     = select(select_odd, prod_odd,shift32Down(prod_even));    // combine high parts of products
   return Vec16i(prods) + Vec16i(min);
}

// 16 integers in the interval min <= x <= max, with extra precision
// (Not guaranteed to give exactly the same sequence as random4ix called 4 times or random8ix called 2 times)
Vec16i Ranvec1::random16ix(int min, int max) {
   if (max <= min) {
      if (max == min) {
          return Vec16i(min);                                                  // Interval has only one value
      }
      else {
          return Vec16i(0x80000000);                                           // Error: interval length is negative
      }
   }
   uint32_t interval = (uint32_t)(max - min) + 1u;
   if (interval == 0) return random512b();                                     // avoid division by 0
   if (interval != randomixInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder >= 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then
      randomixLimit = uint32_t(((uint64_t)1 << 32) / interval) * interval - 1;
      randomixInterval = interval;
   }
   Vec16ui prods1;                                                             // Random number * interval
   Vec16ib reject1 = true;                                                     // Initial values in prods1 are invalid
   static const Vec16ib select_odd(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1);           // Odd elements true
   bool incomplete;                                                            // when not all elements are accepted
   int n = 0;                                                                  // Counter for rejection loop

   // Rejection loop
   do {
       Vec16ui ranbits   = random512b();                                       // Random bits
       Vec16ui prod_even = mulExtended(ranbits, Vec16ui(interval));            // 64-bit product of even numbered elements
       Vec16ui prod_odd  = mulExtended(shift32Down(ranbits),Vec16ui(interval));// 64-bit product of odd  numbered elements
       Vec16ui prods2    = select(select_odd, prod_odd,shift32Down(prod_even));// combine high parts of products
       Vec16ui remainders= select(select_odd, shift32Up(prod_odd), prod_even); // combine low  parts of products
       Vec16ib reject2   = remainders > randomixLimit;                         // Reject these elements
       Vec16ib transfer  = andnot(reject1, reject2);                           // Transfer elements that are missing in prods1 and accepted in prods2
       prods1  = select(transfer, prods2, prods1);
       reject1 = reject1 & reject2;                                            // Indicates elements still missing
       incomplete = horizontal_or(reject1);                                    // True if any elements missing
       n++;
       if (n > 1 && incomplete) {
           // Still some elements missing after 2 or more iterations. Look for accepted elements that can be moved.
           reject2 |= transfer;                                                // These elements have already been used
#if INSTRSET >=  9  // AVX512
           /*
           __m512i collect = _mm512_maskz_compress_epi32(~reject2, prods2);    // Collect accepted elements from prods2
           prods1 = _mm512_mask_expand_epi32(prods1, reject1, collect);        // Distribute these into prods1
           __m512i cmask1 = _mm512_maskz_compress_epi32(~reject2, Vec16i(-1)); // Do the same with dummy 1's
           Vec16i  cmask2 = _mm512_maskz_expand_epi32(reject1, cmask1);        // .. to generate a mask of sucessfully transferred values
           reject1 &= cmask2 == 0;                                             // Remove these from reject1
           */
           // We can rule out prods == -1, so we use this value to indicate rejected values
           __m512i collect = _mm512_mask_compress_epi32(Vec16i(-1), ~reject2, prods2); // Collect accepted elements from prods2
           prods1 = _mm512_mask_expand_epi32(prods1, reject1, collect);        // Distribute these into prods1
           reject1 = prods1 == Vec16i(-1);                                     // Rejected elements are now -1
#else
           int i, j = 0;
           for (i = 0; i < 16; i++) {                                          // Loop through missing elements in prods1
               if (reject1[i]) {
                   while (j < 16) {                                            // Search for accepted unused elements in prods2
                       if (!reject2[j]) {
                           prods1.insert(i, prods2[j]);                        // Move element from prods2[j] to prods1[i]
                           reject1.insert(i, false);                           // Remember this element is now accepted
                           j++;
                           break;
                       }
                       j++;
                   }
               }
           }
#endif
           incomplete = horizontal_or(reject1);                                // Stop if all elements have been fixed now
       }
   }
   while (incomplete);

   return Vec16i(prods1) + Vec16i(min);
}

// 16 floating point numbers with uniform distribution in the interval 0 <= x < 1.
// Extra bit inserted for resolution 2^-24
Vec16f Ranvec1::random16f() {
    Vec16ui const one = 0x3F800000;                                            // Binary representation of 1.0f
    Vec16ui ranbits = random512b();                                            // Get random bits
    Vec16ui r1 = one - ((ranbits >> 8) & 1);                                   // 1.0 if bit8 is 0, or 1.0-2^-24 = 0.99999994f if bit8 is 1
    Vec16ui r2 = (ranbits >> 9) | one;                                         // bits 9 - 31 inserted as mantissa
    return Vec16f(reinterpret_f(r2)) - Vec16f(reinterpret_f(r1));              // Get into interval 0 <= x < 1
}

// 8 doubles with uniform distribution in the interval 0 <= x < 1
// Resolution 2^-52. (One more bit is possible, but this would be overkill)
Vec8d Ranvec1::random8d() {
    Vec8uq const one = 0x3FF0000000000000;                                     // Binary representation of 1.0
    Vec8uq ranbits = Vec8uq(random512b());                                     // Get random bits
    Vec8uq r = (ranbits >> 12) | one;                                          // bits 12 - 63 inserted as mantissa
    return Vec8d(reinterpret_d(r)) - Vec8d(reinterpret_d(one));                // Get into interval 0 <= x < 1
}
#endif
