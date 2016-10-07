#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <stdint.h>
#include <iostream>

NTL_CLIENT


typedef uint4 uint128;
typedef uint3 uint96;
typedef uint2 uint64;
typedef uint uint32;

template <typename T>
static inline ZZ to_ZZ(const T & a){
	return ZZFromBytes((const unsigned char *)&a, sizeof(T));
}
template <typename T>
static inline ZZ_p to_ZZ_p(const T & a){
	return to_ZZ_p(ZZFromBytes((const unsigned char *)&a, sizeof(T)));
}

template <typename T>
static inline T to_uint(const ZZ & a)
{
	T res;
	BytesFromZZ((unsigned char *)&res, a, sizeof(T));
	return res;
}

static inline ZZ to_ZZ(const uint96 & a) { return to_ZZ<uint96>(a); }
static inline ZZ_p to_ZZ_p(const uint96 & a) { return to_ZZ_p<uint96>(a); }

static inline uint96 to_uint96(const ZZ & a) { return to_uint<uint96>(a); }
static inline ZZ to_ZZ(const uint128 & a) { return to_ZZ<uint128>(a); }
static inline ZZ_p to_ZZ_p(const uint128 & a) { return to_ZZ_p<uint128>(a); }
static inline uint128 to_uint128(const ZZ & a) { return to_uint<uint128>(a); }


/*__device__ __forceinline__ void get_q(uint96 & a, const uint96 b, const uint96 c)
{
    uint temp0 = 0, temp1 = 0;

   asm("mul.hi.u32	%0, %4, %8;\n\t"		//(O1 = r0 * s0.hi)
"mad.lo.cc.u32	%0, %6, %8, %0;\n\t"		//(O1 += r1 * s0.lo)
"madc.hi.u32	%1, %6, %8, 0;\n\t"		//(O2 = r1 * s0.hi + c)
"mad.lo.cc.u32	%0, %5, %9, %0;\n\t"		//(O1 += r0 * s1.lo)
"madc.hi.cc.u32	%1, %5, %9, %1;\n\t"		//(O2 += r0 * s1.hi + c)
"madc.lo.u32	%2, %7, %9, 0;\n\t"		//(O3 = r2 * s1.lo + c  )


"mad.lo.cc.u32	%1, %7, %8, %1;\n\t"		//(O2 += r2 * s0.lo)
"madc.hi.cc.u32	%2, %7, %8, %2;\n\t"		//(O3 += r2 * s0.hi + c)
"madc.hi.u32	%3, %7, %9, 0;\n\t"		//(O4 = r2 * s1.hi + c   )
"mad.lo.cc.u32	%1, %6, %9, %1;\n\t"		//(O2 += r1 * s1.lo)
"madc.hi.cc.u32	%2, %6, %9, %2;\n\t"		//(O3 += r1 * s1.hi + c)
"addc.cc.u32	%3, %7, %3;\n\t"		//(O4 += r2 + c  )
"addc.u32	%4, 0, 0;\n\t"		//(O5 = +c)
"add.cc.u32	%1, %5, %1;\n\t"		//(O2 += r0)


"addc.cc.u32	%2, %6, %2;\n\t"		//(O3 += r1 + c)

"addc.cc.u32	%3, 0, %3;\n\t"		//(O4 += +c
"addc.u32	%4, 0, %4;\n\t"		//(O5 + =+c)

	 : "=r"(temp0), "=r"(temp1), "=r"(a.x), "=r"(a.y), "=r"(a.z) 
	 : "r"(b.x), "r"(b.y), "r"(b.z),"r"(c.x), "r"(c.y));
}*/

/*__device__ __forceinline__ void get_q(uint96 & a, const uint128 b, const uint96 c)
{

    asm("mul.hi.u32      %1,%3,%5   ;\n\t"  // (a0 * b0)_hi 	*11*

	"mad.lo.cc.u32   %1,%3,%6,%1;\n\t"  // (a0 * b1)_lo 	*11*
	"madc.hi.u32     %2,%3,%6, 0;\n\t"  // (a0 * b1)_hi 	!22!

	"mad.lo.cc.u32   %1,%4,%5,%1;\n\t"  // (a1 * b0)_lo 	*11*
	"madc.hi.cc.u32  %2,%4,%5,%2;\n\t"  // (a1 * b0)_hi 	!22!
	"madc.hi.u32     %0,%3,%7, 0;\n\t"  // (a0 * b2)_hi     <--0

	"mad.lo.cc.u32   %2,%3,%7,%2;\n\t"  // (a0 * b2)_lo 	!22!
	"madc.hi.cc.u32  %0,%4,%6,%0;\n\t"  // (a1 * b1)_hi
	"madc.hi.u32     %1,%4,%7, 0;\n\t"  // (a1 * b2)_hi

	"madc.lo.cc.u32  %2,%4,%6,%5;\n\t"  //			!22!
	"madc.lo.cc.u32  %0,%4,%7,%0;\n\t"  // (a1 * b2)_lo
	"addc.cc.u32     %1,%1,%7   ;\n\t"  // (a2 * b2)_lo
	"addc.u32        %2, 0, 0   ;\n\t"  // propagate carry  <--2

	"add.cc.u32      %0,%0,%6   ;\n\t"  // (a2 * b1)_lo
	"addc.cc.u32     %1,%1, 0   ;\n\t"  // propagate carry
	"addc.u32        %2,%2, 0   ;\n\t"  // propagate carry
	 : "=r"(a.x), "=r"(a.y), "=r"(a.z) 
	 : "r"(c.x), "r"(c.y), "r"(b.x), "r"(b.y), "r"(b.z));
}*/


/*__device__ __forceinline__ void get_q(uint96 & a, const uint128 b, const uint96 c)
{
    uint temp0 = 0, temp1 = 0, temp2=0,temp3=0;

   asm("mul.hi.u32     %0, %7, %11;\n\t"               //(O0 = r0 * s0.hi)

"mad.lo.cc.u32  %0, %8, %11, %0;\n\t"           //(O0 += r1 * s0.lo)
"madc.hi.u32    %1, %8, %11, 0;\n\t"            //(O1 = r1 * s0.hi + c)

"mad.lo.cc.u32  %0, %7, %12, %0;\n\t"           //(O0 += r0 * s1.lo)
"madc.hi.cc.u32 %1, %7, %12, %1;\n\t"           //(O1 += r0 * s1.hi + c)
"madc.lo.u32    %2, %10, %11, 0;\n\t"           //(O2 = r3 * s0.lo + c  )

"mad.lo.cc.u32  %1, %9, %11, %1;\n\t"           //(O1 += r2 * s0.lo)
"madc.hi.cc.u32 %2, %9, %11, %2;\n\t"           //(O2 += r2 * s0.hi + c)
"madc.hi.u32    %3, %10, %11, 0;\n\t"           //(O3 = r3 * s0.hi + c   )

"mad.lo.cc.u32  %1, %8, %12, %1;\n\t"           //(O1 += r1 * s1.lo)
"madc.hi.cc.u32 %2, %8, %12, %2;\n\t"           //(O2 += r1 * s1.hi + c)
"madc.lo.cc.u32 %3, %10, %12, %3;\n\t"          //(O3 += r3 * s1.lo + c  )
"addc.u32       %4, 0, 0;\n\t"          	//(O4 = +c)

"mad.lo.cc.u32  %1, %7, %13, %1;\n\t"           //(O1 += r0 * s2.lo)
"madc.hi.cc.u32 %2, %7, %13, %2;\n\t"           //(O2 += r0 * s2.hi + c)
"madc.hi.u32    %3, %9, %12, %3;\n\t"           //(O3 += r2 * s1.hi + c   )
"addc.u32       %4, %4, 0;\n\t"         	//(O4 += +c)

"mad.lo.cc.u32  %2, %9, %12, %2;\n\t"           //(O2 += r2 * s1.lo)
"madc.hi.cc.u32 %3, %8, %13, %3;\n\t"           //(O3 += r1 * s2.hi + c)
"madc.lo.cc.u32 %4, %10, %13, %4;\n\t"          //(O4 += r3 * s2.lo + c  )
"addc.u32       %5, 0, 0;\n\t"          	//(O5 = +c)

"mad.lo.cc.u32  %2, %8, %13, %2;\n\t"           //(O2 += r1 * s2.lo)
"addc.cc.u32   %3, %7, %3;\n\t"                	//(O3 += r0 + c)
"madc.hi.u32    %4, %10, %12, %4;\n\t"          //(O4 += r3 * s1.hi + c   )
"addc.u32       %5, %5, 0;\n\t"         	//(O5 += +c)

"mad.lo.cc.u32  %3, %9, %13, %3;\n\t"           //(O3 += r2 * s2.lo)
"madc.hi.cc.u32 %4, %9, %13, %4;\n\t"           //(O4 += r2 * s2.hi + c)
"addc.cc.u32   %5, %10, %5;\n\t"               	//(O5 += r3 + c  )
"addc.u32       %6, 0, 0;\n\t"          	//(O6 = +c)

"add.cc.u32     %3, %8, %3;\n\t"                //(O3 += r1)
"addc.cc.u32    %4, %9, %4;\n\t"                //(O4 += r2 + c)
"madc.hi.cc.u32 %5, %10, %13, %5;\n\t"          //(O5 += r3 * s2.hi + c)
"addc.u32       %6, 0, %6;\n\t"         	//(O6 + =+c)
	 : "=r"(temp0), "=r"(temp1),"=r"(temp2), "=r"(temp3), "=r"(a.x), "=r"(a.y), "=r"(a.z) 
	 : "r"(b.x), "r"(b.y), "r"(b.z),"r"(b.w),"r"(c.x), "r"(c.y),"r"(c.z));
printf("\nGPU: %u\t%u\t%u\t%u\t%u\t%u\t%u\n",temp0, temp1, temp2, temp3, a.x, a.y, a.z);
}*/

__device__ __forceinline__ void get_q(uint96 & a, const uint128 b, const uint96 c)
{
    uint temp0 = 0, temp1 = 0, temp2=0,temp3=0;

   asm(
"mul.hi.u32     %0, %7, %11;\n\t"        //(O0 = r0 * s0.hi)

"mad.lo.cc.u32  %0, %8, %11, %0;\n\t"           //(O0 += r1 * s0.lo)
"madc.hi.u32    %1, %8, %11, 0;\n\t"            //(O1 = r1 * s0.hi + c)

"mad.lo.cc.u32  %0, %7, %12, %0;\n\t"           //(O0 += r0 * s1.lo)
"madc.hi.cc.u32 %1, %7, %12, %1;\n\t"           //(O1 += r0 * s1.hi + c)
"madc.lo.u32    %2, %10, %11, 0;\n\t"           //(O2 = r3 * s0.lo + c  )


"mad.lo.cc.u32  %1, %9, %11, %1;\n\t"           //(O1 += r2 * s0.lo)
"madc.hi.cc.u32 %2, %9, %11, %2;\n\t"           //(O2 += r2 * s0.hi + c)
"madc.hi.u32    %3, %10, %11, 0;\n\t"           //(O3 = r3 * s0.hi + c   )

"mad.lo.cc.u32  %1, %8, %12, %1;\n\t"           //(O1 += r1 * s1.lo)
"madc.hi.cc.u32 %2, %8, %12, %2;\n\t"           //(O2 += r1 * s1.hi + c)
"madc.lo.cc.u32 %3, %10, %12, %3;\n\t"          //(O3 += r3 * s1.lo + c  )
"addc.u32       %4, 0, 0;\n\t"          	//(O4 = +c)

"mad.lo.cc.u32  %1, %7, %13, %1;\n\t"           //(O1 += r0 * s2.lo)
"madc.hi.cc.u32 %2, %7, %13, %2;\n\t"           //(O2 += r0 * s2.hi + c)
"madc.lo.cc.u32 %3, %9, %13, %3;\n\t"           //(O3 += r2 * s2.lo + c  )
"addc.u32       %4, 0, 0;\n\t"          	//(O4 = +c)



"mad.lo.cc.u32  %2, %9, %12, %2;\n\t"           //(O2 += r2 * s1.lo)
"madc.hi.cc.u32 %3, %9, %12, %3;\n\t"           //(O3 += r2 * s1.hi + c)
"madc.hi.u32    %4, %10, %12, %4;\n\t"          //(O4 += r3 * s1.hi + c   )
"addc.u32       %5, 0, 0;\n\t"          	//(O5 = +c)

"mad.lo.cc.u32  %2, %8, %13, %2;\n\t"           //(O2 += r1 * s2.lo)
"madc.hi.cc.u32 %3, %8, %13, %3;\n\t"           //(O3 += r1 * s2.hi + c)
"madc.lo.cc.u32 %4, %10, %13, %4;\n\t"          //(O4 += r3 * s2.lo + c  )
"addc.u32       %5, 0, 0;\n\t"          	//(O5 = +c)

"add.cc.u32     %2, %7, %2;\n\t"                //(O2 += r0)
"addc.cc.u32    %3, %8, %3;\n\t"                //(O3 += r1 + c)
"madc.hi.cc.u32 %4, %9, %13, %4;\n\t"           //(O4 += r2 * s2.hi + c)
"madc.hi.u32    %5, %10, %13, %5;\n\t"          //(O5 += r3 * s2.hi + c   )
"addc.u32       %6, 0, 0;\n\t"          	//(O6 = +c)

"add.cc.u32     %4, %9, %4;\n\t"                //(O4 += r2)
"addc.cc.u32    %5, %10, %5;\n\t"               //(O5 += r3 + c)
"addc.u32       %6, 0, %6;\n\t"         	//(O6 + =+c)
	 : "=r"(temp0), "=r"(temp1),"=r"(temp2), "=r"(temp3), "=r"(a.x), "=r"(a.y), "=r"(a.z) 
	 : "r"(b.x), "r"(b.y), "r"(b.z),"r"(b.w),"r"(c.x), "r"(c.y),"r"(c.z));
printf("\nGPU: %u\t%u\t%u\t%u\t%u\t%u\t%u\n",temp0, temp1, temp2, temp3, a.x, a.y, a.z);
}


__global__ void kernel(uint96 *a, const uint128 b, const uint96 c)
{
    get_q(*a, b, c);
}

int main(){
	ZZ modulus = RandomPrime_ZZ(128);
	ZZ_p::init(modulus);
	ZZ b = RandomBnd(modulus);
	ZZ c = RandomBits_ZZ(96);
	cout << b << "<-b, c->" << (c + power2_ZZ(96)) << "\n";
	ZZ a = (b * (c + power2_ZZ(96))) ;

	uint128 b_ = to_uint128(b);
	uint96 c_ = to_uint96(c);
	cout << to_ZZ(b_) << "<-b, c->" << (to_ZZ(c_) + power2_ZZ(96));
	uint96 a_;

	uint96 * a__;
	cudaMalloc((void**) & a__, sizeof(uint96));

	kernel<<< 1, 1 >>>(a__, b_, c_);

	cudaMemcpy(&a_, a__, sizeof(uint96), cudaMemcpyDeviceToHost);
	cudaFree(a__);

	ZZ a2 = to_ZZ(a_);

	cout << "CPU: " << trunc_long(a>>32,32) << "\t" << trunc_long(a>>64,32) << "\t" << trunc_long(a>>96,32) << "\t" << trunc_long(a>>128,32) << "\t" << trunc_long(a>>160,32) << "\t" << trunc_long(a>>192,32) << "\t" << trunc_long(a>>224,32) << "\n";
	//cout << "GPU\t" << a_.x << "\n";
}
