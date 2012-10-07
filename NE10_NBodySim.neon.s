@
@  Copyright 2011-12 ARM Limited
@
@  Licensed under the Apache License, Version 2.0 (the "License");
@  you may not use this file except in compliance with the License.
@  You may obtain a copy of the License at
@
@      http://www.apache.org/licenses/LICENSE-2.0
@
@  Unless required by applicable law or agreed to in writing, software
@  distributed under the License is distributed on an "AS IS" BASIS,
@  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
@  See the License for the specific language governing permissions and
@  limitations under the License.
@


        .text
        .syntax   unified

.include "NE10header.s"




        .align   4
        .global   NBodySim_neon
        .thumb
        .thumb_func

NBodySim_neon:
        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        @
        @ arm_result_t NBodySim_neon(arm_float_t * x,
        @                 arm_float_t * y,
        @                 arm_float_t * z,
		@                 arm_float_t xi,
		@                 arm_float_t yi,
		@                 arm_float_t zi,
		@                 arm_float_t * pax,
		@                 arm_float_t * pay,
		@                 arm_float_t * paz,
		@                 arm_float_t * m,
		@                 arm_float_t eps)
        @
        @  r0: *x & current x entry's address
        @  r1: *y & current y entry's address
        @  r2: *z & current z entry's address
		@  r3: xi
		@  r4: yi
		@  r5: zi
		@  r6: *pax
		@  r7: *pay
		@  r8: *paz
        @  r9: *m
		@  r10: eps
        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

@		push {r4-r10}
		stmdb sp!, {r4-r10}
		ldr				r4, [sp , #(4*7)]
		ldr				r5, [sp , #(4*8)]
		ldr				r6, [sp , #(4*9)]
		ldr				r7, [sp , #(4*10)]
		ldr				r8, [sp , #(4*11)]
		ldr				r9, [sp , #(4*12)]
		ldr				r10, [sp , #(4*13)]
		vdup.f32		q13, r3				@ lanes of xi
		vdup.f32		q14, r4				@ lanes of yi
		vdup.f32		q15, r5				@ lanes of zi
		add             r3, r0, #4000		@ hardcoded loop size: N=1000
		vld1.f32		d6[0], [r6]			@ get ax (put it in q3)
		vld1.f32		d6[1], [r7]			@ get ay (put it in q3)
		vld1.f32		d7[0], [r8]			@ get az (put it in q3)
		
.L_jloop_float:
		vldmia          r0!, {d0-d1}		@ load x[j] to q0
        vldmia          r1!, {d2-d3}		@ load y[j] to q1
		vldmia          r2!, {d4-d5}		@ load z[j] to q2
		vldmia          r9!, {d8-d9}		@ load m[j] to q4
		
		cmp r0, r3							@ get comparison ready for branch
		
        @ calculations in jloop
		vsub.f32		q0, q0, q13			@ dx=x[j]-x[i]
		vsub.f32		q1, q1, q14			@ dy=y[j]-y[i]
		vsub.f32		q2, q2, q15			@ dz=z[j]-z[i]
		vdup.f32		q5, r10				@ lanes of eps		
		vmla.f32		q5, q0, q0			@ eps + (dx*dx)
		vmla.f32		q5, q1, q1			@ $ + (dy*dy)
		vmla.f32		q5, q2, q2			@ $ + (dz*dz)
		vrsqrte.f32		q5, q5				@ 1/sqrt($)
		vmul.f32		q6, q5, q5			@ invr_sq = invr*invr
		vmul.f32		q6, q6, q5			@ invr3 = invr_sq*invr
		vmul.f32		q7, q4, q6			@ f = m[j]*invr3
		vmul.f32		q0, q7, q0			@ ax[] = f*dx
		vmul.f32		q1, q7, q1			@ ay[] = f*dy
		vmul.f32		q2, q7, q0			@ az[] = f*dx
@		vmov.i32		q5, #0				@ clear q5 to be used in pairwise add for az
@		vpadd.f32		q0, q0, q1			@ q0[0]=ax0+ax1, q0[1]=ax2+ax3, q0[2]=ay0+ay1, q0[3]=ay2+ay3
@		vpadd.f32		q2, q2, q5			@ q2[0]=az0+az1, q2[1]=az2+az3, q2[2,3]=0
@		vpadd.f32		q0, q0, q2			@ q0[0]=ax0+ax1+ax2+ax3, q0[1]=ay0+ay1+ay2+ay3, q0[2]=az0+az1+az2+az3, q0[3]=0

		vpadd.f32		d0, d0, d1			@ q0[0]=ax0+ax1, q0[1]=ax2+ax3 
		vpadd.f32		d1, d2, d3			@ q0[2]=ay0+ay1, q0[3]=ay2+ay3
		vpadd.f32		d2, d4, d5			@ q2[0]=az0+az1, q2[1]=az2+az3, q2[2,3]=0
		vpadd.f32		d0, d0, d1			@ q0[0]=ax0+ax1+ax2+ax3, q0[1]=ay0+ay1+ay2+ay3
		vadd.f32		s2, s4, s5			@ q0[0]=ax0+ax1+ax2+ax3, q0[1]=ay0+ay1+ay2+ay3, q0[2]=az0+az1+az2+az3, q0[3]=0
		vadd.f32		q3, q3, q0			@ ax+=, ay+=, az+=

		bne             .L_jloop_float          @ loop if r0 != r0+4000


.L_return_float:
		vst1.f32		d6[0], [r6]			@ store ax
		vst1.f32		d6[1], [r7]			@ store ay
		vst1.f32		d7[0], [r8]			@ store az
        @ return
@        pop {r4-r10}
		ldmia sp!, {r4-r10}
        mov               r0, #0
        bx                lr

