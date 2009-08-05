/*
 * Copyright (c) 2003, 2007-8 Matteo Frigo
 * Copyright (c) 2003, 2007-8 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Mon Feb  9 19:52:58 EST 2009 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_twiddle_c -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -twiddle-log3 -precompute-twiddles -no-generate-bytw -n 4 -name t3fv_4 -include t3f.h */

/*
 * This function contains 12 FP additions, 10 FP multiplications,
 * (or, 10 additions, 8 multiplications, 2 fused multiply/add),
 * 16 stack variables, 0 constants, and 8 memory accesses
 */
#include "t3f.h"

static void t3fv_4(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     INT m;
     R *x;
     x = ri;
     for (m = mb, W = W + (mb * ((TWVL / VL) * 4)); m < me; m = m + VL, x = x + (VL * ms), W = W + (TWVL * 4), MAKE_VOLATILE_STRIDE(rs)) {
	  V T2, T3, T1, Ta, T5, T8;
	  T2 = LDW(&(W[0]));
	  T3 = LDW(&(W[TWVL * 2]));
	  T1 = LD(&(x[0]), ms, &(x[0]));
	  Ta = LD(&(x[WS(rs, 3)]), ms, &(x[WS(rs, 1)]));
	  T5 = LD(&(x[WS(rs, 2)]), ms, &(x[0]));
	  T8 = LD(&(x[WS(rs, 1)]), ms, &(x[WS(rs, 1)]));
	  {
	       V T4, Tb, T9, T6;
	       T4 = VZMULJ(T2, T3);
	       Tb = VZMULJ(T3, Ta);
	       T9 = VZMULJ(T2, T8);
	       T6 = VZMULJ(T4, T5);
	       {
		    V Tc, Te, T7, Td;
		    Tc = VSUB(T9, Tb);
		    Te = VADD(T9, Tb);
		    T7 = VSUB(T1, T6);
		    Td = VADD(T1, T6);
		    ST(&(x[0]), VADD(Td, Te), ms, &(x[0]));
		    ST(&(x[WS(rs, 2)]), VSUB(Td, Te), ms, &(x[0]));
		    ST(&(x[WS(rs, 3)]), VFMAI(Tc, T7), ms, &(x[WS(rs, 1)]));
		    ST(&(x[WS(rs, 1)]), VFNMSI(Tc, T7), ms, &(x[WS(rs, 1)]));
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     VTW(0, 1),
     VTW(0, 3),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 4, "t3fv_4", twinstr, &GENUS, {10, 8, 2, 0}, 0, 0, 0 };

void X(codelet_t3fv_4) (planner *p) {
     X(kdft_dit_register) (p, t3fv_4, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_twiddle_c -simd -compact -variables 4 -pipeline-latency 8 -twiddle-log3 -precompute-twiddles -no-generate-bytw -n 4 -name t3fv_4 -include t3f.h */

/*
 * This function contains 12 FP additions, 8 FP multiplications,
 * (or, 12 additions, 8 multiplications, 0 fused multiply/add),
 * 16 stack variables, 0 constants, and 8 memory accesses
 */
#include "t3f.h"

static void t3fv_4(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     INT m;
     R *x;
     x = ri;
     for (m = mb, W = W + (mb * ((TWVL / VL) * 4)); m < me; m = m + VL, x = x + (VL * ms), W = W + (TWVL * 4), MAKE_VOLATILE_STRIDE(rs)) {
	  V T2, T3, T4;
	  T2 = LDW(&(W[0]));
	  T3 = LDW(&(W[TWVL * 2]));
	  T4 = VZMULJ(T2, T3);
	  {
	       V T1, Tb, T6, T9, Ta, T5, T8;
	       T1 = LD(&(x[0]), ms, &(x[0]));
	       Ta = LD(&(x[WS(rs, 3)]), ms, &(x[WS(rs, 1)]));
	       Tb = VZMULJ(T3, Ta);
	       T5 = LD(&(x[WS(rs, 2)]), ms, &(x[0]));
	       T6 = VZMULJ(T4, T5);
	       T8 = LD(&(x[WS(rs, 1)]), ms, &(x[WS(rs, 1)]));
	       T9 = VZMULJ(T2, T8);
	       {
		    V T7, Tc, Td, Te;
		    T7 = VSUB(T1, T6);
		    Tc = VBYI(VSUB(T9, Tb));
		    ST(&(x[WS(rs, 1)]), VSUB(T7, Tc), ms, &(x[WS(rs, 1)]));
		    ST(&(x[WS(rs, 3)]), VADD(T7, Tc), ms, &(x[WS(rs, 1)]));
		    Td = VADD(T1, T6);
		    Te = VADD(T9, Tb);
		    ST(&(x[WS(rs, 2)]), VSUB(Td, Te), ms, &(x[0]));
		    ST(&(x[0]), VADD(Td, Te), ms, &(x[0]));
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     VTW(0, 1),
     VTW(0, 3),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 4, "t3fv_4", twinstr, &GENUS, {12, 8, 0, 0}, 0, 0, 0 };

void X(codelet_t3fv_4) (planner *p) {
     X(kdft_dit_register) (p, t3fv_4, &desc);
}
#endif				/* HAVE_FMA */
