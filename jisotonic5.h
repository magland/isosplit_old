/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef jisotonic5_h
#define jisotonic5_h

void jisotonic5(int N, float* BB, float* MSE, float* AA, float* WW);
void jisotonic5_updown(int N, float* out, float* in, float* weights);
void jisotonic5_downup(int N, float* out, float* in, float* weights);
void jisotonic5_sort(int N, float* out, const float* in);

#endif
