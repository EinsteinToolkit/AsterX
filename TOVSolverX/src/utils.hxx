#include <loop.hxx>

namespace TOVSolverX {

using namespace Loop;
/* - utility routine
   - fills an real-array 'var' of size 'i' with value 'r' */
void TOVX_C_fill(CCTK_REAL *var, CCTK_INT i, CCTK_REAL r) {
  for (i--; i >= 0; i--)
    var[i] = r;
}

void TOVX_Copy(CCTK_INT size, CCTK_REAL *var_p, CCTK_REAL *var) {
#pragma omp parallel for
  for (int i = 0; i < size; i++)
    var_p[i] = var[i];
}

template <typename T>
CCTK_DEVICE CCTK_HOST T calc_avg_c2v(const GF3D2<T> &gf, const PointDesc &p) {
  //  constexpr auto DI = PointDesc::DI;
  CCTK_REAL gf_avg = 0.0;
  for (int dk = 0; dk < 2; ++dk) {
    for (int dj = 0; dj < 2; ++dj) {
      for (int di = 0; di < 2; ++di) {
        gf_avg += gf(p.I - p.DI[0] * di - p.DI[1] * dj - p.DI[2] * dk);
      }
    }
  }
  return gf_avg / 8.0;
}

} // namespace TOVSolverX
