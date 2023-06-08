#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <algorithm>
#include <cmath>

#include "estimate_error.hxx"

namespace AsterX {
using namespace std;

enum class regrid_t {
  first_deriv,
  second_deriv,
  second_deriv_norm,
  first_grad
};
regrid_t regridMethod;
vector<int> regridVarsI;
vector<array<int, Loop::dim> > indextypeRegridVars;

extern "C" void AsterX_EstimateError_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  vector<string> regridGroups;
  istringstream groupStream(regrid_groups);
  read_stream(regridGroups, groupStream);
  regridVarsI.clear();
  indextypeRegridVars.clear();
  for (size_t g = 0; g < regridGroups.size(); g++) {
    const int firstVarI = CCTK_FirstVarIndex(regridGroups[g].c_str());
    const int numVars = CCTK_NumVarsInGroup(regridGroups[g].c_str());
    const int gi = CCTK_GroupIndex(regridGroups[g].c_str());
    assert(firstVarI >= 0);
    assert(numVars >= 0);
    for (int i = 0; i < numVars; i++) {
      regridVarsI.push_back(firstVarI + i);
      indextypeRegridVars.push_back(get_group_indextype(gi));
    }
  }

  /* check refinement method */
  if (CCTK_EQUALS(regrid_method, "first derivative"))
    regridMethod = regrid_t::first_deriv;
  else if (CCTK_EQUALS(regrid_method, "second derivative"))
    regridMethod = regrid_t::second_deriv;
  else if (CCTK_EQUALS(regrid_method, "second derivative norm"))
    regridMethod = regrid_t::second_deriv_norm;
  else if (CCTK_EQUALS(regrid_method, "first gradient"))
    regridMethod = regrid_t::first_grad;
  else
    CCTK_ERROR("Unknown value for parameter \"regrid_method\"");

  /* print infos about regird error calculation */
  CCTK_VINFO("Regrid error will be calculated based on %ld vars (method %d):",
             regridVarsI.size(), int(regridMethod));
  ostringstream buf;
  for (size_t i = 0; i < regridVarsI.size(); i++) {
    buf << "  " << CCTK_VarName(regridVarsI[i]);
    CCTK_INFO(buf.str().c_str());
    buf.str("");
    buf.clear();
  }
}

extern "C" void AsterX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  Loop::GF3D2<const CCTK_REAL> *regridVars =
      (Loop::GF3D2<const CCTK_REAL> *)amrex::The_Arena()->alloc(
          regridVarsI.size() * sizeof(Loop::GF3D2<const CCTK_REAL>));

  const int tl = 0;
  const int maxNregridVars = 10;
  // array<Loop::GF3D2<const CCTK_REAL>, maxNregridVars> regridVars = {
  //    rho, rho, rho, rho, rho, rho, rho, rho, rho, rho};
  const int NregridVars = regridVarsI.size();
  assert(NregridVars < maxNregridVars);
  for (int i = 0; i < NregridVars; i++) {
    const Loop::GF3D2layout layout(cctkGH, indextypeRegridVars[i]);
    const Loop::GF3D2<const CCTK_REAL> gf(
        layout,
        static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, regridVarsI[i])));
    regridVars[i] = gf;
  }

  auto regridMethod_local = regridMethod;
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        switch (regridMethod_local) {

        case regrid_t::first_deriv: {
          regrid_error(p.I) = 0.0;
          for (int i = 0; i < NregridVars; i++) {
            regrid_error(p.I) =
                max({regrid_error(p.I), calc_deriv_1st(regridVars[i], p)});
          }
          break;
        }

        case regrid_t::second_deriv: {
          regrid_error(p.I) = 0.0;
          for (int i = 0; i < NregridVars; i++) {
            regrid_error(p.I) =
                max({regrid_error(p.I), calc_deriv_2nd(regridVars[i], p)});
          }
          break;
        }

        case regrid_t::second_deriv_norm: {
          regrid_error(p.I) = 0.0;
          for (int i = 0; i < NregridVars; i++) {
            regrid_error(p.I) =
                max({regrid_error(p.I),
                     calc_deriv_2nd_norm(regridVars[i], p, epsilon_err)});
          }
          break;
        }

        case regrid_t::first_grad: {
          regrid_error(p.I) = 0.0;
          for (int i = 0; i < NregridVars; i++) {
            regrid_error(p.I) =
                max({regrid_error(p.I), calc_grad_1st(regridVars[i], p)});
          }
          break;
        }

        default:
          assert(0);
        }
      });

  amrex::Gpu::Device::streamSynchronize();
  amrex::The_Arena()->free(regridVars);
}

} // namespace AsterX
