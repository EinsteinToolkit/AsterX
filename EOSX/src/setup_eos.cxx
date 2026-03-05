#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>

#include <setup_eos.hxx>

#include <mpi.h>
#include <hdf5.h>

#include "eos_3p_tabulated3d/eos_readtable_scollapse.hxx"
#include "eos_3p_tabulated3d/eos_readtable_compose.hxx"

namespace EOSX {

using namespace amrex;

enum class eos_1param { Polytropic, PWPolytropic };
enum class eos_3param { IdealGas, Hybrid, Tabulated };

// initial data EOS
eos_1p_polytropic *global_eos_1p_poly = nullptr;
eos_1p_piecewise_polytropic *global_eos_1p_pwpoly = nullptr;

// evolution EOS
eos_3p_idealgas *global_eos_3p_ig = nullptr;
eos_3p_hybrid_poly *global_eos_3p_hyb_poly = nullptr;
eos_3p_hybrid_pwpoly *global_eos_3p_hyb_pwpoly = nullptr;
eos_3p_tabulated3d *global_eos_3p_tab3d = nullptr;

enum class eos_table_format { StellarCollapse = 0, Compose = 1 };

static inline eos_table_format
detect_table_format_rank0(const std::string &filename) {
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file_id >= 0);

  // StellarCollapse "signature"
  const bool has_pointsrho = (H5Lexists(file_id, "pointsrho", H5P_DEFAULT) > 0);
  const bool has_pointstemp =
      (H5Lexists(file_id, "pointstemp", H5P_DEFAULT) > 0);
  const bool has_pointsye = (H5Lexists(file_id, "pointsye", H5P_DEFAULT) > 0);

  // CompOSE "signature"
  const bool has_parameters =
      (H5Lexists(file_id, "/Parameters", H5P_DEFAULT) > 0);
  const bool has_thermo = (H5Lexists(file_id, "/Thermo_qty", H5P_DEFAULT) > 0);

  H5Fclose(file_id);

  if (has_pointsrho && has_pointstemp && has_pointsye)
    return eos_table_format::StellarCollapse;
  if (has_parameters && has_thermo)
    return eos_table_format::Compose;

  CCTK_ERROR(
      "Could not auto-detect EOS table format. "
      "Set EOSX::EOSTable_format to \"StellarCollapse\" or \"Compose\".");
  return eos_table_format::StellarCollapse;
}

CCTK_HOST void eos_readtable(const std::string &filename,
                             eos_tabulated3d_raw_table_t &tab) {
  DECLARE_CCTK_PARAMETERS;

  eos_table_format fmt;

  if (CCTK_EQUALS(EOSTable_format, "StellarCollapse")) {
    fmt = eos_table_format::StellarCollapse;
  } else if (CCTK_EQUALS(EOSTable_format, "Compose")) {
    fmt = eos_table_format::Compose;
  } else if (CCTK_EQUALS(EOSTable_format, "Auto")) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int ifmt = 0;
    if (rank == 0) {
      fmt = detect_table_format_rank0(filename);
      ifmt = (int)fmt;
    }
    MPI_Bcast(&ifmt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    fmt = (eos_table_format)ifmt;
  } else {
    CCTK_ERROR("Unknown value for parameter \"EOSTable_format\"");
    fmt = eos_table_format::StellarCollapse;
  }

  switch (fmt) {
  case eos_table_format::StellarCollapse:
    eos_readtable_scollapse(filename, tab);
    break;
  case eos_table_format::Compose:
    eos_readtable_compose(filename, tab);
    break;
  default:
    assert(false);
  }
}

extern "C" void EOSX_Setup_EOSID(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  eos_1param eos_1p_type;
  if (CCTK_EQUALS(initial_data_eos, "Polytropic")) {
    eos_1p_type = eos_1param::Polytropic;
  } else if (CCTK_EQUALS(initial_data_eos, "PWPolytropic")) {
    eos_1p_type = eos_1param::PWPolytropic;
  } else {
    CCTK_ERROR("Unknown value for parameter \"initial_data_eos\"");
  }

  switch (eos_1p_type) {
  case eos_1param::Polytropic: {
    CCTK_INFO("Setting initial data EOS to Polytropic");
    global_eos_1p_poly = (eos_1p_polytropic *)The_Managed_Arena()->alloc(
        sizeof *global_eos_1p_poly);
    assert(global_eos_1p_poly);
    new (global_eos_1p_poly) eos_1p_polytropic;
    global_eos_1p_poly->init(poly_gamma, poly_k, rho_max);
    break;
  }
  case eos_1param::PWPolytropic: {
    CCTK_INFO("Setting initial data EOS to Piecewise Polytropic");

    global_eos_1p_pwpoly =
        (eos_1p_piecewise_polytropic *)The_Managed_Arena()->alloc(
            sizeof *global_eos_1p_pwpoly);
    assert(global_eos_1p_pwpoly);
    new (global_eos_1p_pwpoly) eos_1p_piecewise_polytropic;

    EOSX::validate_pwpoly_params(pwpoly_nsegm, pwpoly_segm_bound,
                                 pwpoly_segm_gamma);
    global_eos_1p_pwpoly->init(pwpoly_rho_p0, pwpoly_nsegm, pwpoly_segm_bound,
                               pwpoly_segm_gamma, rho_max);
    break;
  }
  default:
    assert(0);
  }
}

extern "C" void EOSX_Setup_EOS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  eos_3param eos_3p_type;
  eos_3p::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
      rgye(ye_min, ye_max);

  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eos_3p_type = eos_3param::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eos_3p_type = eos_3param::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    eos_3p_type = eos_3param::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  switch (eos_3p_type) {
  case eos_3param::IdealGas: {
    CCTK_INFO("Setting evolution EOS to Ideal Gas");
    global_eos_3p_ig =
        (eos_3p_idealgas *)The_Managed_Arena()->alloc(sizeof *global_eos_3p_ig);
    assert(global_eos_3p_ig);
    new (global_eos_3p_ig) eos_3p_idealgas;
    global_eos_3p_ig->init(gl_gamma, particle_mass, rgeps, rgrho, rgye);
    break;
  }
  case eos_3param::Hybrid: {
    CCTK_INFO("Setting evolution EOS to Hybrid");

    if (global_eos_1p_pwpoly) {
      global_eos_3p_hyb_pwpoly =
          (eos_3p_hybrid_pwpoly *)The_Managed_Arena()->alloc(
              sizeof *global_eos_3p_hyb_pwpoly);
      assert(global_eos_3p_hyb_pwpoly);
      new (global_eos_3p_hyb_pwpoly) eos_3p_hybrid_pwpoly(
          global_eos_1p_pwpoly, gamma_th, rgeps, rgrho, rgye);
    } else {
      global_eos_3p_hyb_poly = (eos_3p_hybrid_poly *)The_Managed_Arena()->alloc(
          sizeof *global_eos_3p_hyb_poly);
      assert(global_eos_3p_hyb_poly);
      new (global_eos_3p_hyb_poly)
          eos_3p_hybrid_poly(global_eos_1p_poly, gamma_th, particle_mass, rgeps, rgrho, rgye);
    }

    break;
  }
  case eos_3param::Tabulated: {
    CCTK_INFO("Setting evolution EOS to Tabulated3D");
    const string eos_filename = EOSTable_filename;
    global_eos_3p_tab3d = (eos_3p_tabulated3d *)The_Managed_Arena()->alloc(
        sizeof *global_eos_3p_tab3d);
    assert(global_eos_3p_tab3d);
    new (global_eos_3p_tab3d) eos_3p_tabulated3d;
    global_eos_3p_tab3d->init(eos_filename, rgeps, rgrho, rgye);
    break;
  }
  default:
    assert(0);
  }
}

} // namespace EOSX
