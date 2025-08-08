#ifndef EOS_3P_TABULATED3D_HXX
#define EOS_3P_TABULATED3D_HXX

#include <cmath>
#include <cassert>
#include <limits>
#include <string>
#include <array>
#include <mpi.h>
#include <hdf5.h>

#include "../eos_3p.hxx"
#include "../utils/eos_brent.hxx" // zero_brent
#include "../utils/eos_linear_interp_ND.hxx"

#define NTABLES 19

namespace EOSX {
using namespace std;
using namespace eos_constants;

class eos_3p_tabulated3d : public eos_3p {
public:
  // must match order of HDF5 datasets below
  enum EV {
    PRESS = 0,   // "logpress"
    EPS = 1,     // "logenergy"
    S = 2,       // "entropy"
    MUNU = 3,    // "munu"
    CS2 = 4,     // "cs2"
    DEDT = 5,    // "dedt"
    DPDRHOE = 6, // "dpdrhoe"
    DPDERHO = 7, // "dpderho"
    MUHAT = 8,   // "muhat"
    MU_E = 9,    // "mu_e"
    MU_P = 10,   // "mu_p"
    MU_N = 11,   // "mu_n"
    XA = 12,     // "Xa"
    XH = 13,     // "Xh"
    XN = 14,     // "Xn"
    XP = 15,     // "Xp"
    ABAR = 16,   // "Abar"
    ZBAR = 17,   // "Zbar"
    GAMMA = 18,  // "Gamma"
    NUM_VARS
  };

  CCTK_REAL gamma; // table Γ
  CCTK_REAL *energy_shift;
  range rgeps;
  linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> *interptable;

  CCTK_HOST void init(const std::string &filename, range &rgeps_out,
                      const range &, const range &) {
    CCTK_VINFO("Reading Stellar Collapse EOS table '%s'", filename.c_str());

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    assert(fapl_id >= 0);

    hid_t file_id = -1;
#ifdef H5_HAVE_PARALLEL
    CHECK_ERROR(H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL));
    CHECK_ERROR(H5Pset_all_coll_metadata_ops(fapl_id, true));
    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl_id);
    assert(file_id >= 0);
#else
    if (rank == 0) {
      file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      assert(file_id >= 0);
    }
#endif

    int nrho, ntemp, nye;
    if (rank == 0) {
      get_hdf5_int_dset(file_id, "pointsrho", 1, &nrho);
      get_hdf5_int_dset(file_id, "pointstemp", 1, &ntemp);
      get_hdf5_int_dset(file_id, "pointsye", 1, &nye);
    }
    MPI_Bcast(&nrho, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ntemp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nye, 1, MPI_INT, 0, MPI_COMM_WORLD);

    const int npoints = nrho * ntemp * nye;

    CCTK_REAL *logrho = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
        nrho * sizeof(CCTK_REAL));
    CCTK_REAL *logtemp = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
        ntemp * sizeof(CCTK_REAL));
    CCTK_REAL *yes =
        (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(nye * sizeof(CCTK_REAL));
    CCTK_REAL *alltables_tmp = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
        npoints * NTABLES * sizeof(CCTK_REAL));
    CCTK_REAL *alltables = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
        npoints * NTABLES * sizeof(CCTK_REAL));
    energy_shift =
        (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(sizeof(CCTK_REAL));

    static const char *dnames[NTABLES] = {
        "logpress", "logenergy", "entropy", "munu", "cs2",  "dedt", "dpdrhoe",
        "dpderho",  "muhat",     "mu_e",    "mu_p", "mu_n", "Xa",   "Xh",
        "Xn",       "Xp",        "Abar",    "Zbar", "gamma"};

#ifdef H5_HAVE_PARALLEL
    get_hdf5_real_dset(file_id, "logrho", nrho, logrho);
    get_hdf5_real_dset(file_id, "logtemp", ntemp, logtemp);
    get_hdf5_real_dset(file_id, "ye", nye, yes);
    get_hdf5_real_dset(file_id, "energy_shift", 1, energy_shift);
    for (int iv = 0; iv < NTABLES; iv++) {
      get_hdf5_real_dset(file_id, dnames[iv], npoints,
                         &alltables_tmp[iv * npoints]);
    }

    // change ordering of alltables array so that
    // the table kind is the fastest changing index
    for (int iv = 0; iv < NTABLES; iv++)
      for (int k = 0; k < nye; k++)
        for (int j = 0; j < ntemp; j++)
          for (int i = 0; i < nrho; i++) {
            int indold = i + nrho * (j + ntemp * (k + nye * iv));
            int indnew = iv + NTABLES * (i + nrho * (j + ntemp * k));
            alltables[indnew] = alltables_tmp[indold];
          }
    amrex::The_Managed_Arena()->free(alltables_tmp);

    CHECK_ERROR(H5Fclose(file_id));
    CHECK_ERROR(H5Pclose(fapl_id));
#else
    if (rank == 0) {
      get_hdf5_real_dset(file_id, "logrho", nrho, logrho);
      get_hdf5_real_dset(file_id, "logtemp", ntemp, logtemp);
      get_hdf5_real_dset(file_id, "ye", nye, yes);
      get_hdf5_real_dset(file_id, "energy_shift", 1, energy_shift);
      for (int iv = 0; iv < NTABLES; iv++) {
        get_hdf5_real_dset(file_id, dnames[iv], npoints,
                           &alltables_tmp[iv * npoints]);
      }
      CHECK_ERROR(H5Fclose(file_id));

      // change ordering of alltables array so that
      // the table kind is the fastest changing index
      for (int iv = 0; iv < NTABLES; iv++)
        for (int k = 0; k < nye; k++)
          for (int j = 0; j < ntemp; j++)
            for (int i = 0; i < nrho; i++) {
              int indold = i + nrho * (j + ntemp * (k + nye * iv));
              int indnew = iv + NTABLES * (i + nrho * (j + ntemp * k));
              alltables[indnew] = alltables_tmp[indold];
            }
    }
    amrex::The_Managed_Arena()->free(alltables_tmp);

    MPI_Bcast(logrho, nrho, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(logtemp, ntemp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(yes, nye, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(energy_shift, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(alltables, npoints * NTABLES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    *energy_shift *= EPSGF;
    const CCTK_REAL ln10 = log(10.0);
    const CCTK_REAL inv_time2 = 1 / (TIMEGF * TIMEGF);
    for (int i = 0; i < nrho; i++)
      logrho[i] = logrho[i] * ln10 + log(RHOGF);
    for (int i = 0; i < ntemp; i++)
      logtemp[i] *= ln10;

    for (size_t idx = 0; idx < npoints; idx++) {
      size_t b = idx * NTABLES;
      alltables[b + PRESS] = alltables[b + PRESS] * ln10 + log(PRESSGF);
      alltables[b + EPS] = alltables[b + EPS] * ln10 + log(EPSGF);
      alltables[b + CS2] = alltables[b + CS2] * LENGTHGF * LENGTHGF * inv_time2;
      alltables[b + DEDT] = alltables[b + DEDT] * EPSGF;
      alltables[b + DPDRHOE] = alltables[b + DPDRHOE] * PRESSGF / RHOGF;
      alltables[b + DPDERHO] = alltables[b + DPDERHO] * PRESSGF / EPSGF;
    }

    std::array<size_t, 3> dims{size_t(nrho), size_t(ntemp), size_t(nye)};
    interptable = (decltype(interptable))amrex::The_Managed_Arena()->alloc(
        sizeof(*interptable));
    new (interptable) linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES>(
        alltables, dims, logrho, logtemp, yes);

    set_range_rho(
        range{exp(interptable->xmin<0>()), exp(interptable->xmax<0>())});
    set_range_temp(
        range{exp(interptable->xmin<1>()), exp(interptable->xmax<1>())});
    set_range_ye(range{interptable->xmin<2>(), interptable->xmax<2>()});

    rgeps = compute_eps_range_full_table();
    rgeps_out = rgeps;
  }

  // Device‐callable routines

  //// Table inversions
  template<size_t var>
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  logtemp_from_valid_rho_var_ye(const CCTK_REAL rho, CCTK_REAL &invar,
                                const CCTK_REAL ye) const {
    // bound inputs
    CCTK_REAL r = std::fmin(std::fmax(rho, rgrho.min), rgrho.max);
    CCTK_REAL lrho = std::log(r);

    // table‐edge clamp
    auto vmin =
        interptable->interpolate<var>(lrho, interptable->xmin<1>(), ye)[0];
    auto vmax =
        interptable->interpolate<var>(lrho, interptable->xmax<1>(), ye)[0];
    if (invar <= vmin) {
      // eps = exp(vmin) - *energy_shift;
      invar = vmin;
      return interptable->xmin<1>();
    }
    if (invar >= vmax) {
      // eps = exp(vmax) - *energy_shift;
      invar = vmax;
      return interptable->xmax<1>();
    }

    // root‐find for logtemp
    auto func = [&](CCTK_REAL &lt) {
      CCTK_REAL val = interptable->interpolate<var>(lrho, lt, ye)[0];
      return invar - val;
    };
    return zero_brent(interptable->xmin<1>(), interptable->xmax<1>(), 1.e-14,
                      func);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  logtemp_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                                const CCTK_REAL ye) const {
    // bound inputs
    eps = std::fmax(eps, rgeps.min);
    CCTK_REAL leps = std::log(eps + *energy_shift);
    CCTK_REAL lt = logtemp_from_valid_rho_var_ye<EV::EPS>(rho, leps, ye);
    eps = exp(leps) - *energy_shift;
    return lt;
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  temp_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                                const CCTK_REAL ye) const {
    CCTK_REAL lt = logtemp_from_valid_rho_eps_ye(rho, eps, ye);
    return exp(lt);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  temp_from_valid_rho_entropy_ye(const CCTK_REAL rho, CCTK_REAL &ent,
                                const CCTK_REAL ye) const {
    CCTK_REAL lt = logtemp_from_valid_rho_var_ye<EV::S>(rho, ent, ye);
    return exp(lt);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  press_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                               const CCTK_REAL ye) const {
    // bound
    CCTK_REAL r = std::fmin(std::fmax(rho, rgrho.min), rgrho.max);
    CCTK_REAL t = std::fmin(std::fmax(temp, rgtemp.min), rgtemp.max);
    CCTK_REAL lr = std::log(rho), lt = std::log(temp);
    CCTK_REAL v = interptable->interpolate<EV::PRESS>(lr, lt, ye)[0];
    return exp(v);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  press_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                              const CCTK_REAL ye) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL lt = logtemp_from_valid_rho_eps_ye(rho, eps, ye);
    CCTK_REAL v = interptable->interpolate<EV::PRESS>(lr, lt, ye)[0];
    return exp(v);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  eps_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                             const CCTK_REAL ye) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL lt = std::log(std::fmin(std::fmax(temp, rgtemp.min), rgtemp.max));
    CCTK_REAL v = interptable->interpolate<EV::EPS>(lr, lt, ye)[0];
    return exp(v) - *energy_shift;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_press_ye(const CCTK_REAL rho, const CCTK_REAL press,
                              const CCTK_REAL ye) const {

    printf(
        "This routine should not be used. There is no monotonicity condition "
        "to enforce a succesfull inversion from eps(press). So you better "
        "rewrite your code to not require this call. \n");
    assert(false);
    return 0;
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  csnd_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                              const CCTK_REAL ye) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL lt = std::log(std::fmin(std::fmax(temp, rgtemp.min), rgtemp.max));
    CCTK_REAL v = interptable->interpolate<EV::CS2>(lr, lt, ye)[0];
    assert(v >= 0);
    return sqrt(v);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  csnd_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                             const CCTK_REAL ye) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL lt = logtemp_from_valid_rho_eps_ye(rho, eps, ye);
    CCTK_REAL v = interptable->interpolate<EV::CS2>(lr, lt, ye)[0];
    assert(v >= 0);
    return sqrt(v);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  press_derivs_from_valid_rho_eps_ye(CCTK_REAL &press, CCTK_REAL &dpdrho,
                                     CCTK_REAL &dpdeps, const CCTK_REAL rho,
                                     const CCTK_REAL eps,
                                     const CCTK_REAL ye) const {
    printf("press_derivs_from_valid_rho_eps_ye is not supported for now! \n");
    assert(false);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  entropy_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                                 const CCTK_REAL ye) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL lt = std::log(std::fmin(std::fmax(temp, rgtemp.min), rgtemp.max));
    return interptable->interpolate<EV::S>(lr, lt, ye)[0];
  }

  CCTK_HOST CCTK_DEVICE inline void
  mu_pne_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                                 const CCTK_REAL ye, CCTK_REAL &mup, CCTK_REAL &mun, CCTK_REAL &mue) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL lt = std::log(std::fmin(std::fmax(temp, rgtemp.min), rgtemp.max));
    mup = interptable->interpolate<EV::MU_P>(lr, lt, ye)[0];
    mun = interptable->interpolate<EV::MU_N>(lr, lt, ye)[0];
    mue = interptable->interpolate<EV::MU_E>(lr, lt, ye)[0];
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_valid_rho_kappa_ye(const CCTK_REAL rho,
                                const CCTK_REAL kappa, // kappa=entropy
                                const CCTK_REAL ye) const {
    printf(
        "press_from_valid_rho_kappa_ye is not supported for tabulated EOS! \n");
    assert(false);
    return 0.0;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_kappa_ye(const CCTK_REAL rho,
                              const CCTK_REAL kappa, // kappa=entropy
                              const CCTK_REAL ye) const {
    printf(
        "eps_from_valid_rho_kappa_ye is not supported for tabulated EOS! \n");
    assert(false);
    return 0.0;
  };

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  kappa_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                              const CCTK_REAL ye) const {
    return entropy_from_valid_rho_temp_ye(
        rho, temp_from_valid_rho_eps_ye(rho, eps, ye), ye);
  }

  CCTK_HOST CCTK_DEVICE inline range
  range_eps_from_valid_rho_ye(const CCTK_REAL rho, const CCTK_REAL ye) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL vmin =
        interptable->interpolate<EV::EPS>(lr, interptable->xmin<1>(), ye)[0];
    CCTK_REAL vmax =
        interptable->interpolate<EV::EPS>(lr, interptable->xmax<1>(), ye)[0];
    return range{exp(vmin) - *energy_shift, exp(vmax) - *energy_shift};
  }

  CCTK_HOST CCTK_DEVICE inline range compute_eps_range_full_table() const {
    size_t n0 = interptable->num_points[0];
    size_t n1 = interptable->num_points[1];
    size_t n2 = interptable->num_points[2];
    size_t total = n0 * n1 * n2;
    const CCTK_REAL *eps_data = interptable->y + EV::EPS * total;
    CCTK_REAL eps_min = std::numeric_limits<CCTK_REAL>::max();
    CCTK_REAL eps_max = std::numeric_limits<CCTK_REAL>::lowest();
    for (size_t i = 0; i < total; i++) {
      CCTK_REAL val = exp(eps_data[i]) - *energy_shift;
      eps_min = std::fmin(eps_min, val);
      eps_max = std::fmax(eps_max, val);
    }
    eps_min = std::fmax(eps_min, 1e-15); // Force eps positive
    return range{eps_min, eps_max};
  }

}; // class eos_3p_tabulated3d

} // namespace EOSX

#endif // EOS_3P_TABULATED3D_HXX
