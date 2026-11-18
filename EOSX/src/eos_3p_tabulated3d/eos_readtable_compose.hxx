#ifndef EOS_READTABLE_COMPOSE_HXX
#define EOS_READTABLE_COMPOSE_HXX

#include <cassert>
#include <cmath>
#include <limits>
#include <string>
#include <mpi.h>
#include <hdf5.h>

#include <cctk.h>
#include <cctk_Parameters.h>

#include "eos_3p_tabulated3d.hxx"

#ifndef NTABLES
#define NTABLES 19
#endif

namespace EOSX {
using namespace std;
using namespace eos_constants;

// Catch HDF5 errors
#define HDF5_ERROR(fn_call)                                                    \
  do {                                                                         \
    auto _error_code = (fn_call);                                              \
    if (_error_code < 0) {                                                     \
      printf("HDF5 call '%s' returned error code %lld\n", #fn_call,            \
             (long long)(_error_code));                                        \
      assert(false);                                                           \
    }                                                                          \
  } while (0)

#define READ_EOS_HDF5_COMPOSE(GROUP, NAME, VAR, TYPE, MEM)                     \
  do {                                                                         \
    hid_t dataset;                                                             \
    HDF5_ERROR(dataset = H5Dopen2(GROUP, NAME, H5P_DEFAULT));                  \
    HDF5_ERROR(H5Dread(dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR));        \
    HDF5_ERROR(H5Dclose(dataset));                                             \
  } while (0)

#define READ_ATTR_HDF5_COMPOSE(GROUP, NAME, VAR, TYPE)                         \
  do {                                                                         \
    hid_t attr;                                                                \
    HDF5_ERROR(attr = H5Aopen(GROUP, NAME, H5P_DEFAULT));                      \
    HDF5_ERROR(H5Aread(attr, TYPE, VAR));                                      \
    HDF5_ERROR(H5Aclose(attr));                                                \
  } while (0)

static inline int hdf5_link_exists(hid_t loc, const char *path) {
  const htri_t ex = H5Lexists(loc, path, H5P_DEFAULT);
  return (ex > 0);
}

CCTK_HOST inline void eos_readtable_compose(const std::string &filename,
                                            eos_tabulated3d_raw_table_t &tab) {

  CCTK_VINFO("Reading COMPOSE EOS table '%s'", filename.c_str());

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Pick the right MPI datatype for CCTK_REAL
  MPI_Datatype mpi_real = MPI_DOUBLE;
  if (sizeof(CCTK_REAL) == sizeof(float)) {
    mpi_real = MPI_FLOAT;
  } else if (sizeof(CCTK_REAL) == sizeof(double)) {
    mpi_real = MPI_DOUBLE;
  } else {
    assert(false && "Unsupported CCTK_REAL size for MPI_Bcast");
  }

  // Physical constants needed for conversions
  const double MeV_to_erg = 1.602176634e-6; // 1 MeV in erg
  const double cm3_to_fm3 = 1.0e39;         // (cm^-3) = (fm^-3) * 1e39
  const double c_cgs = 2.99792458e10;       // cm/s
  const double c2_cgs = c_cgs * c_cgs;
  const double m_neutron_MeV = 939.56542052; // neutron mass energy [MeV]
  const double baryon_mass = m_neutron_MeV * MeV_to_erg / c2_cgs; // g

  hid_t file = -1;
  hid_t parameters = -1;
  hid_t thermo_id = -1;

#ifdef H5_HAVE_PARALLEL
  hid_t fapl_id = -1;
#endif

  int nrho = 0, ntemp = 0, nye = 0;
  int nthermo = 0;

  auto *host_arena = amrex::The_Pinned_Arena();

#ifdef H5_HAVE_PARALLEL
  // ALL ranks do the same HDF5 calls to avoid deadlocks
  HDF5_ERROR(fapl_id = H5Pcreate(H5P_FILE_ACCESS));
  HDF5_ERROR(H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL));
  HDF5_ERROR(H5Pset_all_coll_metadata_ops(fapl_id, true));

  HDF5_ERROR(file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl_id));
  HDF5_ERROR(parameters = H5Gopen2(file, "/Parameters", H5P_DEFAULT));

  READ_ATTR_HDF5_COMPOSE(parameters, "pointsnb", &nrho, H5T_NATIVE_INT);
  READ_ATTR_HDF5_COMPOSE(parameters, "pointst", &ntemp, H5T_NATIVE_INT);
  READ_ATTR_HDF5_COMPOSE(parameters, "pointsyq", &nye, H5T_NATIVE_INT);

  HDF5_ERROR(thermo_id = H5Gopen2(file, "/Thermo_qty", H5P_DEFAULT));
  READ_ATTR_HDF5_COMPOSE(thermo_id, "pointsqty", &nthermo, H5T_NATIVE_INT);
#else
  if (rank == 0) {
    HDF5_ERROR(file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    HDF5_ERROR(parameters = H5Gopen2(file, "/Parameters", H5P_DEFAULT));

    READ_ATTR_HDF5_COMPOSE(parameters, "pointsnb", &nrho, H5T_NATIVE_INT);
    READ_ATTR_HDF5_COMPOSE(parameters, "pointst", &ntemp, H5T_NATIVE_INT);
    READ_ATTR_HDF5_COMPOSE(parameters, "pointsyq", &nye, H5T_NATIVE_INT);

    HDF5_ERROR(thermo_id = H5Gopen2(file, "/Thermo_qty", H5P_DEFAULT));
    READ_ATTR_HDF5_COMPOSE(thermo_id, "pointsqty", &nthermo, H5T_NATIVE_INT);
  }
#endif

  MPI_Bcast(&nrho, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ntemp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nye, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nthermo, 1, MPI_INT, 0, MPI_COMM_WORLD);

  const int npoints = nrho * ntemp * nye;

  CCTK_REAL *logrho_h =
      (CCTK_REAL *)host_arena->alloc(nrho * sizeof(CCTK_REAL));
  CCTK_REAL *logtemp_h =
      (CCTK_REAL *)host_arena->alloc(ntemp * sizeof(CCTK_REAL));
  CCTK_REAL *yes_h = (CCTK_REAL *)host_arena->alloc(nye * sizeof(CCTK_REAL));
  CCTK_REAL *alltables_h = (CCTK_REAL *)host_arena->alloc(
      (size_t)npoints * NTABLES * sizeof(CCTK_REAL));
  CCTK_REAL *energy_shift_h = (CCTK_REAL *)host_arena->alloc(sizeof(CCTK_REAL));

#ifdef H5_HAVE_PARALLEL
  {
    double *nb = (double *)host_arena->alloc(nrho * sizeof(double));
    double *t = (double *)host_arena->alloc(ntemp * sizeof(double));
    double *yq = (double *)host_arena->alloc(nye * sizeof(double));

    READ_EOS_HDF5_COMPOSE(parameters, "nb", nb, H5T_NATIVE_DOUBLE, H5S_ALL);
    READ_EOS_HDF5_COMPOSE(parameters, "t", t, H5T_NATIVE_DOUBLE, H5S_ALL);
    READ_EOS_HDF5_COMPOSE(parameters, "yq", yq, H5T_NATIVE_DOUBLE, H5S_ALL);

    for (int i = 0; i < nrho; i++) {
      const double rho_cgs = nb[i] * baryon_mass * cm3_to_fm3; // g/cm^3
      logrho_h[i] = (CCTK_REAL)log10(rho_cgs);
    }
    for (int i = 0; i < ntemp; i++) {
      logtemp_h[i] = (CCTK_REAL)log10(t[i]); // MeV
    }
    for (int i = 0; i < nye; i++) {
      yes_h[i] = (CCTK_REAL)yq[i];
    }

    int *thermo_index = (int *)host_arena->alloc(nthermo * sizeof(int));
    READ_EOS_HDF5_COMPOSE(thermo_id, "index_thermo", thermo_index,
                          H5T_NATIVE_INT, H5S_ALL);

    double *thermo_table =
        (double *)host_arena->alloc((size_t)nthermo * npoints * sizeof(double));
    READ_EOS_HDF5_COMPOSE(thermo_id, "thermo", thermo_table, H5T_NATIVE_DOUBLE,
                          H5S_ALL);

    constexpr int PRESS_C = 1;
    constexpr int S_C = 2;
    constexpr int MUN_C = 3;
    constexpr int MUP_C = 4;
    constexpr int MUE_C = 5;
    constexpr int EPS_C = 7;  // dimensionless in CompOSE
    constexpr int CS2_C = 12; // dimensionless

    auto const find_index = [&](int const &index) {
      for (int i = 0; i < nthermo; ++i) {
        if (thermo_index[i] == index)
          return i;
      }
      return -1;
    };

    const int iP = find_index(PRESS_C);
    const int iE = find_index(EPS_C);
    const int iS = find_index(S_C);
    const int iCS2 = find_index(CS2_C);
    const int iMUE = find_index(MUE_C);
    const int iMUP = find_index(MUP_C);
    const int iMUN = find_index(MUN_C);

    if (iP < 0 || iE < 0 || iS < 0 || iCS2 < 0 || iMUE < 0 || iMUP < 0 ||
        iMUN < 0) {
      CCTK_VError(
          __LINE__, __FILE__, CCTK_THORNSTRING,
          "Could not find required Thermo_qty indices in COMPOSE table");
    }

    // eps is dimensionless -> compute shift in dimensionless units
    double eps_min_dimless = std::numeric_limits<double>::max();
    for (int k = 0; k < nye; k++)
      for (int j = 0; j < ntemp; j++)
        for (int i = 0; i < nrho; i++) {
          const int idx = i + nrho * (j + ntemp * (k + nye * iE));
          eps_min_dimless = std::min(eps_min_dimless, thermo_table[idx]);
        }

    double energy_shift_dimless = 0.0;
    if (eps_min_dimless < 0.0) {
      energy_shift_dimless = -2.0 * eps_min_dimless;
    }

    // store shift in cgs erg/g so common post-processing can apply *EPSGF
    *energy_shift_h = (CCTK_REAL)(energy_shift_dimless * c2_cgs);

    for (size_t q = 0; q < (size_t)npoints * NTABLES; q++)
      alltables_h[q] = 0.0;

    for (int k = 0; k < nye; k++)
      for (int j = 0; j < ntemp; j++)
        for (int i = 0; i < nrho; i++) {

          const int ip = i + nrho * (j + ntemp * k);
          const size_t b = (size_t)ip * NTABLES;

          // pressure (MeV/fm^3 -> erg/cm^3 -> log10)
          {
            const int idx = i + nrho * (j + ntemp * (k + nye * iP));
            const double press_cgs =
                thermo_table[idx] * MeV_to_erg * cm3_to_fm3;
            alltables_h[b + eos_3p_tabulated3d::EV::PRESS] =
                (CCTK_REAL)log10(std::max(press_cgs, 1e-300));
          }

          // eps is dimensionless -> (eps + shift)*c^2 in erg/g -> log10
          {
            const int idx = i + nrho * (j + ntemp * (k + nye * iE));
            const double e_cgs =
                (thermo_table[idx] + energy_shift_dimless) * c2_cgs;
            alltables_h[b + eos_3p_tabulated3d::EV::EPS] =
                (CCTK_REAL)log10(std::max(e_cgs, 1e-300));
          }

          // entropy
          {
            const int idx = i + nrho * (j + ntemp * (k + nye * iS));
            alltables_h[b + eos_3p_tabulated3d::EV::S] =
                (CCTK_REAL)thermo_table[idx];
          }

          // cs2 (dimensionless -> cgs cm^2/s^2)
          {
            const int idx = i + nrho * (j + ntemp * (k + nye * iCS2));
            double cs2 = thermo_table[idx];
            if (cs2 < 0)
              cs2 = 0;
            if (cs2 > 0.9999999)
              cs2 = 0.9999999;
            alltables_h[b + eos_3p_tabulated3d::EV::CS2] =
                (CCTK_REAL)(cs2 * c2_cgs);
          }

          // chemical potentials mapping
          {
            const int idx_b = i + nrho * (j + ntemp * (k + nye * iMUN));
            const int idx_q = i + nrho * (j + ntemp * (k + nye * iMUP));
            const int idx_le = i + nrho * (j + ntemp * (k + nye * iMUE));

            const double mu_b = thermo_table[idx_b];
            const double mu_q = thermo_table[idx_q];
            const double mu_le = thermo_table[idx_le];

            const double mu_n = mu_b;
            const double mu_p = mu_b + mu_q;
            const double mu_e = mu_le - mu_q;

            alltables_h[b + eos_3p_tabulated3d::EV::MU_N] = (CCTK_REAL)mu_n;
            alltables_h[b + eos_3p_tabulated3d::EV::MU_P] = (CCTK_REAL)mu_p;
            alltables_h[b + eos_3p_tabulated3d::EV::MU_E] = (CCTK_REAL)mu_e;

            alltables_h[b + eos_3p_tabulated3d::EV::MUHAT] =
                (CCTK_REAL)(mu_n - mu_p);
          }
        }

    // Optional composition groups
    int ncomp = 0;
    if (hdf5_link_exists(file, "/Composition_pairs")) {
      hid_t comp_id;
      HDF5_ERROR(comp_id = H5Gopen2(file, "/Composition_pairs", H5P_DEFAULT));
      READ_ATTR_HDF5_COMPOSE(comp_id, "pointspairs", &ncomp, H5T_NATIVE_INT);

      int *index_yi = (int *)host_arena->alloc(ncomp * sizeof(int));
      READ_EOS_HDF5_COMPOSE(comp_id, "index_yi", index_yi, H5T_NATIVE_INT,
                            H5S_ALL);

      double *yi_table =
          (double *)host_arena->alloc((size_t)ncomp * npoints * sizeof(double));
      READ_EOS_HDF5_COMPOSE(comp_id, "yi", yi_table, H5T_NATIVE_DOUBLE,
                            H5S_ALL);

      auto const find_index_yi = [&](int const &index) {
        for (int i = 0; i < ncomp; ++i) {
          if (index_yi[i] == index)
            return i;
        }
        return -1;
      };

      const int i_n = find_index_yi(10);
      const int i_p = find_index_yi(11);
      const int i_a = find_index_yi(4002);

      for (int k = 0; k < nye; k++)
        for (int j = 0; j < ntemp; j++)
          for (int i = 0; i < nrho; i++) {
            const int ip = i + nrho * (j + ntemp * k);
            const size_t b = (size_t)ip * NTABLES;

            if (i_n >= 0)
              alltables_h[b + eos_3p_tabulated3d::EV::XN] =
                  (CCTK_REAL)yi_table[ip + (size_t)npoints * i_n];
            if (i_p >= 0)
              alltables_h[b + eos_3p_tabulated3d::EV::XP] =
                  (CCTK_REAL)yi_table[ip + (size_t)npoints * i_p];
            if (i_a >= 0)
              alltables_h[b + eos_3p_tabulated3d::EV::XA] =
                  (CCTK_REAL)(4.0 * yi_table[ip + (size_t)npoints * i_a]);
          }

      host_arena->free(index_yi);
      host_arena->free(yi_table);
      HDF5_ERROR(H5Gclose(comp_id));
    }

    int nav = 0;
    if (hdf5_link_exists(file, "/Composition_quadrupels")) {
      hid_t av_id;
      HDF5_ERROR(av_id =
                     H5Gopen2(file, "/Composition_quadrupels", H5P_DEFAULT));
      READ_ATTR_HDF5_COMPOSE(av_id, "pointsav", &nav, H5T_NATIVE_INT);

      if (nav > 0) {
        assert(nav == 1 &&
               "nav != 1 in this table, so there is none or more than "
               "one definition of an average nucleus."
               "Please check and generalize accordingly.");

        double *zav =
            (double *)host_arena->alloc((size_t)npoints * sizeof(double));
        double *yav =
            (double *)host_arena->alloc((size_t)npoints * sizeof(double));
        double *aav =
            (double *)host_arena->alloc((size_t)npoints * sizeof(double));

        READ_EOS_HDF5_COMPOSE(av_id, "zav", zav, H5T_NATIVE_DOUBLE, H5S_ALL);
        READ_EOS_HDF5_COMPOSE(av_id, "yav", yav, H5T_NATIVE_DOUBLE, H5S_ALL);
        READ_EOS_HDF5_COMPOSE(av_id, "aav", aav, H5T_NATIVE_DOUBLE, H5S_ALL);

        for (int k = 0; k < nye; k++)
          for (int j = 0; j < ntemp; j++)
            for (int i = 0; i < nrho; i++) {
              const int ip = i + nrho * (j + ntemp * k);
              const size_t b = (size_t)ip * NTABLES;

              alltables_h[b + eos_3p_tabulated3d::EV::ABAR] =
                  (CCTK_REAL)aav[ip];
              alltables_h[b + eos_3p_tabulated3d::EV::ZBAR] =
                  (CCTK_REAL)zav[ip];
              alltables_h[b + eos_3p_tabulated3d::EV::XH] =
                  (CCTK_REAL)(aav[ip] * yav[ip]);
            }

        host_arena->free(zav);
        host_arena->free(yav);
        host_arena->free(aav);
      }

      HDF5_ERROR(H5Gclose(av_id));
    }

    HDF5_ERROR(H5Gclose(thermo_id));
    HDF5_ERROR(H5Gclose(parameters));
    HDF5_ERROR(H5Fclose(file));
    HDF5_ERROR(H5Pclose(fapl_id));

    host_arena->free(nb);
    host_arena->free(t);
    host_arena->free(yq);
    host_arena->free(thermo_index);
    host_arena->free(thermo_table);
  }
#else
  if (rank == 0) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Serial-HDF5 path for COMPOSE reader is not enabled in this "
                "header version.");
  }
#endif

  // Broadcast tables
  MPI_Bcast(logrho_h, nrho, mpi_real, 0, MPI_COMM_WORLD);
  MPI_Bcast(logtemp_h, ntemp, mpi_real, 0, MPI_COMM_WORLD);
  MPI_Bcast(yes_h, nye, mpi_real, 0, MPI_COMM_WORLD);
  MPI_Bcast(energy_shift_h, 1, mpi_real, 0, MPI_COMM_WORLD);
  MPI_Bcast(alltables_h, npoints * NTABLES, mpi_real, 0, MPI_COMM_WORLD);

  // Final device/managed storage
  CCTK_REAL *logrho =
      (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(nrho * sizeof(CCTK_REAL));
  CCTK_REAL *logtemp =
      (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(ntemp * sizeof(CCTK_REAL));
  CCTK_REAL *yes =
      (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(nye * sizeof(CCTK_REAL));
  CCTK_REAL *alltables = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
      (size_t)npoints * NTABLES * sizeof(CCTK_REAL));
  CCTK_REAL *energy_shift =
      (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(sizeof(CCTK_REAL));

  amrex::Gpu::copy(amrex::Gpu::hostToDevice, logrho_h, logrho_h + nrho, logrho);
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, logtemp_h, logtemp_h + ntemp,
                   logtemp);
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, yes_h, yes_h + nye, yes);
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, energy_shift_h, energy_shift_h + 1,
                   energy_shift);
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, alltables_h,
                   alltables_h + (size_t)npoints * NTABLES, alltables);

  amrex::Gpu::synchronize();

  host_arena->free(logrho_h);
  host_arena->free(logtemp_h);
  host_arena->free(yes_h);
  host_arena->free(alltables_h);
  host_arena->free(energy_shift_h);

  tab.nrho = nrho;
  tab.ntemp = ntemp;
  tab.nye = nye;
  tab.npoints = npoints;

  // CompOSE cs2 is already relativistic
  tab.have_rel_cs2 = 1;

  tab.logrho = logrho;
  tab.logtemp = logtemp;
  tab.yes = yes;
  tab.alltables = alltables;
  tab.energy_shift = energy_shift;

  return;
}

} // namespace EOSX

#endif // EOS_READTABLE_COMPOSE_HXX
