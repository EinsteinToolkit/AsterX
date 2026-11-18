#ifndef EOS_READTABLE_SCOLLAPSE_HXX
#define EOS_READTABLE_SCOLLAPSE_HXX

#include <cassert>
#include <string>
#include <mpi.h>
#include <hdf5.h>

#include "eos_3p_tabulated3d.hxx"

#ifndef NTABLES
#define NTABLES 19
#endif

namespace EOSX {
using namespace std;
using namespace eos_constants;

CCTK_HOST inline void
eos_readtable_scollapse(const std::string &filename,
                        eos_tabulated3d_raw_table_t &tab) {

  CCTK_VINFO("Reading Stellar Collapse EOS table '%s'", filename.c_str());

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

#ifdef H5_HAVE_PARALLEL
  get_hdf5_int_dset(file_id, "pointsrho", 1, &nrho);
  get_hdf5_int_dset(file_id, "pointstemp", 1, &ntemp);
  get_hdf5_int_dset(file_id, "pointsye", 1, &nye);
#else
  if (rank == 0) {
    get_hdf5_int_dset(file_id, "pointsrho", 1, &nrho);
    get_hdf5_int_dset(file_id, "pointstemp", 1, &ntemp);
    get_hdf5_int_dset(file_id, "pointsye", 1, &nye);
  }
  MPI_Bcast(&nrho, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ntemp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nye, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  const int npoints = nrho * ntemp * nye;

  int have_rel_cs2 = 0;
#ifdef H5_HAVE_PARALLEL
  have_rel_cs2 = (H5Lexists(file_id, "/have_rel_cs2", H5P_DEFAULT) > 0);
#else
  if (rank == 0) {
    have_rel_cs2 = (H5Lexists(file_id, "/have_rel_cs2", H5P_DEFAULT) > 0);
  }
  MPI_Bcast(&have_rel_cs2, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  // Read and communicate tables using host memory
  auto *host_arena = amrex::The_Pinned_Arena();

  CCTK_REAL *logrho_h =
      (CCTK_REAL *)host_arena->alloc(nrho * sizeof(CCTK_REAL));
  CCTK_REAL *logtemp_h =
      (CCTK_REAL *)host_arena->alloc(ntemp * sizeof(CCTK_REAL));
  CCTK_REAL *yes_h = (CCTK_REAL *)host_arena->alloc(nye * sizeof(CCTK_REAL));
  CCTK_REAL *alltables_tmp_h = (CCTK_REAL *)host_arena->alloc(
      (size_t)npoints * NTABLES * sizeof(CCTK_REAL));
  CCTK_REAL *alltables_h = (CCTK_REAL *)host_arena->alloc(
      (size_t)npoints * NTABLES * sizeof(CCTK_REAL));
  CCTK_REAL *energy_shift_h = (CCTK_REAL *)host_arena->alloc(sizeof(CCTK_REAL));

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

  static const char *dnames[NTABLES] = {
      "logpress", "logenergy", "entropy", "munu", "cs2",  "dedt", "dpdrhoe",
      "dpderho",  "muhat",     "mu_e",    "mu_p", "mu_n", "Xa",   "Xh",
      "Xn",       "Xp",        "Abar",    "Zbar", "gamma"};

#ifdef H5_HAVE_PARALLEL
  // ALL ranks participate in reads in the same order
  get_hdf5_real_dset(file_id, "logrho", nrho, logrho_h);
  get_hdf5_real_dset(file_id, "logtemp", ntemp, logtemp_h);
  get_hdf5_real_dset(file_id, "ye", nye, yes_h);
  get_hdf5_real_dset(file_id, "energy_shift", 1, energy_shift_h);
  for (int iv = 0; iv < NTABLES; iv++) {
    get_hdf5_real_dset(file_id, dnames[iv], npoints,
                       &alltables_tmp_h[(size_t)iv * npoints]);
  }

  // change ordering of alltables array so that
  // the table kind is the fastest changing index
  for (int iv = 0; iv < NTABLES; iv++)
    for (int k = 0; k < nye; k++)
      for (int j = 0; j < ntemp; j++)
        for (int i = 0; i < nrho; i++) {
          int indold = i + nrho * (j + ntemp * (k + nye * iv));
          int indnew = iv + NTABLES * (i + nrho * (j + ntemp * k));
          alltables_h[indnew] = alltables_tmp_h[indold];
        }
  host_arena->free(alltables_tmp_h);

  CHECK_ERROR(H5Fclose(file_id));
  CHECK_ERROR(H5Pclose(fapl_id));
#else
  if (rank == 0) {
    get_hdf5_real_dset(file_id, "logrho", nrho, logrho_h);
    get_hdf5_real_dset(file_id, "logtemp", ntemp, logtemp_h);
    get_hdf5_real_dset(file_id, "ye", nye, yes_h);
    get_hdf5_real_dset(file_id, "energy_shift", 1, energy_shift_h);
    for (int iv = 0; iv < NTABLES; iv++) {
      get_hdf5_real_dset(file_id, dnames[iv], npoints,
                         &alltables_tmp_h[(size_t)iv * npoints]);
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
            alltables_h[indnew] = alltables_tmp_h[indold];
          }
  }
  host_arena->free(alltables_tmp_h);

  MPI_Bcast(logrho_h, nrho, mpi_real, 0, MPI_COMM_WORLD);
  MPI_Bcast(logtemp_h, ntemp, mpi_real, 0, MPI_COMM_WORLD);
  MPI_Bcast(yes_h, nye, mpi_real, 0, MPI_COMM_WORLD);
  MPI_Bcast(energy_shift_h, 1, mpi_real, 0, MPI_COMM_WORLD);
  MPI_Bcast(alltables_h, npoints * NTABLES, mpi_real, 0, MPI_COMM_WORLD);
#endif

  // Copy host -> managed
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

  tab.have_rel_cs2 = have_rel_cs2;

  tab.logrho = logrho;
  tab.logtemp = logtemp;
  tab.yes = yes;
  tab.alltables = alltables;
  tab.energy_shift = energy_shift;

  return;
}

} // namespace EOSX

#endif // EOS_READTABLE_SCOLLAPSE_HXX
