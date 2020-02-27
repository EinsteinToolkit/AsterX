#include "driver.hxx"
#include "reduction.hxx"
#include "schedule.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

namespace CarpetX {
using namespace amrex;
using namespace std;

namespace {

// Interpolate a grid function at one point, dimensionally recursive
template <typename T, int order> struct interpolator {
  const Array4<const T> &vars;
  const int vi;
  const vect<int, dim> &derivs;
  const T *restrict const dx;

#warning "TODO: Check whether interpolated variables are valid"

  // Base case: only access a grid point
  template <int dir>
  enable_if_t<(dir == -1), T>
  interpolate(const vect<int, dim> &i, const vect<CCTK_REAL, dim> &di) const {
    return vars(i[0], i[1], i[2], vi);
  }

  // General case: interpolate in one direction, then recurse
  template <int dir>
  enable_if_t<(dir >= 0), T> interpolate(const vect<int, dim> &i,
                                         const vect<CCTK_REAL, dim> &di) const {
    const auto DI = vect<int, dim>::unit(dir);
    switch (order) {
    case 0: {
      const T x = di[dir];
#ifdef CCTK_DEBUG
      assert(fabs(x) <= T(0.5));
#endif
      const T y0 = interpolate<dir - 1>(i, di);
      switch (derivs[dir]) {
      case 0:
        return y0;
      case 1:
        return 0;
      case 2:
        return 0;
      }
    }
    case 1: {
      const T x = di[dir];
#ifdef CCTK_DEBUG
      assert(x >= T(0) && x <= T(1));
#endif
      const T y0 = interpolate<dir - 1>(i, di);
      const T y1 = interpolate<dir - 1>(i + DI, di);
      switch (derivs[dir]) {
      case 0:
        return (1 - x) * y0 + x * y1;
      case 1:
        return (-y0 + y1) / dx[dir];
      case 2:
        return 0;
      }
    }
    case 2: {
      const T x = di[dir] - order / T(2);
#ifdef CCTK_DEBUG
      assert(fabs(x) <= T(0.5));
#endif
      const T y0 = interpolate<dir - 1>(i, di);
      const T y1 = interpolate<dir - 1>(i + DI, di);
      const T y2 = interpolate<dir - 1>(i + 2 * DI, di);
      switch (derivs[dir]) {
      case 0:
        return (-1 / T(2) * x + 1 / T(2) * pow(x, 2)) * y0 +
               (1 - pow(x, 2)) * y1 +
               (1 / T(2) * x + 1 / T(2) * pow(x, 2)) * y2;
      case 1:
        return ((-1 / T(2) + x) * y0 - 2 * x * y1 + (1 / T(2) + x) * y2) /
               dx[dir];
      case 2:
        return (y0 - 2 * y1 + y2) / pow(dx[dir], 2);
      }
    }
    case 3: {
      const T x = di[dir] - order / T(2);
#ifdef CCTK_DEBUG
      assert(fabs(x) <= T(0.5));
#endif
      const T y0 = interpolate<dir - 1>(i, di);
      const T y1 = interpolate<dir - 1>(i + DI, di);
      const T y2 = interpolate<dir - 1>(i + 2 * DI, di);
      const T y3 = interpolate<dir - 1>(i + 3 * DI, di);
      switch (derivs[dir]) {
      case 0:
        return (-1 / T(16) + 1 / T(24) * x + 1 / T(4) * pow(x, 2) -
                1 / T(6) * pow(x, 3)) *
                   y0 +
               (9 / T(16) - 9 / T(8) * x - 1 / T(4) * pow(x, 2) +
                1 / T(2) * pow(x, 3)) *
                   y1 +
               (9 / T(16) + 9 / T(8) * x - 1 / T(4) * pow(x, 2) -
                1 / T(2) * pow(x, 3)) *
                   y2 +
               (-1 / T(16) - 1 / T(24) * x + 1 / T(4) * pow(x, 2) +
                1 / T(6) * pow(x, 3)) *
                   y3;
      case 1:
        return ((1 / T(24) + 1 / T(2) * x - 1 / T(2) * pow(x, 2)) * y0 +
                (-9 / T(8) - 1 / T(2) * x + 3 / T(2) * pow(x, 2)) * y1 +
                (9 / T(8) - 1 / T(2) * x - 3 / T(2) * pow(x, 2)) * y2 +
                (-1 / T(24) + 1 / T(2) * x + 1 / T(2) * pow(x, 2)) * y3) /
               dx[dir];
      case 2:
        return ((1 / T(2) - x) * y0 + (-1 / T(2) + 3 * x) * y1 +
                (-1 / T(2) - 3 * x) * y2 + (1 / T(2) + x) * y3) /
               pow(dx[dir], 2);
      }
    }
    case 4: {
      const T x = di[dir] - order / T(2);
#ifdef CCTK_DEBUG
      assert(fabs(x) <= T(0.5));
#endif
      const T y0 = interpolate<dir - 1>(i, di);
      const T y1 = interpolate<dir - 1>(i + DI, di);
      const T y2 = interpolate<dir - 1>(i + 2 * DI, di);
      const T y3 = interpolate<dir - 1>(i + 3 * DI, di);
      const T y4 = interpolate<dir - 1>(i + 4 * DI, di);
      switch (derivs[dir]) {
      case 0:
        return (1 / T(12) * x - 1 / T(24) * pow(x, 2) - 1 / T(12) * pow(x, 3) +
                1 / T(24) * pow(x, 4)) *
                   y0 +
               (-2 / T(3) * x + 2 / T(3) * pow(x, 2) + 1 / T(6) * pow(x, 3) -
                1 / T(6) * pow(x, 4)) *
                   y1 +
               (1 - 5 / T(4) * pow(x, 2) + 1 / T(4) * pow(x, 4)) * y2 +
               (2 / T(3) * x + 2 / T(3) * pow(x, 2) - 1 / T(6) * pow(x, 3) -
                1 / T(6) * pow(x, 4)) *
                   y3 +
               (-1 / T(12) * x - 1 / T(24) * pow(x, 2) + 1 / T(12) * pow(x, 3) +
                1 / T(24) * pow(x, 4)) *
                   y4;
      case 1:
        return ((1 / T(12) - 1 / T(12) * x - 1 / T(4) * pow(x, 2) +
                 1 / T(6) * pow(x, 3)) *
                    y0 +
                (-2 / T(3) + 4 / T(3) * x + 1 / T(2) * pow(x, 2) -
                 2 / T(3) * pow(x, 3)) *
                    y1 +
                (-5 / T(2) * x + pow(x, 3)) * y2 +
                (2 / T(3) + 4 / T(3) * x - 1 / T(2) * pow(x, 2) -
                 2 / T(3) * pow(x, 3)) *
                    y3 +
                (-1 / T(12) - 1 / T(12) * x + 1 / T(4) * pow(x, 2) +
                 1 / T(6) * pow(x, 3)) *
                    y4) /
               dx[dir];
      case 2:
        return ((-1 / T(12) - 1 / T(2) * x + 1 / T(2) * pow(x, 2)) * y0 +
                (4 / T(3) + x - 2 * pow(x, 2)) * y1 +
                (-5 / T(2) + 3 * pow(x, 2)) * y2 +
                (4 / T(3) - x - 2 * pow(x, 2)) * y3 +
                (-1 / T(12) + 1 / T(2) * x + 1 / T(2) * pow(x, 2)) * y4) /
               pow(dx[dir], 2);
      }
    }
    default:
      assert(0);
    }
  }
};

} // namespace

extern "C" void CarpetX_Interpolate(const CCTK_POINTER_TO_CONST cctkGH_,
                                    const CCTK_INT npoints,
                                    const CCTK_REAL *restrict const coordsx,
                                    const CCTK_REAL *restrict const coordsy,
                                    const CCTK_REAL *restrict const coordsz,
                                    const CCTK_INT nvars,
                                    const CCTK_INT *restrict const varinds,
                                    const CCTK_INT *restrict const operations,
                                    const CCTK_POINTER resultptrs_) {
  DECLARE_CCTK_PARAMETERS;
  const cGH *restrict const cctkGH = static_cast<const cGH *>(cctkGH_);
  assert(in_global_mode(cctkGH));

  // Create particle container
  typedef AmrParticleContainer<0, 2> Container;
  Container container(ghext->amrcore.get());

  // Set particle positions
  {
    const int level = 0;
    const auto &restrict leveldata = ghext->leveldata.at(level);
    const MFIter mfi(*leveldata.mfab0);
    assert(mfi.isValid());
    auto &particles = container.GetParticles(
        level)[make_pair(mfi.index(), mfi.LocalTileIndex())];
    // particles.GetArrayOfStructs()().reserve(npoints);
    for (int n = 0; n < npoints; ++n) {
      Particle<0, 2> p;
      p.id() = Container::ParticleType::NextID();
      p.cpu() = ParallelDescriptor::MyProc();
      p.pos(0) = coordsx[n];
      p.pos(1) = coordsy[n];
      p.pos(2) = coordsz[n];
      p.idata(0) = ParallelDescriptor::MyProc(); // source process
      p.idata(1) = n;                            // source index
      particles.push_back(p);
    }
  }

  // Send particles to interpolation points
  container.Redistribute();

  // Define result variables
  map<int, vector<CCTK_REAL> > results;

  // Interpolate
  constexpr int tl = 0;
  struct givi_t {
    int gi, vi;
  };
  vector<givi_t> givis(nvars);
  for (int v = 0; v < nvars; ++v) {
    int gi = CCTK_GroupIndexFromVarI(varinds[v]);
    assert(gi >= 0);
    assert(gi < CCTK_NumGroups());
    int vi = varinds[v] - CCTK_FirstVarIndexI(gi);
    assert(vi >= 0);
    assert(vi < CCTK_NumVarsInGroupI(gi));
    givis.at(v) = {gi, vi};
  }

  for (const auto &leveldata : ghext->leveldata) {
    const int level = leveldata.level;
#warning "TODO: use OpenMP"
    for (ParIter<0, 2> pti(container, level); pti.isValid(); ++pti) {
      const Geometry &geom = ghext->amrcore->Geom(level);
      const CCTK_REAL *restrict const x0 = geom.ProbLo();
      const CCTK_REAL *restrict const dx = geom.CellSize();

      const int np = pti.numParticles();
      const auto &particles = pti.GetArrayOfStructs();

      vector<vector<CCTK_REAL> > varresults(nvars);

      for (int v = 0; v < nvars; ++v) {
        const int gi = givis.at(v).gi;
        const int vi = givis.at(v).vi;
        const auto &restrict groupdata = *leveldata.groupdata.at(gi);
        // Ensure interpolated variables are vertex centred
        // TODO: Generalize this
        assert((groupdata.indextype == array<int, dim>{0, 0, 0}));
        const Array4<const CCTK_REAL> &vars = groupdata.mfab.at(tl)->array(pti);
        vect<int, dim> derivs;
        int op = operations[v];
        while (op > 0) {
          const int dir = op % 10 - 1;
          if (dir >= 0) {
            assert(dir >= 0 && dir < dim);
            ++derivs[dir];
          }
          op /= 10;
        }
        auto &varresult = varresults.at(v);
        varresult.resize(np);

        switch (interpolation_order) {
        case 0: {
          constexpr int order = 0;
#pragma omp simd
          for (int n = 0; n < np; ++n) {
            vect<int, dim> i;
            vect<CCTK_REAL, dim> di;
            for (int d = 0; d < dim; ++d) {
              CCTK_REAL x = particles[n].pos(d);
              CCTK_REAL ri = (x - x0[d]) / dx[d];
              i[d] = lrint(ri - (order / CCTK_REAL(2)));
              di[d] = ri - i[d];
            }
            const interpolator<CCTK_REAL, order> interp{vars, vi, derivs, dx};
            varresult.at(n) = interp.interpolate<dim - 1>(i, di);
          }
          break;
        }
        case 1: {
          constexpr int order = 1;
#pragma omp simd
          for (int n = 0; n < np; ++n) {
            vect<int, dim> i;
            vect<CCTK_REAL, dim> di;
            for (int d = 0; d < dim; ++d) {
              CCTK_REAL x = particles[n].pos(d);
              CCTK_REAL ri = (x - x0[d]) / dx[d];
              i[d] = lrint(ri - (order / CCTK_REAL(2)));
              di[d] = ri - i[d];
            }
            const interpolator<CCTK_REAL, order> interp{vars, vi, derivs, dx};
            varresult.at(n) = interp.interpolate<dim - 1>(i, di);
          }
          break;
        }
        case 2: {
          constexpr int order = 2;
#pragma omp simd
          for (int n = 0; n < np; ++n) {
            vect<int, dim> i;
            vect<CCTK_REAL, dim> di;
            for (int d = 0; d < dim; ++d) {
              CCTK_REAL x = particles[n].pos(d);
              CCTK_REAL ri = (x - x0[d]) / dx[d];
              i[d] = lrint(ri - (order / CCTK_REAL(2)));
              di[d] = ri - i[d];
            }
            const interpolator<CCTK_REAL, order> interp{vars, vi, derivs, dx};
            varresult.at(n) = interp.interpolate<dim - 1>(i, di);
          }
          break;
        }
        case 3: {
          constexpr int order = 3;
#pragma omp simd
          for (int n = 0; n < np; ++n) {
            vect<int, dim> i;
            vect<CCTK_REAL, dim> di;
            for (int d = 0; d < dim; ++d) {
              CCTK_REAL x = particles[n].pos(d);
              CCTK_REAL ri = (x - x0[d]) / dx[d];
              i[d] = lrint(ri - (order / CCTK_REAL(2)));
              di[d] = ri - i[d];
            }
            const interpolator<CCTK_REAL, order> interp{vars, vi, derivs, dx};
            varresult.at(n) = interp.interpolate<dim - 1>(i, di);
          }
          break;
        }
        case 4: {
          constexpr int order = 4;
#pragma omp simd
          for (int n = 0; n < np; ++n) {
            vect<int, dim> i;
            vect<CCTK_REAL, dim> di;
            for (int d = 0; d < dim; ++d) {
              CCTK_REAL x = particles[n].pos(d);
              CCTK_REAL ri = (x - x0[d]) / dx[d];
              i[d] = lrint(ri - (order / CCTK_REAL(2)));
              di[d] = ri - i[d];
            }
            const interpolator<CCTK_REAL, order> interp{vars, vi, derivs, dx};
            varresult.at(n) = interp.interpolate<dim - 1>(i, di);
          }
          break;
        }
        default:
          assert(0);
        }
      }

      for (int n = 0; n < np; ++n) {
        const int proc = particles[n].idata(0);
        const int id = particles[n].idata(1);
        if (!results.count(proc))
          results[proc];
        auto &result = results.at(proc);
        result.push_back(id);
        for (int v = 0; v < nvars; ++v)
          result.push_back(varresults.at(v).at(n));
      }
    }
  }

  // Collect particles back
  const int nprocs = ParallelDescriptor::NProcs();
  const MPI_Comm comm = ParallelDescriptor::Communicator();
  const MPI_Datatype datatype = mpi_datatype<CCTK_REAL>::value;

  vector<int> sendcounts(nprocs, 0);
  for (const auto &proc_result : results) {
    const int p = proc_result.first;
    const auto &result = proc_result.second;
    sendcounts.at(p) = result.size();
  }
  vector<int> senddispls(nprocs);
  int sendcount = 0;
  for (int p = 0; p < nprocs; ++p) {
    senddispls.at(p) = sendcount;
    sendcount += sendcounts.at(p);
  }
  vector<int> recvcounts(nprocs);
  MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT,
               comm);
  vector<int> recvdispls(nprocs);
  int recvcount = 0;
  for (int p = 0; p < nprocs; ++p) {
    recvdispls.at(p) = recvcount;
    recvcount += recvcounts.at(p);
  }
  assert(recvcount == (nvars + 1) * npoints);
  vector<CCTK_REAL> sendbuf(sendcount);
  for (const auto &proc_result : results) {
    const int p = proc_result.first;
    const auto &result = proc_result.second;
    copy(result.begin(), result.end(), &sendbuf.at(senddispls.at(p)));
  }
  vector<CCTK_REAL> recvbuf(recvcount);
  MPI_Alltoallv(sendbuf.data(), sendcounts.data(), senddispls.data(), datatype,
                recvbuf.data(), recvcounts.data(), recvdispls.data(), datatype,
                comm);
#ifdef CCTK_DEBUG
  // Check consistency of received ids
  vector<bool> idxs(npoints, false);
  for (int n = 0; n < npoints; ++n) {
    const int offset = (nvars + 1) * n;
    const int idx = int(recvbuf.at(offset));
    assert(!idxs.at(idx));
    idxs.at(idx) = true;
  }
  for (int n = 0; n < npoints; ++n)
    assert(idxs.at(n));
#endif

  // Set result
  CCTK_REAL *const restrict *const restrict resultptrs =
      static_cast<CCTK_REAL *const *const restrict>(resultptrs_);
  for (int n = 0; n < npoints; ++n) {
    const int offset = (nvars + 1) * n;
    const int idx = int(recvbuf.at(offset));
    for (int v = 0; v < nvars; ++v)
      resultptrs[v][idx] = recvbuf.at(offset + 1 + v);
  }
}
} // namespace CarpetX
