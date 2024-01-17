#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <algorithm>
#include <cmath>

namespace AsterX {
using namespace std;

/* calculate max of d(var)/(dx) * hx in all dirs */
template <typename T>
CCTK_DEVICE inline T calc_grad_1st(const Loop::GF3D2<const T> &gf,
                                   const Loop::PointDesc &p) {
  CCTK_REAL err{0}, errp{0}, errm{0};
  for (int d = 0; d < Loop::dim; ++d) {
    auto varm = gf(p.I - p.DI[d]);
    auto var0 = gf(p.I);
    auto varp = gf(p.I + p.DI[d]);
    errp += (varp - var0) * (varp - var0);
    errm += (var0 - varm) * (var0 - varm);
  }
  err = max({errp, errm});
  return sqrt(err);
}

/* calculate max of d(var)/(dx) * hx in all dirs */
template <typename T>
CCTK_DEVICE inline T calc_deriv_1st(const Loop::GF3D2<const T> &gf,
                                    const Loop::PointDesc &p) {
  CCTK_REAL err{0};
  for (int d = 0; d < Loop::dim; ++d) {
    auto varm = gf(p.I - p.DI[d]);
    auto var0 = gf(p.I);
    auto varp = gf(p.I + p.DI[d]);
    err = max({err, fabs(var0 - varm), fabs(varp - var0)});
  }
  return err;
}

/* calculate max of d^2(var)/(dx^2) * hx^2 in all dirs */
template <typename T>
CCTK_DEVICE inline T calc_deriv_2nd(const Loop::GF3D2<const T> &gf,
                                    const Loop::PointDesc &p) {
  CCTK_REAL err{0};
  for (int d = 0; d < Loop::dim; ++d) {
    auto varm = gf(p.I - p.DI[d]);
    auto var0 = gf(p.I);
    auto varp = gf(p.I + p.DI[d]);
    err = max({err, fabs(varp + varm - 2 * var0)});
  }
  return err;
}

/* calculate max of second derivative error norm in all dirs
 * based on eq. 50 of Mignone+2011 (https://arxiv.org/abs/1110.0740) */
template <typename T>
CCTK_DEVICE inline T calc_deriv_2nd_norm(const Loop::GF3D2<const T> &gf,
                                         const Loop::PointDesc &p,
                                         CCTK_REAL epsilon_err) {
  CCTK_REAL err{0};
  for (int d = 0; d < Loop::dim; ++d) {
    auto varm = gf(p.I - p.DI[d]);
    auto var0 = gf(p.I);
    auto varp = gf(p.I + p.DI[d]);

    auto diffp0 = varp - var0;
    auto diff0m = var0 - varm;

    auto varsigma = fabs(varp) + 2 * fabs(var0) + fabs(varm);

    auto diffpm2 = fabs(diffp0 - diff0m) * fabs(diffp0 - diff0m);
    auto diffpm_denom = (fabs(diffp0) + fabs(diff0m) + epsilon_err * varsigma) *
                        (fabs(diffp0) + fabs(diff0m) + epsilon_err * varsigma);

    auto chi = sqrt(diffpm2 / diffpm_denom);
    err = max({err, chi});
  }
  return err;
}

/* functions for interpreting pars */
void read_stream(vector<string> &groups, istringstream &groupstream) {
  string delim = " \n", group;
  while (getline(groupstream, group)) {
    size_t prev = 0, pos;
    while ((pos = group.find_first_of(delim, prev)) != string::npos) {
      if (pos > prev)
        groups.push_back(group.substr(prev, pos - prev));
      prev = pos + 1;
    }
    if (prev < group.length())
      groups.push_back(group.substr(prev, string::npos));
  }
}

array<int, Loop::dim> get_group_indextype(const int gi) {
  DECLARE_CCTK_PARAMETERS;

  assert(gi >= 0);
  const int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  array<CCTK_INT, Loop::dim> index;
  int iret = Util_TableGetIntArray(tags, Loop::dim, index.data(), "index");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    // Fallback: use centering table
    const int centering = CCTK_GroupCenteringTableI(gi);
    assert(centering >= 0);
    iret =
        Util_TableGetIntArray(centering, Loop::dim, index.data(), "centering");
  }
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    // Default: vertex-centred
    index = {0, 0, 0};
  } else if (iret >= 0) {
    assert(iret == Loop::dim);
  } else {
    assert(0);
  }

  // Convert to index type
  array<int, Loop::dim> indextype;
  for (int d = 0; d < Loop::dim; ++d)
    indextype[d] = index[d];

  return indextype;
}

} // namespace AsterX
