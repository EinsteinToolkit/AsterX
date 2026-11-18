#include <cctype>
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

#include "AMReX.H"

// EOSX Stuff
#include "setup_eos.hxx"
namespace ID_TabEOS_HydroQuantities {

using namespace amrex;
using namespace std;

class Ye_reader {
private:
  double *Ye_rho_arr;
  double *logrho_arr;
  FILE *in1D;
  int nrho;

public:
  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  init(int init_nrho, FILE *init_in1D, const double *init_logrho_arr) {
    nrho = init_nrho;
    in1D = init_in1D;

    if (!(logrho_arr =
              (double *)The_Managed_Arena()->alloc(nrho * sizeof(double)))) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Cannot allocate memory for EOS table");
    }
    for (int i = 0; i < nrho; i++)
      logrho_arr[i] = init_logrho_arr[i];

    read_1dfile__set_array();
  }

  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void release() {
    The_Managed_Arena()->free(Ye_rho_arr);
    The_Managed_Arena()->free(logrho_arr);
    Ye_rho_arr = nullptr;
    logrho_arr = nullptr;
  }

  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  read_1dfile__set_array(const int num_header_lines = 0) {
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int which_line = 0;
    int which_value = 0;

    if (!(Ye_rho_arr =
              (double *)The_Managed_Arena()->alloc(nrho * sizeof(double)))) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Cannot allocate memory for EOS table");
    }

    while ((read = getline(&line, &len, in1D)) != -1) {
      if (which_line++ < num_header_lines)
        continue;

      if (which_value >= nrho)
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "Y_e file contains more than %d entries", nrho);

      char *end;
      const double value = strtod(line, &end);
      while (*end != '\0' &&
             std::isspace(static_cast<unsigned char>(*end)))
        ++end;

      if (end == line || *end != '\0' || !std::isfinite(value))
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "Invalid Y_e value on data line %d", which_value + 1);

      Ye_rho_arr[which_value++] = value;
    }
    free(line);

    if (which_value != nrho)
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Y_e file contains %d entries; expected %d", which_value,
                  nrho);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  interpolate_1d_quantity_as_function_of_rho(const CCTK_REAL rho) const {
    const CCTK_REAL logrho = log(rho);

    if (logrho <= logrho_arr[0])
      return Ye_rho_arr[0];
    if (logrho >= logrho_arr[nrho - 1])
      return Ye_rho_arr[nrho - 1];

    int lower = 0;
    int upper = nrho - 1;
    while (upper - lower > 1) {
      const int middle = (lower + upper) / 2;
      if (logrho < logrho_arr[middle])
        upper = middle;
      else
        lower = middle;
    }

    const CCTK_REAL weight =
        (logrho - logrho_arr[lower]) /
        (logrho_arr[upper] - logrho_arr[lower]);
    return Ye_rho_arr[lower] * (1.0 - weight) +
           Ye_rho_arr[upper] * weight;
  }
};

} // namespace ID_TabEOS_HydroQuantities
