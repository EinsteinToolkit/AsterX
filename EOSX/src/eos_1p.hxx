#ifndef EOS1P_HXX
#define EOS1P_HXX
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "eos_utils.hxx"

namespace EOSX {

/// Abstract EOS class for implementations of degenerate (cold) equations of
/// state.
// sed = specific energy density
// rmd = rest mass density

/**
\par
This holds a pointer with shared ownership to the actual implementation.
Therefore it can be copied around by value without wasting resources
(e.g. in case of tabulated EOS) and without the need to manage resources.
\par Notation
It is assumed that the unit sytem is geometric, i.e. \f$ c=G=1 \f$,
but the exact choice is left to the implementation and the code using it.

<table border>
<tr>
<td> \f$ P \f$      </td><td>Pressure              </td><td><tt>p</tt></td>
</tr><tr>
<td> \f$ \rho_E \f$ </td><td>Total energy density  </td><td><tt>ed</tt></td>
</tr><tr>
<td> \f$ n_B \f$    </td><td>Baryon number density </td><td> </td>
</tr><tr>
<td> \f$ \rho = m_B n_B \f$ </td><td>Formal rest mass density
</td><td><tt>rmd</tt></td>
</tr><tr>
<td> \f$ m_B \f$    </td><td>Formal baryon mass    </td><td></td>
</tr><tr>
<td> \f$ \epsilon = \rho_E / \rho -1 \f$ </td>
<td>Specific internal energy</td><td><tt>sed</tt></td>
</tr><tr>
<td> \f$ \rho_I = \rho_E - \rho \f$ </td>
<td>Internal energy density</td><td><tt>ied</tt></td>
</tr><tr>
<td> \f$ h = 1+ \epsilon + \frac{P}{\rho} \f$ </td>
<td>Specific enthalpy including restmass</td>
<td><tt>h</tt></td>
</tr><tr>
<td> \f$ h-1 \f$ </td>
<td>Specific enthalpy excluding restmass</td>
<td><tt>hm1</tt></td>
</tr><tr>
<td>\f$ g = \exp\left(\int_{P(\rho=0)}^P \frac{dP'}{\rho_E(P') + P'}\right)
\f$</td> <td>Hydrostatic potential</td><td><tt>g</tt></td>
</tr><tr>
<td>\f$ g - 1 \f$</td>
<td></td><td><tt>gm1</tt></td>
</tr><tr>
<td>\f$ c_s^2 \f$</td><td>Squared adiabatic
soundspeed</td><td><tt>csnd2</tt></td>
</tr>
</table>

\par General EOS properties
The EOSs are parametrized by density and alternatively by the quantity
\f$ g-1 \f$ (which is useful to compute hydrostatic equilibrium).
At zero density, \f$ g-1=0 \f$.
For isentropic EOSs, \f$ g=h \f$.
No assumption is made that the EOS is isentropic, one can have e.g. isothermal
as well. Still, the implementation needs to provide that information.
\par
All implementations have to define a validity range. Calls to the EOS
with input outside the valid range throws an exception.
The mass constant \f$ m_B \f$ is not specified. The EOS only deals with
rest mass density, not baryon number density.
Note the specific internal energy and enthalpy are formal quantities
depending on this choice, only \f$ \rho h, \rho_E \f$ have physical meaning.
*/

class eos_1p {
public:
  typedef eos_range range;

private:
  bool isentropic;  ///< Whether EOS is isentropic
  bool temp_avail;  ///< Whether EOS has temperature information
  bool efrac_avail; ///< Whether EOS has electron fraction information
  range rg_rmd;     ///< Density range in which EOS is valid.
  range rg_gm1;     ///< Range of \f$ g-1 \f$ in which EOS is valid.
  range rg_p;       ///< Pressure range in which EOS is valid.

public:
  /// Set valid density range. Has to be called by constructor of each
  /// implementation.
  CCTK_DEVICE CCTK_HOST void set_ranges(const range &rg_rmd_);
  /// Check if density is in valid range, else throw exception.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  check_rmd_valid(CCTK_REAL rmd) const;
  /// Check if \f$ g-1 \f$ is in valid range, else throw exception.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  check_gm1_valid(CCTK_REAL gm1) const;
  /// Check if pressure is in valid range, else throw exception.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  check_p_valid(CCTK_REAL p) const;

  /// Whether EOS is isentropic
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_isentropic() const {
    return isentropic;
  }
  /// Whether EOS can compute temperature
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  has_temp() const {
    return temp_avail;
  }
  /// Whether EOS can compute electron fraction
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  has_efrac() const {
    return efrac_avail;
  }
  /// Returns range of validity for density
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline const range &
  range_rmd() const {
    return rg_rmd;
  }
  /// Returns range of validity for \f$ g-1 \f$
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline const range &
  range_gm1() const {
    return rg_gm1;
  }
  /// Returns range of validity for pressure
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline const range &
  range_p() const {
    return rg_p;
  }
  /// Whether given density is in the valid range of the EOS
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_rmd_valid(CCTK_REAL rmd) const {
    return rg_rmd.contains(rmd);
  }
  /// Whether given \f$ g-1 \f$ is in the valid range of the EOS
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_gm1_valid(CCTK_REAL gm1) const {
    return rg_gm1.contains(gm1);
  }
  /// Whether given pressure is in the valid range of the EOS
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_p_valid(CCTK_REAL p) const {
    return rg_p.contains(p);
  }
};

} // namespace EOSX

#endif
