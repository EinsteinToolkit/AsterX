#ifndef ASTERX_SYNC_HXX
#define ASTERX_SYNC_HXX

#include <cctk.h>

#include <vector>

// Thin wrappers around the CarpetX driver's level window and its
// restriction / prolongation entry points. They are defined in sync.cxx, the
// one AsterX file that calls CarpetX's public API (carpetx.hxx), so other
// files (in particular gauge_register.cxx) can express "every time-aligned
// level pair" without reaching into the driver themselves.
//
// Under subcycling CarpetX calls global functions with
// CarpetX::active_window() set to the window [min_level, max_level) of levels
// that are at the current time. Without subcycling every level is in the
// window.

namespace AsterX {

// The active window reaches level 0: every level is at the same time.
bool full_cascade();

// Restrict `groups` fine-to-coarse onto every level that has an aligned
// child (a child in the active window), finest pair first, without poisoning
// the restricted region. The driver's restrict tag does not apply: a
// restrict="no" group (IGr) is restricted too.
void RestrictFromAlignedChildren(const std::vector<int> &groups);

// Same-level ghost exchange (plus outer boundary conditions) of `groups` on
// every active level. No prolongation.
void SyncGhostsOnly(const std::vector<int> &groups);

} // namespace AsterX

#endif // #ifndef ASTERX_SYNC_HXX
