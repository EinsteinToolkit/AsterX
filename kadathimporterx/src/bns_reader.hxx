#ifndef _BNS_READER_H
#define _BNS_READER_H
#include "Configurator/config_binary.hpp"
#include "Solvers/bns_xcts/bns_exporter.hpp"
#include "kadath_bin_ns.hpp"

#include <cmath>
#include <memory>

#include <loop.hxx>
#include <vec.hxx>
#include <vect.hxx>
#include <mat.hxx>
#include <omp.h>
#include <mutex>
#ifdef _STANDALONE_
#include <iostream>
#else
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#endif

namespace BNS_READER {
  using config_t = Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BIN_INFO>;
  using reader_t = Kadath::FUKA_Solvers::CFMS_BNS_Exporter;

  void malloc_bns_reader(std::string & fn);
  void free_bns_reader();
}
#endif
