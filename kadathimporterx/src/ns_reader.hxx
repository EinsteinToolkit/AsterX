#ifndef _NS_READER_H
#define _NS_READER_H
#include "Configurator/config_bco.hpp"
#include "Solvers/ns_3d_xcts/ns_exporter.hpp"
#include "kadath_adapted.hpp"
#include "kadath_adapted_bh.hpp"

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

namespace NS_READER {
  using config_t = Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BCO_NS_INFO>;
  using reader_t = Kadath::FUKA_Solvers::CFMS_NS_Exporter;

  void malloc_ns_reader(std::string & fn);
  void free_ns_reader();
}
#endif
