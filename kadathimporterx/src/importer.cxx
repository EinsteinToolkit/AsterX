#include <cmath>
#include <functional>
#include <iostream>
#include <memory>

#include <loop.hxx>
#include <vec.hxx>
#include <vect.hxx>
#include <mat.hxx>
#include <omp.h>
#include <mutex>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "Solvers/exporter.hpp"
#include "bns_reader.hxx"
#include "bh_reader.hxx"
#include "ns_reader.hxx"

namespace KadathImporterX {

    extern "C" void KadathImporterX_SetImporter(CCTK_ARGUMENTS) {
        DECLARE_CCTK_ARGUMENTS;
        DECLARE_CCTK_PARAMETERS;
        CCTK_INFO("Setting up KADATH initial data");
        
        std::string id_type{type};
        std::string fn{filename};
        std::string ss{"Initializing reader for: "+id_type + " " + fn};
        CCTK_INFO(ss.c_str());
        
        if(id_type == "BH") {
          BH_READER::malloc_bh_reader(fn);
        } else if(id_type == "NS") {
          NS_READER::malloc_ns_reader(fn);
        } else if(id_type == "BNS") {
          BNS_READER::malloc_bns_reader(fn);
        } else {
          CCTK_PARAMWARN ("id_type not yet implemented");
        }
    }
    extern "C" void KadathImporterX_FreeImporter(CCTK_ARGUMENTS) {
        DECLARE_CCTK_ARGUMENTS;
        DECLARE_CCTK_PARAMETERS;
        CCTK_INFO("Cleaning up KADATH initial data importer");
        
        std::string id_type{type};
        if(id_type == "BH") {
          BH_READER::free_bh_reader();
        } else if(id_type == "NS") {
          NS_READER::free_ns_reader();
        } else if(id_type == "NS") {
          BNS_READER::free_bns_reader();
        } else {
          CCTK_PARAMWARN ("id_type not yet implemented");
        }
    }    
}
