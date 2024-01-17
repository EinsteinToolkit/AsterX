#include "ns_reader.hxx"

namespace NS_READER {
    NS_READER::reader_t* ns_reader{nullptr};

    void malloc_ns_reader(std::string & fn) {
        ns_reader = new NS_READER::reader_t(fn);
    }
    void free_ns_reader() {
        delete ns_reader;
    }
    extern "C" void KadathImporterNS(CCTK_ARGUMENTS) { 
        DECLARE_CCTK_ARGUMENTS;
        DECLARE_CCTK_PARAMETERS;
        using namespace Loop;
        
        using importer_t = NS_READER::reader_t;
        static std::mutex copy_mutex;
        // Copy ID object per thread (hopefully)        
        copy_mutex.lock();
        importer_t importer(*ns_reader);
        copy_mutex.unlock();
        
        const GF3D<CCTK_REAL, 0, 0, 0> gxx_(cctkGH, gxx);
        const GF3D<CCTK_REAL, 0, 0, 0> gxy_(cctkGH, gxy);
        const GF3D<CCTK_REAL, 0, 0, 0> gxz_(cctkGH, gxz);
        const GF3D<CCTK_REAL, 0, 0, 0> gyy_(cctkGH, gyy);
        const GF3D<CCTK_REAL, 0, 0, 0> gyz_(cctkGH, gyz);
        const GF3D<CCTK_REAL, 0, 0, 0> gzz_(cctkGH, gzz);

        const GF3D<CCTK_REAL, 0, 0, 0> Kxx_(cctkGH, kxx);
        const GF3D<CCTK_REAL, 0, 0, 0> Kxy_(cctkGH, kxy);
        const GF3D<CCTK_REAL, 0, 0, 0> Kxz_(cctkGH, kxz);
        const GF3D<CCTK_REAL, 0, 0, 0> Kyy_(cctkGH, kyy);
        const GF3D<CCTK_REAL, 0, 0, 0> Kyz_(cctkGH, kyz);
        const GF3D<CCTK_REAL, 0, 0, 0> Kzz_(cctkGH, kzz);

        const GF3D<CCTK_REAL, 0, 0, 0> alp_(cctkGH, alp);

        const GF3D<CCTK_REAL, 0, 0, 0> betax_(cctkGH, betax);
        const GF3D<CCTK_REAL, 0, 0, 0> betay_(cctkGH, betay);
        const GF3D<CCTK_REAL, 0, 0, 0> betaz_(cctkGH, betaz);

        const GF3D<CCTK_REAL, 0, 0, 0> dtalp_(cctkGH, dtalp);

        const GF3D<CCTK_REAL, 0, 0, 0> dtbetax_(cctkGH, dtbetax);
        const GF3D<CCTK_REAL, 0, 0, 0> dtbetay_(cctkGH, dtbetay);
        const GF3D<CCTK_REAL, 0, 0, 0> dtbetaz_(cctkGH, dtbetaz);

        loop_all<0, 0, 0>(cctkGH, [&, ns_reader](const PointDesc &p) {
            double const xx = p.x;
            double const yy = p.y;
            double const zz = p.z;

            auto vars = std::move(importer.export_pointwise_spacetime_vars(xx, yy, zz));

            gxx_(p.I) = vars[importer_t::OUTPUT_VARS::GXX];
            gxy_(p.I) = vars[importer_t::OUTPUT_VARS::GXY];
            gxz_(p.I) = vars[importer_t::OUTPUT_VARS::GXZ];
            gyy_(p.I) = vars[importer_t::OUTPUT_VARS::GYY];
            gyz_(p.I) = vars[importer_t::OUTPUT_VARS::GYZ];
            gzz_(p.I) = vars[importer_t::OUTPUT_VARS::GZZ];

            Kxx_(p.I) = vars[importer_t::OUTPUT_VARS::KXX];
            Kxy_(p.I) = vars[importer_t::OUTPUT_VARS::KXY];
            Kxz_(p.I) = vars[importer_t::OUTPUT_VARS::KXZ];
            Kyy_(p.I) = vars[importer_t::OUTPUT_VARS::KYY];
            Kyz_(p.I) = vars[importer_t::OUTPUT_VARS::KYZ];
            Kzz_(p.I) = vars[importer_t::OUTPUT_VARS::KZZ];
            alp_(p.I) = vars[importer_t::OUTPUT_VARS::ALPHA];
            betax_(p.I) = 0;
            betay_(p.I) = 0;
            betaz_(p.I) = 0;
            dtalp_(p.I) = 0;
            dtbetax_(p.I) = 0;
            dtbetay_(p.I) = 0;
            dtbetaz_(p.I) = 0;
        });

        const GF3D<CCTK_REAL, 1, 1, 1> rho_(cctkGH, rho);
        const GF3D<CCTK_REAL, 1, 1, 1> eps_(cctkGH, eps);
        const GF3D<CCTK_REAL, 1, 1, 1> press_(cctkGH, press);
        const GF3D<CCTK_REAL, 1, 1, 1> velx_(cctkGH, velx);
        const GF3D<CCTK_REAL, 1, 1, 1> vely_(cctkGH, vely);
        const GF3D<CCTK_REAL, 1, 1, 1> velz_(cctkGH, velz);
        
        loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {

            double const xx = p.x;
            double const yy = p.y;
            double const zz = p.z;
            
            auto vars = std::move(importer.export_pointwise_fluid_vars(xx, yy, zz));

            rho_(p.I)   = vars[importer_t::OUTPUT_VARS::RHO];
            eps_(p.I)   = vars[importer_t::OUTPUT_VARS::EPS];
            press_(p.I) = vars[importer_t::OUTPUT_VARS::PRESS];
            velx_(p.I)  = vars[importer_t::OUTPUT_VARS::VELX];
            vely_(p.I)  = vars[importer_t::OUTPUT_VARS::VELY];
            velz_(p.I)  = vars[importer_t::OUTPUT_VARS::VELZ];

        });
    }
}
