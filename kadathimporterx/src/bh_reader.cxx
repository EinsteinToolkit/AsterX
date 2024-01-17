#include "bh_reader.hxx"

namespace BH_READER {
    BH_READER::reader_t* bh_reader{nullptr};

    void malloc_bh_reader(std::string & fn) {
        bh_reader = new BH_READER::reader_t(fn);
    }
    void free_bh_reader() {
        delete bh_reader;
    }
    extern "C" void KadathImporterBH(CCTK_ARGUMENTS) { 
        DECLARE_CCTK_ARGUMENTS;
        DECLARE_CCTK_PARAMETERS;
        using namespace Loop;
        
        using importer_t = BH_READER::reader_t;
        static std::mutex copy_mutex;
        // Copy ID object per thread (hopefully)        
        copy_mutex.lock();
        importer_t importer(*bh_reader);
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

        loop_all<0, 0, 0>(cctkGH, [&, bh_reader](const PointDesc &p) {
            double const xx = p.x;
            double const yy = p.y;
            double const zz = p.z;

            auto vars = std::move(importer.export_pointwise(
              xx, yy, zz, 
              interpolation_offset, interp_order, delta_r_rel
              )
            );

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
    }
}
