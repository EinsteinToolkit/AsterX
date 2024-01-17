# Install FUKA with AsterX (Frontera)

FUKA thorns are available (privately) at:

* https://bitbucket.org/fukaws/kadaththorn/
* https://bitbucket.org/fukaws/kadathimporterx/

Ask Samuel Tootle for access.
To load the FUKA library, source the FUKA environment with the following command before compiling:

* `source /work2/09803/sdtootle/fuka_export/env.sh`

Correctly set the KADATH, BLAS and LAPACK DIR in the optionlist:
```
KADATH_DIR=${HOME_KADATH}
BLAS_DIR=${MKLROOT}/lib
BLAS_LIBS = -Wl,--start-group  mkl_gf_lp64 mkl_sequential mkl_core mkl_blacs_intelmpi_lp64 -Wl,--end-group mkl_blas95_lp64  pthread m mkl_scalapack_lp64 dl
LAPACK_DIR=${MKLROOT}/lib
LAPACK_LIBS= mkl_lapack95_lp64
```

Add FUKA thorns to the thornlist:
```
#KADATH thorns

!TARGET   = $ARR
!TYPE     = git
!URL      = https://bitbucket.org/fukaws/$2
!REPO_BRANCH = $ET_RELEASE
!REPO_PATH=../$2
!CHECKOUT =
Fuka/KadathImporterX
Fuka/KadathThorn

# FUKA initial data thorns
# the crazy path works around a bug in GetComponents that does not handle
# symbolic links and ".." correctly
!TARGET = $ARR/Fuka/KadathThorn/src/fuka/../../../../repos/KadathThorn/src/fuka
!TYPE   = git
!URL    = https://bitbucket.org/fukaws/fuka
!CHECKOUT = Cmake build_debug build_release codes eos include install_par.sh install_seq.sh src src_par src_seq
```

Build using simfactory:
* `sim build -j28 asterxFUKA --thornlist=./thornlists/asterx.th --machine=frontera_Intel19 SILENT=NO 2>&1 | tee make_out`

Rebuild to get `SILENT=NO` option activated using `make`:
* `make -j28 asterxFUKA-rebuild SILENT=NO 2>&1 | tee make_out`

To correct for the linking of the mkl libraries, manually copy the link command (last line) from `make_out` to another file, say `make_out_mod`, and execute using the command:
* `sh make_out_mod`

