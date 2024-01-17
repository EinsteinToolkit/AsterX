# Install CarpetX with Spack (Frontera)

Use interactive session

* Compile CPU version: `idev -m 120`

* Compile GPU version: `idev -p rtx-dev -m 120`

Download spack

* `git clone -c feature.manyFiles=true https://github.com/spack/spack.git`

* `git checkout relesases/v0.19`

* `. share/spack/setup-env.sh`

Install gcc@11.2.0

* `spack compiler find`

* `spack install gcc@11.2.0 %gcc@4.8.5`

* `spack compiler add ...` (`...` is the last line of previous command)

## Download CarpetX

```
curl -kLO https://raw.githubusercontent.com/gridaphobe/CRL/ET_2022_05/GetComponents
chmod a+x GetComponents
./GetComponents --root Cactus --parallel --no-shallow https://bitbucket.org/eschnett/cactusamrex/raw/rhaas-thornlist/manifest/carpetx.th
```

## Install the CPU version

Create a dir where you want put `view` in (say `/home1/.../username/Cactus_gcc`)

* replace the last line of `CPU/spack_yaml` with your dir (say `/home1/.../username/Cactus_gcc/view`)

* replace the dir `/home1/08708/liwei/Cactus_gcc/view` (with say `/home1/.../username/Cactus_gcc/view`)
in `config_frontera_gcc-11.2.0.cfg`

Install other required packages

* `env TMPDIR=$WORK/tmp spack --env-dir ./CPU compiler find`

* `env TMPDIR=$WORK/tmp spack --env-dir ./CPU concretize --force`

* `env TMPDIR=$WORK/tmp spack --env-dir ./CPU install --fail-fast`

Install CarpetX

* `spack load gcc@11.2.0`

* `cd Cactus`

* `gmake CarpetX-gcc options=config_frontera_gcc-11.2.0.cfg`

* `cp repos/cactusamrex/azure-pipelines/carpetx.th configs/CarpetX-gcc/ThornList`

* `gmake -j16 CarpetX-gcc`


## Install the GPU version

Create a dir where you want put `view` in (say `/home1/.../username/Cactus_cuda`)

* replace the last line of `GPU/spack_yaml` with your dir (say `/home1/.../username/Cactus_cuda/view`)

* replace the dir `/home1/08708/liwei/Cactus_cuda/view` (with say `/home1/.../username/Cactus_cuda/view`)
in `config_frontera_cuda-11.5.2.cfg`

Install other required packages

* `env TMPDIR=$WORK/tmp spack --env-dir ./GPU compiler find view-cuda-compilers`

* `env TMPDIR=$WORK/tmp spack --env-dir ./GPU concretize --force`

* `env TMPDIR=$WORK/tmp spack --env-dir ./GPU install --fail-fast`

Install CarpetX

* `spack load gcc@11.2.0`

* `spack load cuda@11.5.2`

* `cd Cactus`

* `gmake CarpetX-cuda options=config_frontera_cuda-11.5.2.cfg`

* `cp repos/cactusamrex/azure-pipelines/carpetx.th configs/CarpetX-cuda/ThornList`

* `gmake -j16 CarpetX-cuda`



