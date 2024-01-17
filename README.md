<img align="top" src="Docs/figures/asterx.png" width="140">

**AsterX** is a GPU-accelerated GRMHD code for dynamical spacetimes, written in C++. It is built upon the [CarpetX](https://github.com/eschnett/CarpetX) driver, which is intended for the [Einstein Toolkit](https://einsteintoolkit.org/). **CarpetX** is based on [AMReX](https://amrex-codes.github.io), a software framework for block-structured AMR (adaptive mesh refinement).

Full documentation will soon be available at [asterx.readthedocs.io](https://asterx.readthedocs.io/en/latest/#).

* [![GitHub CI](https://github.com/jaykalinani/AsterX/workflows/CI/badge.svg)](https://github.com/jaykalinani/AsterX/actions)  [![Documentation Status](https://readthedocs.org/projects/asterx/badge/?version=latest)](https://asterx.readthedocs.io/en/latest/?badge=latest) [![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://github.com/jaykalinani/AsterX/blob/main/LICENSE.md)

## Overview 

* Heavily derived from the GRMHD code [Spritz](https://zenodo.org/record/4350072).
* Solves the GRMHD equations in 3D Cartesian coordinates and on dynamical spacetimes using high-resolution shock capturing (HRSC) schemes. 
* Based on the flux-conservative Valencia formulation.
* Directly evolves the staggered vector potential.

## Available modules 

* `AsterX` - the core GRMHD module
* `AsterSeeds` - initial data module
* `Con2PrimFactory` - module providing different conservative-to-primitive variable recovery routines
* `EOSX` - equation of state driver
* `ReconX` - provider of different reconstruction schemes
* `TOVSolverX` - a modified version of the publicly available TOVSolver thorn used within the Einstein Toolkit

## Getting started ##

Instructions for downloading and building the Einstein Toolkit including
CarpetX can be found [here](https://github.com/eschnett/CarpetX).

Details for building and running AsterX along with CarpetX will be added to [asterx.readthedocs.io](https://asterx.readthedocs.io/en/latest/#) soon..

## Related talks and tutorials ##

* "[Using CarpetX: A Guide for Early Adopters](http://einsteintoolkit.org/seminars/2021_03_18/index.html)". 
Recorded seminar talk by Erik Schnetter, providing an overview of the current capabilities of CarpetX.
* "[Tutorial: GPUs and the Einstein Toolkit](https://einsteintoolkit.github.io/et2022uidaho/lectures/38-Tutorial8/index.html)". 
Recorded tutorial by Lorenzo Ennoggi, Jay Kalinani and Federico Lopez Armengol during the North American Einstein Toolkit workshop 2022, presenting a brief overview on AsterX, followed by a hands-on session.
* "[AsterX: a new open-source GPU-accelerated GRMHD code for dynamical spacetimes](https://drive.google.com/file/d/1Z4i--W56mxeNIu598LQTpEEowX56FOoD/view?usp=sharing)". 
Slides based on the talk by Jay Kalinani at the APS April Meeting 2023.
