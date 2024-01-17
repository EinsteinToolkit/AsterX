# RITAsterX


## ToDo

* [ ] Cell-centered metric

    - pros: less interpolation, better interface to SphericalNR, easier to differentiate

    - cons: prolongation operation, will be public supported

    - links:

        - https://github.com/lwJi/SpacetimeX/tree/cell-centered-metric

        - https://github.com/lwJi/AsterX/tree/cell-centered-metric

* [ ] Verify if the order of gridient of metrics in the source term matters (Spritz can do both 2nd and 4th order derivatives)

* [ ] WenoZ and MP5 (debug, see CarpetX)

* [ ] Neutrino Transport

    - Neutrino Leakage

        - Port ZelmaniLeak thorn

    - M1

        - Port David's thorn

    - Full Boltzmann equation
