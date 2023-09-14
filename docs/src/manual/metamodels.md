# Metamodels

## Design Of Experiments

Design Of Experiments (DOE) offers various designs that can be used for creating a model of a given system. The core idea is to evaluate significant points of the system in order to obtain a sufficient model while keeping the effort to achieve this relatively low. Depending on the parameters, their individual importance and interconnections, different designs may be adequate.

The ones implemented here are `TwoLevelFactorial`, `FullFactorial`, `FractionalFactorial`, `CentralComposite`, `BoxBehnken` and `PlackettBurman`.

## Response Surface

A Response Surface is a simple polynomial surrogate model. It can be trained by providing it with evaluated points of a function or any of the aforementioned experimental designs.