# Metamodels
## Design Of Experiments

Design Of Experiments (DOE) offers various designs that can be used for creating a model of a given system. The core idea is to evaluate significant points of the system in order to obtain a sufficient model while keeping the effort to achieve this relatively low.
Depending on the parameters, their individual importance and interconnections, different designs may be adequate.

The ones implemented here are two level factorial, full factorial, fractional factorial and box behnken. (TODO: link to)

## Response Surface

A Response Surface is a structure used for modeling. 
It can be trained by providing it with evaluated points of a function.
It will then, using polynomial regression, compute a model of that function.

## Example

In this example, we will model the following variation of the ishigami function in the interval $[-\pi, \pi]$.

$f(x) = sin(x_1) + a * sin(x_2) ^2 + b * x_3 ^4 * sin(x_1)$

