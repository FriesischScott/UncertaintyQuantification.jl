# Metamodels

## Design Of Experiments

Design Of Experiments (DOE) offers various designs that can be used for creating a model of a given system. The core idea is to evaluate significant points of the system in order to obtain a sufficient model while keeping the effort to achieve this relatively low. Depending on the parameters, their individual importance and interconnections, different designs may be adequate.

The ones implemented here are `TwoLevelFactorial`, `FullFactorial`, `FractionalFactorial`, `CentralComposite`, `BoxBehnken` and `PlackettBurman`.

## Response Surface

A Response Surface is a simple polynomial surrogate model. It can be trained by providing it with evaluated points of a function or any of the aforementioned experimental designs.

### Example

In this example, we will model the following test function (known as Himmelblau's function) in the range ``x1, x2 ∈ [-5, 5]``. It is defined as

``f(x1, x2) = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2``.

![Image](metamodels_himmelblau.png)

We start by creating the random variables representing the inputs. These also define the domain over which the experimental designs are spread. We also define the original model using Himmelblau's function

```@example himmelblau
using UncertaintyQuantification # hide
x = RandomVariable.(Uniform(-5, 5), [:x1, :x2])

himmelblau = Model(df -> (df.x1 .^2 .+ df.x2 .- 11) .^2 .+ (df.x1 .+ df.x2 .^2 .- 7) .^2, :y)

nothing # hide
```

Our next step is to choose the design we want to use. In this example, we are applying a `FullFactorial` design with 5 levels in each dimension:

```@example himmelblau
design = FullFactorial([5, 5])
```

After that, we call the sample function with our design. This produces a `DataFrame` containing the points of our design spread over the range defined by the `RandomVariable`s. We then evaluate the model at these points and use the resulting data to train the `ResponseSurface`. The `ResponseSurface` uses least squares regression to fit a polynomial model to the given data points. The maximum polynomial degree is set in the `ResponseSurface` constructor.

> **_NOTE:_** The choice of the degree and the design and its parameters may be crucial to obtaining a sufficient model.

```@example himmelblau
training_data = sample(x, design)
evaluate!(himmelblau, training_data)
rs = ResponseSurface(training_data, :y, 4)
```

To evaluate the `ResponseSurface` run `evaluate!(rs::ResponseSurface, data::DataFrame)` with the `Dataframe` containing the points you want to evaluate. We validate the model by sampling 200 random points from the inputs and computing the mean squared error between the true model and the `ResponseSurface`.

```@example himmelblau
validation_data = sample(x, 200)
evaluate!(rs, validation_data)

y_h = himmelblau(validation_data)

mse = mean((y_h .- validation_data.y).^2)
```

The mse of ``1e26`` proves the accuracy of the metamodel. The following Figure compares the `ResponseSurface` to the true model over the entire domain.

![Image](metamodels_rs.png)
