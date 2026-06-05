# mvpt

An R package for learning about path value sensitivity across multiple, automatically generated SEMs. You can use mvpt() by supplying a lavaan-formatted model, a path to be tested, and data. 

Output includes the number of compared models (your model + all auto-generated models), the path value estimated in each model, and whether path value estimation significantly changed across these models. If no change is detected, path value estimation can be interpreted as robust despite the model specification differences among the compared models. But if change is detected, this suggests that a valuable model specification change may be possible. 

## Installation

Install the development version from GitHub:

```R
install.packages("remotes")
remotes::install_github("rdyeeflores/mvpt")

library(mvpt)
```

## Usage

To get an idea of how mvpt works using lavaan-formatted SEMs, consider the following example:

```R

data("UnfairApprais")

lavaan_input <-
  "
  ## latent variables
  Rumi =~ rumi1 + rumi2
  Angr =~ anger1 + anger1
  UnApp =~ unfair1 + unfair2
  ## regressions
  Aggr ~ Angr + Rumi
  Rumi ~ Angr + UnApp
  Angr ~ UnApp
  "

## path test 1
mvpt_path1 <- mvpt(lavaan_input, 
                    path = "Rumi~UnApp", 
                    data = UnfairApprais, 
                    showplots = TRUE)
mvpt_path1 

## path test 2
mvpt_path2 <- mvpt(lavaan_input, 
                    path = "Aggr~Rumi", 
                    data = UnfairApprais, 
                    showplots = TRUE)
mvpt_path2 

```