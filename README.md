# mvpt

An R package for learning about path value sensitivity across multiple, automatically generated SEMs. You can use mvpt() by supplying a lavaan-formatted model, a path to be tested, and data. 

Ouput will show how many models were compared (your model + all auto-generated models), the path value estimated in each model, and whether path value estimation significantly changed across these models. If no change is detected, path value estimation can be interpreted as robust despite the model specification differences among the compared models. But if change is detected, this is indicative that a valuable model specification change could be possible. 

## Installation

Install the development version from GitHub:

```R
# install.packages("remotes")
remotes::install_github("rdyeeflores/mvpt")
```

## Usage

To an idea of how mvpt works using lavaan-formatted SEMs, consider the following example using the Political Democracy Example from 
Bollen's 1989 book:

```R
library(lavaan)

lavaan_model <- 
"
   # latent variables
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4
     dem65 =~ y5 + y6 + y7 + y8
   # regressions
     dem60 ~ ind60
     dem65 ~ ind60 + dem60
"
path <- "dem60 ~ ind60"
mvpt <- mvpt(lavaan_model, path, data = PoliticalDemocracy)
mvpt

```