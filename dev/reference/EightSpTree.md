# Functional tree with 8 species.

This dataset is a leight-weight example.

## Usage

``` r
data(Paracou618)
```

## Format

An object of class
[`phylog`](https://adeverse.github.io/ade4/reference/phylog.html)
containing a functional tree.

## Examples

``` r
data(Paracou618)
# Preprocess the tree to be able to plot it 
# without loading ade4 package
plot(Preprocess.Tree(EightSpTree), hang=-0.01)
```
