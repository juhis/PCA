# PCA
Principal component analysis

## Description

Calculates

* eigenvector decomposition for the given symmetric matrix or
* principal component scores based on given eigenvector matrix and original
matrix.

This program uses MTJ (Matrix Toolkits Java) which tries to utilize Java
Native Interface to run machine-optimized code. MTJ tells on runtime whether
using natives succeeds or not. It's a LOT faster with natives.

[MTJ](https://github.com/fommil/matrix-toolkits-java/)

[Benchmarks](http://lessthanoptimal.github.io/Java-Matrix-Benchmark/runtime/2013_10_Corei7v2600/)

## Requirements

Java 8

## Eigenvector decomposition

```
java -jar PCA.jar evd symmetricmatrixfile
```

The first command line argument should be "evd". The second command line
argument should be a path to a tab-delimited text file containing a symmetric
matrix. The first row in the file, as well as the first column on each row,
should contain headers. The input matrix is assumed to be ok.

Writes:

* eigenvectors.txt (a tab-delimited matrix: all right eigenvectors on rows)
* eigenvaluesReal.txt (real parts of eigenvalues, one per line)
* eigenvaluesImaginary.txt (imaginary parts of eigenvalues, one per line -- zero when run on a correlation matrix)

## Principal component scores

```
java -jar PCA.jar scores eigenvectorfile originalfile
```

The first command line argument should be "scores". The second command line
argument should be a path to a tab-delimited text file containing
eigenvectors on rows. The first row in the file, as well as the first column
on each row, should contain headers.

The third command line argument should be a path to a tab-delimited text file
containing the original data. The orientation is detected automatically so
that either rows or columns should correspond to the columns in the
eigenvector file. The first row in the file, as well as the first column on
each row, should contain headers.

The input matrices are assumed to be ok.

Writes:

* scores.txt (a tab-delimited matrix: all principal component scores on rows)
* cronbachsAlpha.txt (Cronbach's alpha for each component, one per line)
