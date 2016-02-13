# PCA
Principal component analysis

## Description

* Center, scale and transpose data matrices
* Calculate covariance and correlation matrices
* Calculate eigenvector decomposition for symmetric matrices (such as covariance and correlation matrices)
* Calculate principal component scores based on an eigenvector matrix and an original data matrix.
* Change orientation of PCA results: Calculate "transposed eigenvectors" under certain circumstances

This program uses MTJ (Matrix Toolkits Java) which tries to utilize Java Native Interface to run machine-optimized code. MTJ tells on runtime whether using natives succeeds or not. It's a LOT faster with natives.

[MTJ](https://github.com/fommil/matrix-toolkits-java/)

[Benchmarks of Java matrix libraries](http://lessthanoptimal.github.io/Java-Matrix-Benchmark/runtime/2013_10_Corei7v2600/)

## Requirements

Java 8

## Input data format

All commands use one or more input matrices. These should be tab-delimited text files. They can be plain text or gzipped. The first row in a data file, as well as the first column on each row, should contain headers. The input matrices are assumed to be ok: No missing values, no unparseable values, all rows of the same length.

## Output data format

Output matrices have the same format as input matrices. Output matrices are written in plain text.

## Data centering

```
java -jar PCA.jar center-rows file
```

Writes:

* file.centered-rows.txt (a tab-delimited matrix with zero mean for every row)

```
java -jar PCA.jar center-columns file
```

Writes:

* file.centered-columns.txt (a tab-delimited matrix with zero mean for every column)

## Data scaling

```
java -jar PCA.jar scale-rows file
```

Writes:

* file.scaled-rows.txt (a tab-delimited matrix with standard deviation one for every row)

```
java -jar PCA.jar scale-columns file
```

Writes:

* file.scaled-columns.txt (a tab-delimited matrix with standard deviation one for every column)

## Transposition

```
java -jar PCA.jar transpose file
```

The first command line argument should be "transpose". The second command line argument should be a path to a tab-delimited text file containing a matrix.

Writes:

* file.transposed.txt (a tab-delimited matrix with the data transposed)

## Covariance

```
java -jar PCA.jar covariance-rows file
```

Writes:

* file.covariance-rows.txt (a tab-delimited symmetric covariance matrix with covariance calculated for each pair of rows)

```
java -jar PCA.jar covariance-columns file
```

Writes:

* file.covariance-columns.txt (a tab-delimited symmetric covariance matrix with covariance calculated for each pair of columns)

## Correlation

```
java -jar PCA.jar correlation-rows file
```

Writes:

* file.correlation-rows.txt (a tab-delimited symmetric correlation matrix with correlation calculated for each pair of rows)

```
java -jar PCA.jar correlation-columns file
```

Writes:

* file.correlation-columns.txt (a tab-delimited symmetric correlation matrix with correlation calculated for each pair of columns)

## Eigenvector decomposition

```
java -jar PCA.jar evd symmetricmatrixfile
```

Writes:

* file.eigenvectors.txt (a tab-delimited matrix: all right eigenvectors on rows)
* file.eigenvalues.txt (real parts of eigenvalues, one per line)

## Principal component scores

```
java -jar PCA.jar scores eigenvectorfile originalfile
```

eigenvectorfile should contain eigenvectors on rows.
The orientation of originalfile is determined automatically.

Writes:

* originalfile.scores.txt (a tab-delimited matrix: all principal component scores on rows)
* originalfile.cronbachsAlpha.txt (Cronbach's alpha for each component, one per line)

## Transformation of orientation: "Transposed eigenvectors" from principal component scores and eigenvalues

```
java -jar PCA.jar transform scorefile eigenvaluefile
```

scorefile should contain principal component scores on rows (i.e. each row corresponds to a component).
eigenvalues in eigenvaluefile should correspond to the eigenvectors based on which the principal component scores have been calculated.

*Note: This operation works in the following scenario:*

* There is an original input matrix X
* Each column of X has been centered (zero mean for each column)
* A covariance matrix has been calculated over rows (for pairs of rows) of X
* Eigenvector decomposition has been done for this covariance matrix
* Principal component scores have been calculated based on the resulting eigenvectors

Then, this method calculates "transposed eigenvectors": eigenvectors as they would be if PCA had been done this way:

* Original input matrix Y = X' (X transposed)
* Each column of Y has been centered (zero mean for each column)
* A covariance matrix has been calculated over rows (for pairs of rows) of Y
* Eigenvector decomposition has been done for this covariance matrix

Writes:

* scorefile.transformed.txt (a tab-delimited matrix: all "transposed eigenvectors" on rows)
