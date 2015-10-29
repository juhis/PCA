package pca;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.Date;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;

/**
 *
 * This program uses MTJ (Matrix Toolkits Java) which tries to utilize Java
 * Native Interface to run machine-optimized code. MTJ tells on runtime whether
 * using natives succeeds or not. It's a LOT faster with natives.
 *
 * https://github.com/fommil/matrix-toolkits-java/
 * http://lessthanoptimal.github.io/Java-Matrix-Benchmark/runtime/2013_10_Corei7v2600/
 *
 * 1) Eigenvector decomposition
 *
 * java -jar PCA.jar evd symmetricmatrixfile
 *
 * The first command line argument should be "evd". The second command line
 * argument should be a path to a tab-delimited text file containing a symmetric
 * matrix. The first row in the file, as well as the first column on each row,
 * should contain headers. The input matrix is assumed to be ok.
 *
 * Writes:
 *
 * eigenvectors.txt (a tab-delimited matrix: all right eigenvectors on rows)
 * eigenvaluesReal.txt (real parts of eigenvalues, one per line)
 * eigenvaluesImaginary.txt (imaginary parts of eigenvalues, one per line)
 *
 * 2) Principal component scores
 *
 * java -jar PCA.jar scores eigenvectorfile originalfile
 *
 * The first command line argument should be "scores". The second command line
 * argument should be a path to a tab-delimited text file containing
 * eigenvectors on rows. The first row in the file, as well as the first column
 * on each row, should contain headers.
 *
 * The third command line argument should be a path to a tab-delimited text file
 * containing the original data. The orientation is detected automatically so
 * that either rows or columns should correspond to the columns in the
 * eigenvector file. The first row in the file, as well as the first column on
 * each row, should contain headers.
 *
 * The input matrices are assumed to be ok.
 *
 * Writes:
 *
 * scores.txt (a tab-delimited matrix: all principal component scores on rows)
 * cronbachsAlpha.txt (Cronbach's alpha for each component, one per line)
 *
 * @author juha
 */
public class PCA {

    static final Format dateFormat = new SimpleDateFormat("yyyy-MM-dd");
    static final Format timeFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        if (args.length < 1) {
            printUsage();
            System.exit(1);
        }

        if ("evd".equals(args[0]) && args.length != 2) {
            printUsage();
            System.exit(1);
        }

        if ("scores".equals(args[0]) && args.length != 3) {
            printUsage();
            System.exit(1);
        }

        if ("transpose".equals(args[0]) && args.length != 2) {
            printUsage();
            System.exit(1);
        }

        if ("center".equals(args[0]) && args.length != 2) {
            printUsage();
            System.exit(1);
        }

        if ("scale".equals(args[0]) && args.length != 2) {
            printUsage();
            System.exit(1);
        }

        if ("covariance".equals(args[0]) && args.length != 2) {
            printUsage();
            System.exit(1);
        }

        if ("correlation".equals(args[0]) && args.length != 2) {
            printUsage();
            System.exit(1);
        }

        if ("transform".equals(args[0]) && args.length != 4) {
            printUsage();
            System.exit(1);
        }

        Path path1, path2, path3;
        try {
            switch (args[0]) {
                case "evd":
                    path1 = Paths.get(args[1]);
                    evd(path1);
                    break;
                case "scores":
                    path1 = Paths.get(args[1]);
                    path2 = Paths.get(args[2]);
                    scores(path1, path2);
                    break;
                case "transpose":
                    path1 = Paths.get(args[1]);
                    transpose(path1);
                    break;
                case "center":
                    path1 = Paths.get(args[1]);
                    center(path1);
                    break;
                case "scale":
                    path1 = Paths.get(args[1]);
                    scale(path1);
                    break;
                case "covariance":
                    path1 = Paths.get(args[1]);
                    covariance(path1);
                    break;
                case "correlation":
                    path1 = Paths.get(args[1]);
                    correlation(path1);
                    break;
                case "transform":
                    path1 = Paths.get(args[1]);
                    path2 = Paths.get(args[2]);
                    path3 = Paths.get(args[3]);
                    transform(path1, path2, path3);
                    break;
                default:
                    printUsage();
                    System.exit(1);
            }
        } catch (IOException ex) {
            System.err.println("Could not read or write file: " + ex.getMessage());
        } catch (NotConvergedException ex) {
            System.err.println("Eigenvector decomposition did not converge: " + ex.getMessage());
        }
    }

    private static void printUsage() {

        System.out.println("Usages:");
        System.out.println("java -jar PCA.jar center file");
        System.out.println("java -jar PCA.jar scale file");
        System.out.println("java -jar PCA.jar transpose file");
        System.out.println("java -jar PCA.jar covariance file");
        System.out.println("java -jar PCA.jar correlation file");
        System.out.println("java -jar PCA.jar evd symmetricmatrixfile");
        System.out.println("java -jar PCA.jar scores eigenvectorfile originalfile");
        System.out.println("java -jar PCA.jar transform eigenvectorfile eigenvaluefile originalfile");
    }

    private static void evd(Path path) throws IOException, NotConvergedException {

        log("Reading data");

        String[] headers = FileUtil.readColumnHeaders(path);
        int matrixSize = headers.length;
        log("Matrix size " + matrixSize + " x " + matrixSize);

        String[] pcHeaders = new String[matrixSize];
        for (int i = 0, len = pcHeaders.length; i < len; i++) {
            pcHeaders[i] = "PC" + (i + 1);
        }

        DenseMatrix matrix = new DenseMatrix(matrixSize, matrixSize);
        FileUtil.readMatrix(path, matrix);

        log("Data read");

        log("Calculating eigenvectors using MTJ");
        EVD evd = new EVD(matrixSize, false, true);
        evd.factor(matrix);
        log("Eigenvectors calculated");

        DenseMatrix rightEigenvectors = evd.getRightEigenvectors();
        log("Writing eigenvectors");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".eigenvectors.txt")) {
            FileUtil.writeMatrix(fw, rightEigenvectors, pcHeaders, headers, true);
        }

        log("Writing eigenvalues");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".eigenvalues.txt")) {
            FileUtil.writeArray(fw, "eigenvalues", evd.getRealEigenvalues());
        }

        log("Done");
    }

    private static void scores(Path evPath, Path originalPath) throws IOException {

        log("Reading eigenvectors");
        int numRows = FileUtil.readColumnHeaders(evPath).length;
        log("Matrix size " + numRows + " x " + numRows);
        DenseMatrix evMatrix = new DenseMatrix(numRows, numRows);
        FileUtil.readMatrix(evPath, evMatrix);
        log("Eigenvectors read");

        log("Reading original data");
        String[] headers = FileUtil.readColumnHeaders(originalPath);
        int numCols = headers.length;
        DenseMatrix originalMatrix;
        if (numCols == numRows) { // transpose data
            headers = FileUtil.readRowHeaders(originalPath);
            numCols = headers.length;
            log("Matrix size (transposed) " + numRows + " x " + numCols);
            originalMatrix = new DenseMatrix(numRows, numCols);
            FileUtil.readMatrix(originalPath, originalMatrix, true);
        } else {
            log("Matrix size " + numRows + " x " + numCols);
            originalMatrix = new DenseMatrix(numRows, numCols);
            FileUtil.readMatrix(originalPath, originalMatrix);
        }
        log("Original data read");

        log("Centering each row of original data");
        for (int r = 0; r < numRows; r++) {
            double mean = getMean(originalMatrix, r);
            for (int c = 0; c < numCols; c++) {
                originalMatrix.add(r, c, -mean);
            }
        }
        log("Rows centered");

        log("Calculating principal component scores");
        DenseMatrix scoreMatrix = new DenseMatrix(numRows, numCols);
        evMatrix.mult(originalMatrix, scoreMatrix);
        log("Principal component scores calculated");

        String[] pcHeaders = new String[numRows];
        for (int i = 0, len = pcHeaders.length; i < len; i++) {
            pcHeaders[i] = "PC" + (i + 1);
        }

        log("Writing scores");
        try (FileWriter fw = new FileWriter(originalPath.toString().replace(".gz", "").replace(".txt", "") + ".scores.txt")) {
            FileUtil.writeMatrix(fw, scoreMatrix, pcHeaders, headers);
        }

        log("Calculating Cronbach's alpha for each component");
        double[] alphas = cronbachsAlpha(evMatrix, scoreMatrix);
        log("Cronbach's alphas calculated");

        log("Writing Cronbach's alphas");
        try (FileWriter fw = new FileWriter(originalPath.toString().replace(".gz", "").replace(".txt", "") + ".cronbachsAlpha.txt")) {
            FileUtil.writeArray(fw, "cronbachsAlpha", alphas);
        }

        log("Done");
    }

    private static double getMean(DenseMatrix matrix, int row) {

        double mean = 0;
        for (int c = 0, numCols = matrix.numColumns(); c < numCols; c++) {
            mean += matrix.get(row, c) / numCols;
        }
        return mean;
    }

    private static double getVariance(DenseMatrix matrix, double mean, int row) {

        double variance = 0;
        for (int c = 0, numCols = matrix.numColumns(); c < numCols; c++) {
            variance += (matrix.get(row, c) - mean) * (matrix.get(row, c) - mean) / (numCols - 1);
        }
        return variance;
    }

    private static void center(Path path) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        log("Centering each row of data");
        for (int r = 0; r < numRows; r++) {
            double mean = getMean(matrix, r);
            for (int c = 0; c < numCols; c++) {
                matrix.add(r, c, -mean);
            }
        }
        log("Rows centered");

        log("Writing centered data");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".centered.txt")) {
            FileUtil.writeMatrix(fw, matrix, rowHeaders, colHeaders);
        }

        log("Done");
    }

    private static void scale(Path path) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        log("Scaling each row of data to a standard deviation of one");
        for (int r = 0; r < numRows; r++) {
            double mean = getMean(matrix, r);
            double variance = getVariance(matrix, mean, r);
            double stDev = Math.sqrt(variance);
            for (int c = 0; c < numCols; c++) {
                matrix.set(r, c, matrix.get(r, c) / stDev);
            }
        }
        log("Rows scaled");

        log("Writing scaled data");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".scaled.txt")) {
            FileUtil.writeMatrix(fw, matrix, rowHeaders, colHeaders);
        }

        log("Done");
    }

    private static void covariance(Path path) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        DenseMatrix covarianceMatrix = new DenseMatrix(numRows, numRows);
        log("Calculating covariance matrix");
        for (int r1 = 0; r1 < numRows; r1++) {
            double mean1 = getMean(matrix, r1);
            for (int r2 = r1; r2 < numRows; r2++) {
                double mean2 = getMean(matrix, r2);
                double covariance = 0;
                for (int c = 0; c < numCols; c++) {
                    covariance += (matrix.get(r1, c) - mean1) * (matrix.get(r2, c) - mean2) / (numCols - 1);
                }
                covarianceMatrix.set(r1, r2, covariance);
                covarianceMatrix.set(r2, r1, covariance);
            }
        }
        log("Covariance matrix calculated");

        log("Writing covariance matrix");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".covariance.txt")) {
            FileUtil.writeMatrix(fw, covarianceMatrix, rowHeaders, rowHeaders);
        }

        log("Done");
    }

    private static void correlation(Path path) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        DenseMatrix correlationMatrix = new DenseMatrix(numRows, numRows);
        log("Calculating correlation matrix");
        for (int r1 = 0; r1 < numRows; r1++) {
            correlationMatrix.set(r1, r1, 1);
            double mean1 = getMean(matrix, r1);
            double var1 = getVariance(matrix, mean1, r1);
            for (int r2 = r1 + 1; r2 < numRows; r2++) {
                double mean2 = getMean(matrix, r2);
                double var2 = getVariance(matrix, mean2, r2);
                double covariance = 0;
                for (int c = 0; c < numCols; c++) {
                    covariance += (matrix.get(r1, c) - mean1) * (matrix.get(r2, c) - mean2) / (numCols - 1);
                }
                double denom = Math.sqrt(var1 * var2);
                correlationMatrix.set(r1, r2, covariance / denom);
                correlationMatrix.set(r2, r1, covariance / denom);
            }
        }
        log("Correlation matrix calculated");

        log("Writing correlation matrix");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".correlation.txt")) {
            FileUtil.writeMatrix(fw, correlationMatrix, rowHeaders, rowHeaders);
        }

        log("Done");
    }

    private static void transpose(Path path) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        log("Writing transposed data");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".transposed.txt")) {
            FileUtil.writeMatrix(fw, matrix, colHeaders, rowHeaders, true);
        }

        log("Done");
    }

    private static void transform(Path evPath, Path eigenvaluePath, Path originalPath) throws IOException {

        log("Reading eigenvectors");
        int numRows = FileUtil.readColumnHeaders(evPath).length;
        log("Matrix size " + numRows + " x " + numRows);
        DenseMatrix evMatrix = new DenseMatrix(numRows, numRows);
        FileUtil.readMatrix(evPath, evMatrix);
        log("Eigenvectors read");

        log("Reading original data");
        String[] headers = FileUtil.readColumnHeaders(originalPath);
        int numCols = headers.length;
        DenseMatrix originalMatrix;
        if (numCols == numRows) { // transpose data
            headers = FileUtil.readRowHeaders(originalPath);
            numCols = headers.length;
            log("Matrix size (transposed) " + numRows + " x " + numCols);
            originalMatrix = new DenseMatrix(numRows, numCols);
            FileUtil.readMatrix(originalPath, originalMatrix, true);
        } else {
            log("Matrix size " + numRows + " x " + numCols);
            originalMatrix = new DenseMatrix(numRows, numCols);
            FileUtil.readMatrix(originalPath, originalMatrix);
        }
        log("Original data read");

        log("Reading eigenvalues");
        DenseMatrix eigenvalueMatrix = new DenseMatrix(numRows, 1);
        FileUtil.readArray(eigenvaluePath, eigenvalueMatrix, 0);
        log("Eigenvalues read");

        log("Calculating transformation");
        DenseMatrix transformedMatrix = new DenseMatrix(numRows, numCols);
        evMatrix.mult(originalMatrix, transformedMatrix);
        for (int comp = 0; comp < numRows; comp++) {
            for (int i = 0; i < numCols; i++) {
                double score = transformedMatrix.get(comp, i);
                transformedMatrix.set(comp, i, score / (Math.sqrt(2 * eigenvalueMatrix.get(comp, 0) * (numRows - 1) / (numCols - 1))));
            }
        }
        log("Transformation calculated");

        log("Writing transformed matrix");
        String[] pcHeaders = new String[numRows];
        for (int i = 0, len = pcHeaders.length; i < len; i++) {
            pcHeaders[i] = "PC" + (i + 1);
        }

        try (FileWriter fw = new FileWriter(originalPath.toString().replace(".gz", "").replace(".txt", "") + ".transformed.txt")) {
            FileUtil.writeMatrix(fw, transformedMatrix, pcHeaders, headers);
        }

        log("Done");
    }

    private static double[] cronbachsAlpha(DenseMatrix evMatrix, DenseMatrix scoreMatrix) {

        int numComps = evMatrix.numRows();
        int len = evMatrix.numColumns();
        int lenScores = scoreMatrix.numColumns();

        // Cronbach's alpha for each component
        double[] alphas = new double[numComps];
        for (int comp = 0; comp < numComps; comp++) {

            double evSquaredSum = 0;
            for (int i = 0; i < len; i++) {
                evSquaredSum += evMatrix.get(comp, i) * evMatrix.get(comp, i);
            }

            double scoreMean = 0;
            for (int i = 0; i < lenScores; i++) {
                scoreMean += scoreMatrix.get(comp, i) / lenScores;
            }

            double scoreVariance = 0;
            for (int i = 0; i < lenScores; i++) {
                double score = scoreMatrix.get(comp, i);
                scoreVariance += (score - scoreMean) * (score - scoreMean) / (lenScores - 1);
            }

            double alpha = (lenScores / (lenScores - 1d)) * (1d - (evSquaredSum / scoreVariance));
            alphas[comp] = alpha;
        }

        return alphas;
    }

    private static void log(String text) {

        String time = timeFormat.format(new Date());
        System.out.println(time + "\t" + text);
    }
}
