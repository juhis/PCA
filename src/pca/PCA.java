package pca;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;

public class PCA {

    static final Format dateFormat = new SimpleDateFormat("yyyy-MM-dd");
    static final Format timeFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

    /**
     *
     * Run one of the defined matrix operations:
     *
     * center, scale, transpose, covariance, correlation, evd, scores, transform
     *
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

        if ("scores".equals(args[0]) && args.length < 3) {
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

        if ("transform".equals(args[0]) && args.length != 3) {
            printUsage();
            System.exit(1);
        }

        if ("noise".equals(args[0]) && args.length != 2) {
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
                    boolean isEVTransposed = false;
                    if (args.length > 3 && "transpose".equals(args[3])) {
                        isEVTransposed = true;
                    }
                    scores(path1, path2, isEVTransposed);
                    break;
                case "transpose":
                    path1 = Paths.get(args[1]);
                    transpose(path1);
                    break;
                case "center-rows":
                    path1 = Paths.get(args[1]);
                    center(path1, false);
                    break;
                case "center-columns":
                    path1 = Paths.get(args[1]);
                    center(path1, true);
                    break;
                case "scale-rows":
                    path1 = Paths.get(args[1]);
                    scale(path1, false);
                    break;
                case "scale-columns":
                    path1 = Paths.get(args[1]);
                    scale(path1, true);
                    break;
                case "covariance-rows":
                    path1 = Paths.get(args[1]);
                    covcor(path1, false, false);
                    break;
                case "covariance-columns":
                    path1 = Paths.get(args[1]);
                    covcor(path1, false, true);
                    break;
                case "correlation-rows":
                    path1 = Paths.get(args[1]);
                    covcor(path1, true, false);
                    break;
                case "correlation-columns":
                    path1 = Paths.get(args[1]);
                    covcor(path1, true, true);
                    break;
                case "transform":
                    path1 = Paths.get(args[1]);
                    path2 = Paths.get(args[2]);
                    transform(path1, path2);
                    break;
                case "splithalf":
                    path1 = Paths.get(args[1]);
                    path2 = Paths.get(args[2]);
                    isEVTransposed = false;
                    if (args.length > 3 && "transpose".equals(args[3])) {
                        isEVTransposed = true;
                    }
                    splitHalf(path1, path2, isEVTransposed);
                    break;
                case "cronbach":
                    path1 = Paths.get(args[1]);
                    path2 = Paths.get(args[2]);
                    path3 = Paths.get(args[3]);
                    isEVTransposed = false;
                    if (args.length > 4 && "transpose".equals(args[4])) {
                        isEVTransposed = true;
                    }
                    cronbach(path1, path2, path3, isEVTransposed);
                    break;
                case "noise":
                    path1 = Paths.get(args[1]);
                    addNoise(path1);
                    break;
                default:
                    printUsage();
                    System.exit(1);
            }
        } catch (IOException ex) {
            System.err.println("Could not read or write file: " + ex.getMessage());
        } catch (NotConvergedException ex) {
            System.err.println("Eigenvector decomposition did not converge: " + ex.getMessage());
        } catch (InterruptedException ex) {
            System.err.println("Interrupted:" + ex.getMessage());
        }
    }

    /**
     *
     * Prints program usage instructions to stdout.
     *
     */
    private static void printUsage() {

        System.out.println("Usages:");
        System.out.println("java -jar PCA.jar center-rows file");
        System.out.println("java -jar PCA.jar center-columns file");
        System.out.println("java -jar PCA.jar scale-rows file");
        System.out.println("java -jar PCA.jar scale-columns file");
        System.out.println("java -jar PCA.jar transpose file");
        System.out.println("java -jar PCA.jar covariance-rows file");
        System.out.println("java -jar PCA.jar covariance-columns file");
        System.out.println("java -jar PCA.jar correlation-rows file");
        System.out.println("java -jar PCA.jar correlation-columns file");
        System.out.println("java -jar PCA.jar evd symmetricmatrixfile");
        System.out.println("java -jar PCA.jar scores eigenvectorfile originalfile");
        System.out.println("java -jar PCA.jar transform scorefile eigenvaluefile");
    }

    /**
     *
     * Eigenvector decomposition.
     *
     * Calculates right eigenvectors and (real parts of) eigenvalues for a given
     * matrix.
     *
     * The given path has to contain a symmetric matrix. Uses natives and is
     * multi-threaded when natives are available. Single-threaded when natives
     * are not available.
     *
     * Writes two files: path.eigenvectors.txt and path.eigenvalues.txt
     *
     * @param path Path of matrix file to run decomposition on
     * @throws IOException If cannot read from / write to disk
     * @throws NotConvergedException If eigenvector decomposition doesn't
     * converge
     */
    public static void evd(Path path) throws IOException, NotConvergedException {

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

    /**
     *
     * Principal component scores.
     *
     * Calculates principal component scores and Cronbach's alpha values based
     * on a given eigenvector matrix and original data matrix.
     *
     * Orientation of the original matrix is automatically detected.
     *
     * Writes two files: originalPath.scores.txt and
     * originalPath.cronbachsAlpha.txt
     *
     * @param evPath Path to an eigenvector matrix where each row is an
     * eigenvector
     * @param originalPath Path to the original data matrix
     * @param isEVTransposed True if each column is an eigenvector in the
     * eigenvector matrix
     * @throws IOException If cannot read from / write to disk
     */
    public static void scores(Path evPath, Path originalPath, boolean isEVTransposed) throws IOException {

        if (isEVTransposed) {
            log("Reading eigenvectors (transposed)");
        } else {
            log("Reading eigenvectors");
        }
        int numRows = FileUtil.readColumnHeaders(evPath).length;
        log("Matrix size " + numRows + " x " + numRows);
        DenseMatrix evMatrix = new DenseMatrix(numRows, numRows);
        FileUtil.readMatrix(evPath, evMatrix, isEVTransposed);
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
            double mean = getRowMean(originalMatrix, r);
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
        double[] alphas = cronbachsAlpha(originalMatrix, evMatrix, scoreMatrix);
        log("Cronbach's alphas calculated");

        log("Writing Cronbach's alphas");
        try (FileWriter fw = new FileWriter(originalPath.toString().replace(".gz", "").replace(".txt", "") + ".cronbachsAlpha.txt")) {
            FileUtil.writeArray(fw, "cronbachsAlpha", alphas);
        }

        log("Done");
    }

    /**
     *
     * Calculates mean of a given row.
     *
     * @param matrix Data matrix
     * @param row Row index
     * @return Arithmetic mean of the row
     */
    private static double getRowMean(DenseMatrix matrix, int row) {

        double mean = 0;
        for (int c = 0, numCols = matrix.numColumns(); c < numCols; c++) {
            mean += matrix.get(row, c) / numCols;
        }
        return mean;
    }

    /**
     *
     * Calculates mean of a given column.
     *
     * @param matrix Data matrix
     * @param col Column index
     * @return Arithmetic mean of the column
     */
    private static double getColumnMean(DenseMatrix matrix, int col) {

        double mean = 0;
        for (int r = 0, numRows = matrix.numRows(); r < numRows; r++) {
            mean += matrix.get(r, col) / numRows;
        }
        return mean;
    }

    /**
     *
     * Calculates variance for a given row.
     *
     * @param matrix Data matrix
     * @param mean Mean of the row
     * @param row Row index
     * @return Sample variance of the row
     */
    private static double getRowVariance(DenseMatrix matrix, double mean, int row) {

        double variance = 0;
        for (int c = 0, numCols = matrix.numColumns(); c < numCols; c++) {
            variance += (matrix.get(row, c) - mean) * (matrix.get(row, c) - mean) / (numCols - 1);
        }
        return variance;
    }

    /**
     *
     * Calculates variance for a given row.
     *
     * @param matrix Data matrix
     * @param row Row index
     * @return Sample variance of the row
     */
    private static double getRowVariance(DenseMatrix originalMatrix, int row) {

        double mean = getRowMean(originalMatrix, row);
        return getRowVariance(originalMatrix, mean, row);
    }
    
    /**
     *
     * Calculates variance for a given column.
     *
     * @param matrix Data matrix
     * @param mean Mean of the column
     * @param col Column index
     * @return Sample variance of the column
     */
    private static double getColumnVariance(DenseMatrix matrix, double mean, int col) {

        double variance = 0;
        for (int r = 0, numRows = matrix.numRows(); r < numRows; r++) {
            variance += (matrix.get(r, col) - mean) * (matrix.get(r, col) - mean) / (numRows - 1);
        }
        return variance;
    }

    /**
     *
     * Center each row or column in a given matrix.
     *
     * For each element, subtracts the mean of the corresponding row or column,
     * thus centering each row or column to a mean of zero.
     *
     * Writes path.means.txt and path.centered.txt
     *
     * @param path Path to the data matrix
     * @param onColumns true if columns are to be centered instead of rows
     * @throws IOException If cannot read from / write to disk
     */
    public static void center(Path path, boolean onColumns) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        String suffix;
        double[] means;
        if (onColumns) {
            log("Centering each column of data");
            suffix = "columns";
            means = new double[numCols];
            for (int c = 0; c < numCols; c++) {
                double mean = getColumnMean(matrix, c);
                means[c] = mean;
                for (int r = 0; r < numRows; r++) {
                    matrix.add(r, c, -mean);
                }
            }
            log("Columns centered");
        } else {
            log("Centering each row of data");
            suffix = "rows";
            means = new double[numRows];
            for (int r = 0; r < numRows; r++) {
                double mean = getRowMean(matrix, r);
                means[r] = mean;
                for (int c = 0; c < numCols; c++) {
                    matrix.add(r, c, -mean);
                }
            }
            log("Rows centered");
        }

        log("Writing centered data");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".means-" + suffix + ".txt")) {
            FileUtil.writeArray(fw, "means-" + suffix, means);
        }
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".centered-" + suffix + ".txt")) {
            FileUtil.writeMatrix(fw, matrix, rowHeaders, colHeaders);
        }

        log("Done");
    }

    /**
     *
     * Scale each row or column in a given matrix.
     *
     * Divides each element with the standard deviation of the corresponding
     * row/column, thus scaling each row/column to a standard deviation of one.
     *
     * Writes path stdevs.txt and path.scaled.txt
     *
     * @param path Path to the data matrix
     * @param onColumns true if columns are to be centered instead of rows
     * @throws IOException If cannot read from / write to disk
     */
    public static void scale(Path path, boolean onColumns) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        double[] stDevs;
        String suffix;
        if (onColumns) {
            log("Scaling each column of data to a standard deviation of one");
            suffix = "columns";
            stDevs = new double[numCols];
            for (int c = 0; c < numCols; c++) {
                double mean = getColumnMean(matrix, c);
                double variance = getColumnVariance(matrix, mean, c);
                double stDev = Math.sqrt(variance);
                stDevs[c] = stDev;
                for (int r = 0; r < numRows; r++) {
                    matrix.set(r, c, matrix.get(r, c) / stDev);
                }
            }
            log("Columns scaled");
        } else {
            log("Scaling each row of data to a standard deviation of one");
            suffix = "rows";
            stDevs = new double[numRows];
            for (int r = 0; r < numRows; r++) {
                double mean = getRowMean(matrix, r);
                double variance = getRowVariance(matrix, mean, r);
                double stDev = Math.sqrt(variance);
                stDevs[r] = stDev;
                for (int c = 0; c < numCols; c++) {
                    matrix.set(r, c, matrix.get(r, c) / stDev);
                }
            }
            log("Rows scaled");
        }
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".stdevs-" + suffix + ".txt")) {
            FileUtil.writeArray(fw, "stdevs-" + suffix, stDevs);
        }

        log("Writing scaled data");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".scaled-" + suffix + ".txt")) {
            FileUtil.writeMatrix(fw, matrix, rowHeaders, colHeaders);
        }

        log("Done");
    }

    /**
     *
     * Covariance or Pearson correlation over rows or columns.
     *
     * Calculates covariance or Pearson correlation for each pair of
     * rows/columns in a given matrix.
     *
     * Writes the symmetric covariance matrix path.covariance.txt or correlation
     * matrix path.correlation.txt
     *
     * @param path Path to the data matrix
     * @param doCorrelation true if correlation is used and not covariance
     * @param onColumns true if columns are to be centered instead of rows
     * @throws IOException If cannot read from / write to disk
     */
    public static void covcor(Path path, boolean doCorrelation, boolean onColumns) throws IOException, InterruptedException {

        log("Reading data");
        String[] colHeaders;
        String[] rowHeaders;
        String suffix;
        if (onColumns) {
            suffix = "columns";
            colHeaders = FileUtil.readRowHeaders(path);
            rowHeaders = FileUtil.readColumnHeaders(path);
        } else {
            suffix = "rows";
            colHeaders = FileUtil.readColumnHeaders(path);
            rowHeaders = FileUtil.readRowHeaders(path);
        }
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix, onColumns);
        log("Data read");

        DenseMatrix corcovMatrix = new DenseMatrix(numRows, numRows);
        String type = doCorrelation ? "correlation" : "covariance";
        log("Calculating " + type + " matrix over " + suffix);
        double[] means = new double[numRows];
        for (int r = 0; r < numRows; r++) {
            means[r] = getRowMean(matrix, r);
        }
        double[] variances = new double[numRows];
        if (doCorrelation) {
            for (int r = 0; r < numRows; r++) {
                variances[r] = getRowVariance(matrix, means[r], r);
            }
        }

        double[] data = matrix.getData();
        double[] covcorData = corcovMatrix.getData();
        ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        class Task implements Callable {

            int i1, i2;

            Task(int r1, int r2) {
                i1 = r1;
                i2 = r2;
            }

            @Override
            public Void call() {
                for (int r1 = i1; r1 < Math.min(i2, numRows); r1++) {
                    covcorData[r1 * numRows + r1] = 1;
                    for (int r2 = r1; r2 < numRows; r2++) {
                        double covariance = 0;
                        for (int c = 0; c < numCols; c++) {
                            covariance += (data[c * numRows + r1] - means[r1]) * (data[c * numRows + r2] - means[r2]) / (numCols - 1);
                        }
                        if (doCorrelation) {
                            double denom = Math.sqrt(variances[r1] * variances[r2]);
                            covcorData[r1 * numRows + r2] = covariance / denom;
                            covcorData[r2 * numRows + r1] = covariance / denom;
                        } else {
                            covcorData[r1 * numRows + r2] = covariance;
                            covcorData[r2 * numRows + r1] = covariance;
                        }
                    }
                    if (r1 > 0 && r1 % 1000 == 0) {
                        log(r1 + " rows processed");
                    }
                }
                return null;
            }
        }

        List<Callable<Object>> todo = new ArrayList<>();
        for (int row = 0; row < numRows; row += 100) {
            todo.add(new Task(row, row + 100));
        }
        executor.invokeAll(todo);
        executor.shutdown();
        log("Matrix calculated");

        log("Writing " + type + " matrix");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + "." + type + "-" + suffix + ".txt")) {
            FileUtil.writeMatrix(fw, corcovMatrix, rowHeaders, rowHeaders);
        }

        log("Done");
    }

    /**
     *
     * Transpose a given data matrix.
     *
     * Writes path.transposed.txt
     *
     * @param path Path to the data matrix
     * @throws IOException If cannot read from / write to disk
     */
    public static void transpose(Path path) throws IOException {

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

    /**
     *
     * Change orientation of PCA results.
     *
     * Calculates "transposed eigenvectors" from a given principal component
     * score file and eigenvalue file.
     *
     * Note: This works in the following scenario:
     *
     * - There is an original input matrix X
     *
     * - Each column of X has been centered (zero mean for each column)
     *
     * - A covariance matrix has been calculated over rows (for pairs of rows)
     * of X
     *
     * - Eigenvector decomposition has been done for this covariance matrix
     *
     * - Principal component scores have been calculated based on the resulting
     * eigenvectors
     *
     * Then, this method calculates "transposed eigenvectors": eigenvectors as
     * they would be if PCA had been done this way:
     *
     * - Original input matrix Y = X' (X transposed)
     *
     * - Each column of Y has been centered (zero mean for each column)
     *
     * - A covariance matrix has been calculated over rows (for pairs of rows)
     * of Y
     *
     * - Eigenvector decomposition has been done for this covariance matrix
     *
     * Writes scorePath.transformedToEigenvectors.txt
     *
     * @param scorePath Path to a principal component score matrix where each
     * row is a component
     * @param eigenvaluePath Path to a file with eigenvalues (corresponding to
     * the eigenvectors from which the principal component scores have been
     * calculated)
     * @throws IOException If cannot read from / write to disk
     */
    public static void transform(Path scorePath, Path eigenvaluePath) throws IOException {

        log("Reading scores");
        String[] colHeaders = FileUtil.readColumnHeaders(scorePath);
        String[] rowHeaders = FileUtil.readRowHeaders(scorePath);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix scoreMatrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(scorePath, scoreMatrix);
        log("Scores read");

        log("Reading eigenvalues");
        DenseMatrix eigenvalueMatrix = new DenseMatrix(numRows, 1);
        FileUtil.readArray(eigenvaluePath, eigenvalueMatrix, 0);
        log("Eigenvalues read");

        log("Calculating transposed eigenvectors");
        DenseMatrix transposedEVMatrix = new DenseMatrix(numRows, numCols);
        for (int comp = 0; comp < numRows; comp++) {
            for (int i = 0; i < numCols; i++) {
                double score = scoreMatrix.get(comp, i);
                transposedEVMatrix.set(comp, i, score / (Math.sqrt((numCols) * eigenvalueMatrix.get(comp, 0))));
            }
        }
        log("Transformation calculated");

        try (FileWriter fw = new FileWriter(scorePath.toString().replace(".gz", "").replace(".txt", "") + ".transformedToEigenvectors.txt")) {
            FileUtil.writeMatrix(fw, transposedEVMatrix, rowHeaders, colHeaders);
        }

        log("Done");
    }

    /**
     *
     * Calculate Cronbach's alpha for each principal component.
     *
     * @param evMatrix Eigenvector matrix, each row is an eigenvector
     * @param scoreMatrix Principal component score matrix, each row is a
     * component
     * @return An array of Cronbach's alpha values
     */
    private static double[] cronbachsAlpha(DenseMatrix originalMatrix, DenseMatrix evMatrix, DenseMatrix scoreMatrix) {

        int numComps = evMatrix.numRows();
        int lenItems = evMatrix.numColumns();
        int lenScores = scoreMatrix.numColumns();

        double[] variances = new double[lenItems];
        for (int i = 0; i < lenItems; i++) {
            variances[i] = getRowVariance(originalMatrix, i);
        }

        double[] alphas = new double[numComps];
        for (int comp = 0; comp < numComps; comp++) {

            double sumVariance = 0;
            for (int i = 0; i < lenItems; i++) {
                sumVariance += evMatrix.get(i, comp) * evMatrix.get(i, comp) * variances[i];
            }

            double scoreMean = getRowMean(scoreMatrix, comp);
            double scoreVariance = getRowVariance(scoreMatrix, scoreMean, comp);

            double alpha = (lenItems / (lenItems - 1d)) * (1d - (sumVariance / scoreVariance));
            alphas[comp] = alpha;
        }

        return alphas;
    }

    /**
     *
     * Calculate split-half correlations for each principal component.
     *
     * @param evPath Eigenvector matrix, each row is an eigenvector
     * @param originalPath Path to the original data matrix
     * @param isEVTransposed True if each column is an eigenvector in the
     * eigenvector matrix
     * @throws IOException If cannot read from / write to disk
     */
    private static void splitHalf(Path evPath, Path originalPath, boolean isEVTransposed) throws IOException {

        if (isEVTransposed) {
            log("Reading eigenvectors (transposed)");
        } else {
            log("Reading eigenvectors");
        }
        int numRows = FileUtil.readColumnHeaders(evPath).length;
        log("Matrix size " + numRows + " x " + numRows);
        DenseMatrix evMatrix = new DenseMatrix(numRows, numRows);
        FileUtil.readMatrix(evPath, evMatrix, isEVTransposed);
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
            double mean = getRowMean(originalMatrix, r);
            for (int c = 0; c < numCols; c++) {
                originalMatrix.add(r, c, -mean);
            }
        }
        log("Rows centered");

        List<Integer> indexList = new ArrayList<>();
        for (int r = 0; r < numRows; r++) {
            indexList.add(r);
        }
        Collections.shuffle(indexList);
        int[] splitIndices1 = new int[numRows / 2];
        int[] splitIndices2 = new int[numRows / 2];
        for (int r = 0; r < numRows / 2; r++) {
            splitIndices1[r] = indexList.get(r);
            splitIndices2[r] = indexList.get(r + numRows / 2);
        }

        log("Calculating split half correlations");
        for (int comp = 0; comp < numRows; comp++) {

            double[] pcScoresThisCompSplit1 = new double[numCols];
            double[] pcScoresThisCompSplit2 = new double[numCols];
            for (int c = 0; c < numCols; c++) {
                for (int i = 0; i < numRows / 2; i++) {
                    pcScoresThisCompSplit1[c] += originalMatrix.get(splitIndices1[i], c) * evMatrix.get(comp, splitIndices1[i]);
                    pcScoresThisCompSplit2[c] += originalMatrix.get(splitIndices2[i], c) * evMatrix.get(comp, splitIndices2[i]);
                }
            }

            double mean1 = JSci.maths.ArrayMath.mean(pcScoresThisCompSplit1);
            double var1 = JSci.maths.ArrayMath.variance(pcScoresThisCompSplit1);
            double mean2 = JSci.maths.ArrayMath.mean(pcScoresThisCompSplit2);
            double var2 = JSci.maths.ArrayMath.variance(pcScoresThisCompSplit2);
            double covariance = 0;
            for (int c = 0; c < numCols; c++) {
                covariance += (pcScoresThisCompSplit1[c] - mean1) * (pcScoresThisCompSplit2[c] - mean2) / (numCols - 1);
            }
            double denom = Math.sqrt(var1 * var2);

            log((comp + 1) + "\t" + covariance / denom);
        }
        log("Split half correlations calculated");
    }

    /**
     *
     * Print given text to stdout, preceded by a MySQL-style timestamp and a
     * tab.
     *
     * @param text Text to log
     */
    private static void log(String text) {

        String time = timeFormat.format(new Date());
        System.out.println(time + "\t" + text);
    }

    private static void cronbach(Path originalPath, Path evPath, Path scorePath, boolean isEVTransposed) throws IOException {

        if (isEVTransposed) {
            log("Reading eigenvectors (transposed)");
        } else {
            log("Reading eigenvectors");
        }
        int numRows = FileUtil.readColumnHeaders(evPath).length;
        log("Matrix size " + numRows + " x " + numRows);
        DenseMatrix evMatrix = new DenseMatrix(numRows, numRows);
        FileUtil.readMatrix(evPath, evMatrix, isEVTransposed);
        log("Eigenvectors read");

        log("Reading scores");
        String[] colHeaders = FileUtil.readColumnHeaders(scorePath);
        int numCols = colHeaders.length;
        DenseMatrix scoreMatrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(scorePath, scoreMatrix);
        log("Scores read");

        log("Reading original data");
        String[] headers = FileUtil.readColumnHeaders(originalPath);
        numCols = headers.length;
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

        log("Calculating Cronbach's alpha for each component");
        double[] alphas = cronbachsAlpha(originalMatrix, evMatrix, scoreMatrix);
        log("Cronbach's alphas calculated");

        log("Writing Cronbach's alphas");
        try (FileWriter fw = new FileWriter(scorePath.toString().replace(".gz", "").replace(".txt", "") + ".cronbachsAlpha.txt")) {
            FileUtil.writeArray(fw, "cronbachsAlpha", alphas);
        }
    }

    private static void addNoise(Path path) throws IOException {

        log("Reading data");
        String[] colHeaders = FileUtil.readColumnHeaders(path);
        String[] rowHeaders = FileUtil.readRowHeaders(path);
        int numCols = colHeaders.length;
        int numRows = rowHeaders.length;
        DenseMatrix matrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(path, matrix);
        log("Data read");

        log("Adding noise: a random value from the uniform distribution [0, 0.001[");
        for (int r = 0; r < numRows; r++) {
            for (int c = 0; c < numCols; c++) {
                matrix.add(r, c, Math.random() / 1000);
            }
        }
        log("Noise added");

        log("Writing data");
        try (FileWriter fw = new FileWriter(path.toString().replace(".gz", "").replace(".txt", "") + ".noisied.txt")) {
            FileUtil.writeMatrix(fw, matrix, rowHeaders, colHeaders);
        }

        log("Done");
    }
}
