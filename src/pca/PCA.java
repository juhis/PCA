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
 * Calculates 1) eigenvector decomposition for the given symmetric matrix or 2)
 * principal component scores based on given eigenvector matrix and original
 * matrix.
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

        if ("cronbach".equals(args[0]) && args.length != 3) {
            printUsage();
            System.exit(1);
        }
        
        Path path1, path2;
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
                case "cronbach":
                    path1 = Paths.get(args[1]);
                    path2 = Paths.get(args[2]);
                    cronbach(path1, path2);
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
        System.out.println("java -jar PCA.jar evd symmetricmatrixfile");
        System.out.println("java -jar PCA.jar scores eigenvectorfile originalfile");
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
        log("Writing eigenvectors.txt");
        try (FileWriter fw = new FileWriter("eigenvectors.txt")) {
            FileUtil.writeMatrix(fw, rightEigenvectors, pcHeaders, headers, true);
        }

        log("Writing eigenvaluesReal.txt");
        try (FileWriter fw = new FileWriter("eigenvaluesReal.txt")) {
            FileUtil.writeArray(fw, "eigenvaluesReal", evd.getRealEigenvalues());
        }

        log("Writing eigenvaluesImaginary.txt");
        try (FileWriter fw = new FileWriter("eigenvaluesImaginary.txt")) {
            FileUtil.writeArray(fw, "eigenvaluesImaginary", evd.getImaginaryEigenvalues());
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
            double mean = 0;
            for (int c = 0; c < numCols; c++) {
                mean += originalMatrix.get(r, c) / numCols;
            }
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

        log("Writing scores.txt");
        try (FileWriter fw = new FileWriter("scores.txt")) {
            FileUtil.writeMatrix(fw, scoreMatrix, pcHeaders, headers);
        }

        log("Calculating Cronbach's alpha for each component");
        double[] alphas = cronbachsAlpha(evMatrix, scoreMatrix);
        log("Cronbach's alphas calculated");

        log("Writing cronbachsAlpha.txt");
        try (FileWriter fw = new FileWriter("cronbachsAlpha.txt")) {
            FileUtil.writeArray(fw, "cronbachsAlpha", alphas);
        }

        log("Done");
    }

    private static void cronbach(Path evPath, Path scorePath) throws IOException {
        
        log("Reading eigenvectors");
        int numRows = FileUtil.readColumnHeaders(evPath).length;
        log("Matrix size " + numRows + " x " + numRows);
        DenseMatrix evMatrix = new DenseMatrix(numRows, numRows);
        FileUtil.readMatrix(evPath, evMatrix);
        log("Eigenvectors read");

        log("Reading scores");
        int numCols = FileUtil.readColumnHeaders(scorePath).length;
        log("Matrix size " + numRows + " x " + numCols);
        DenseMatrix scoreMatrix = new DenseMatrix(numRows, numCols);
        FileUtil.readMatrix(scorePath, scoreMatrix);
        log("Scores read");
        
        log("Calculating Cronbach's alpha for each component");
        double[] alphas = cronbachsAlpha(evMatrix, scoreMatrix);
        log("Cronbach's alphas calculated");

        log("Writing cronbachsAlpha.txt");
        try (FileWriter fw = new FileWriter("cronbachsAlpha.txt")) {
            FileUtil.writeArray(fw, "cronbachsAlpha", alphas);
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
