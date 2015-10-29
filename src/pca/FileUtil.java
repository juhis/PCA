package pca;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import no.uib.cipr.matrix.DenseMatrix;

/**
 *
 * @author juha
 */
public class FileUtil {

    private static int rowsInserted = 0;

    public static int countLines(Path path) throws IOException {

        int numLines;
        if (path.toString().endsWith(".gz")) {
            numLines = (int) GZipFiles.lines(path).count();
        } else {
            numLines = (int) Files.lines(path).count();
        }
        return numLines;
    }

    public static String[] readRowHeaders(Path path) throws IOException {

        List<String> headers = new ArrayList<>();
        if (path.toString().endsWith(".gz")) {
            GZipFiles.lines(path)
                    .skip(1) // header
                    .map(line -> line.split("\\t"))
                    .forEach(row -> headers.add(row[0]));
        } else {
            Files.lines(path)
                    .skip(1) // header
                    .map(line -> line.split("\\t"))
                    .forEach(row -> headers.add(row[0]));
        }
        return headers.toArray(new String[0]);
    }

    public static String[] readColumnHeaders(Path path) throws IOException {

        String[] headers;
        if (path.toString().endsWith(".gz")) {
            headers = GZipFiles.lines(path).findFirst().toString().split("\\t");
        } else {
            headers = Files.lines(path).findFirst().toString().split("\\t");
        }
        headers = Arrays.copyOfRange(headers, 1, headers.length);

        return headers;
    }

    public static void readMatrix(Path path, DenseMatrix matrix) throws IOException {

        FileUtil.readMatrix(path, matrix, false);
    }

    public static void readMatrix(Path path, DenseMatrix matrix, boolean transpose) throws IOException {

        rowsInserted = 0;
        if (path.toString().endsWith(".gz")) {
            GZipFiles.lines(path)
                    .skip(1) // header
                    .map(line -> line.split("\\t"))
                    .forEach(row -> insert(row, matrix, -1, transpose));
        } else {
            Files.lines(path)
                    .skip(1) // header
                    .map(line -> line.split("\\t"))
                    .forEach(row -> insert(row, matrix, -1, transpose));
        }
    }

    public static void readArray(Path path, DenseMatrix matrix, int column) throws IOException {
        readArray(path, matrix, column, false);
    }
    
    public static void readArray(Path path, DenseMatrix matrix, int column, boolean transpose) throws IOException {

        rowsInserted = 0;
        if (path.toString().endsWith(".gz")) {
            GZipFiles.lines(path)
                    .skip(1) // header
                    .map(line -> line.split("\\t"))
                    .forEach(row -> insert(row, matrix, column, transpose));
        } else {
            Files.lines(path)
                    .skip(1) // header
                    .map(line -> line.split("\\t"))
                    .forEach(row -> insert(row, matrix, column, transpose));
        }
    }

    private static void insert(String[] row, DenseMatrix matrix, int column, boolean transpose) {

        if (column > -1) {
            if (transpose) {
                matrix.set(0, rowsInserted, Double.parseDouble(row[column]));
            } else {
                matrix.set(rowsInserted, 0, Double.parseDouble(row[column]));
            }
        } else {
            for (int i = 1, len = row.length; i < len; i++) {
                if (transpose) {
                    matrix.set(i - 1, rowsInserted, Double.parseDouble(row[i]));
                } else {
                    matrix.set(rowsInserted, i - 1, Double.parseDouble(row[i]));
                }
            }
        }
        ++rowsInserted;
    }

    public static void writeMatrix(FileWriter fw, DenseMatrix matrix, String[] rowHeaders, String[] colHeaders) throws IOException {

        writeMatrix(fw, matrix, rowHeaders, colHeaders, false);
    }

    public static void writeMatrix(FileWriter fw, DenseMatrix matrix, String[] rowHeaders, String[] colHeaders, boolean transpose) throws IOException {

        String date = PCA.dateFormat.format(new Date());

        // write date and column headers to first line
        fw.write(date);
        if (transpose) {
            for (int c = 0, cols = matrix.numRows(); c < cols; c++) {
                fw.write("\t" + colHeaders[c]);
            }
        } else {
            for (int c = 0, cols = matrix.numColumns(); c < cols; c++) {
                fw.write("\t" + colHeaders[c]);
            }
        }
        fw.write("\n");

        // write data with row headers
        if (transpose) {
            for (int c = 0, cols = matrix.numColumns(); c < cols; c++) {
                fw.write(rowHeaders[c]);
                for (int r = 0, rows = matrix.numRows(); r < rows; r++) {
                    fw.write("\t" + matrix.get(r, c));
                }
                fw.write("\n");
            }
        } else {
            for (int r = 0, rows = matrix.numRows(); r < rows; r++) {
                fw.write(rowHeaders[r]);
                for (int c = 0, cols = matrix.numColumns(); c < cols; c++) {
                    fw.write("\t" + matrix.get(r, c));
                }
                fw.write("\n");
            }
        }
    }

    public static void writeArray(FileWriter fw, String name, double[] array) throws IOException {

        fw.write(name + "\n");
        for (int i = 0, len = array.length; i < len; i++) {
            fw.write(array[i] + "\n");
        }
    }

}
