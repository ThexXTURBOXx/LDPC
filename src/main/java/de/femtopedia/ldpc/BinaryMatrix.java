package de.femtopedia.ldpc;

import java.util.Arrays;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.IntBinaryOperator;
import java.util.stream.IntStream;
import lombok.Getter;

/**
 * Represents a Matrix with binary values, i.e. 0 or 1 in each sot.
 */
public class BinaryMatrix {

    /**
     * The matrix's binary data.
     */
    private final boolean[][] data;

    /**
     * The amount of rows cached for better access.
     */
    @Getter
    private final int rows;

    /**
     * The amount of columns cached for better access.
     */
    @Getter
    private final int cols;

    /**
     * Initializes a binary matrix with the given data.
     *
     * @param data The data to initialize with as integers mod 2.
     */
    public BinaryMatrix(int[][] data) {
        rows = data.length;
        cols = data[0].length;
        this.data = new boolean[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                this.data[i][j] = (Math.abs(data[i][j]) % 2) == 1;
            }
        }
    }

    /**
     * Initializes a binary matrix with the given data.
     *
     * @param data The data to initialize with.
     */
    public BinaryMatrix(boolean[][] data) {
        rows = data.length;
        cols = data[0].length;
        this.data = data;
    }

    /**
     * Constructs the identity matrix of the given size.
     *
     * @param n The size to initialize with.
     * @return The identity matrix of size n.
     */
    public static BinaryMatrix eye(int n) {
        boolean[][] data = new boolean[n][n];
        for (int i = 0; i < n; i++) {
            data[i][i] = true;
        }
        return new BinaryMatrix(data);
    }

    /**
     * Generates a binary matrix from the given index-to-value functions.
     *
     * @param m Amount of the matrix's rows.
     * @param n Amount of the matrix's columns.
     * @param f A function from which to generate the matrix from, mapping
     *          indices to values.
     * @return The binary matrix generated from the given index-to-value
     *         function.
     */
    public static BinaryMatrix
    fromFunction(int m, int n, BiFunction<Integer, Integer, Boolean> f) {
        return fromIntFunction(m, n, (i, j) -> f.apply(i, j) ? 1 : 0);
    }

    /**
     * Generates a binary matrix from the given index-to-value functions.
     *
     * @param m Amount of the matrix's rows.
     * @param n Amount of the matrix's columns.
     * @param f A function from which to generate the matrix from, mapping
     *          indices to values.
     * @return The binary matrix generated from the given index-to-value
     *         function.
     */
    public static BinaryMatrix
    fromIntFunction(int m, int n, BiFunction<Integer, Integer, Integer> f) {
        int[][] data = IntStream.range(0, m)
                .mapToObj(i -> IntStream.range(0, n)
                        .map(j -> f.apply(i, j)).toArray())
                .toArray(int[][]::new);
        return new BinaryMatrix(data);
    }

    /**
     * Concatenates matrices' rows to one big matrix.
     *
     * @param matrices The matrices to concatenate.
     * @return The concatenated matrix.
     */
    public static BinaryMatrix horizConcat(BinaryMatrix... matrices) {
        int rows = matrices[0].rows;
        int cols = Arrays.stream(matrices).mapToInt(BinaryMatrix::getCols)
                .reduce(0, Integer::sum);
        boolean[][] newData = new boolean[rows][cols];

        for (int i = 0; i < rows; i++) {
            int j = 0;
            for (BinaryMatrix m : matrices) {
                for (int k = 0; k < m.cols; k++) {
                    newData[i][j++] = m.data[i][k];
                }
            }
        }
        return new BinaryMatrix(newData);
    }

    /**
     * Concatenates matrices' columns to one big matrix.
     *
     * @param matrices The matrices to concatenate.
     * @return The concatenated matrix.
     */
    public static BinaryMatrix vertConcat(BinaryMatrix... matrices) {
        int rows = Arrays.stream(matrices).mapToInt(BinaryMatrix::getRows)
                .reduce(0, Integer::sum);
        int cols = matrices[0].cols;
        boolean[][] newData = new boolean[rows][cols];

        int i = 0;
        for (BinaryMatrix m : matrices) {
            for (int j = 0; j < m.rows; j++) {
                newData[i++] = m.data[j].clone();
            }
        }
        return new BinaryMatrix(newData);
    }

    /**
     * Returns this matrix's binary data as booleans.
     *
     * @return This matrix's binary data as booleans.
     */
    public boolean[][] getData() {
        return Arrays.stream(data).map(boolean[]::clone)
                .toArray(boolean[][]::new);
    }

    /**
     * Returns this matrix's binary data as integers.
     *
     * @return This matrix's binary data as integers.
     */
    public int[][] getDataInt() {
        return Arrays.stream(data)
                .map(b -> IntStream.range(0, cols)
                        .map(i -> b[i] ? 1 : 0).toArray())
                .toArray(int[][]::new);
    }

    /**
     * Returns the entry at the given slot as boolean.
     *
     * @param row The row index.
     * @param col The column index.
     * @return The entry at the given slot as boolean.
     */
    public boolean getEntry(int row, int col) {
        return data[row][col];
    }

    /**
     * Returns the entry at the given slot as integer.
     *
     * @param row The row index.
     * @param col The column index.
     * @return The entry at the given slot as integer.
     */
    public int getEntryInt(int row, int col) {
        return data[row][col] ? 1 : 0;
    }

    /**
     * Permutes the matrix's columns according to the given permutation.
     *
     * @param permutation The column permutation to apply to this matrix.
     * @return A matrix representing this matrix with the given permutation
     *         applied.
     */
    public BinaryMatrix permuteColumns(List<Integer> permutation) {
        int[][] data = IntStream.range(0, rows).mapToObj(i ->
                IntStream.range(0, cols)
                        .map(j -> this.data[i][permutation.get(j)] ? 1 : 0)
                        .toArray())
                .toArray(int[][]::new);
        return new BinaryMatrix(data);
    }

    /**
     * Returns the transposed matrix.
     *
     * @return The transposed matrix.
     */
    public BinaryMatrix transpose() {
        boolean[][] newData = new boolean[cols][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                newData[j][i] = data[i][j];
            }
        }
        return new BinaryMatrix(newData);
    }

    /**
     * Returns the sub-matrix consisting of all columns with indices between the
     * given values.
     *
     * @param start The start column index.
     * @param end   The end column index.
     * @return The sub-matrix consisting of all columns with indices in the
     *         given range.
     */
    public BinaryMatrix getColumns(int start, int end) {
        boolean[][] newData = new boolean[rows][end - start];
        for (int i = 0; i < rows; i++) {
            if (end - start >= 0) {
                System.arraycopy(data[i],
                        start, newData[i],
                        0, end - start);
            }
        }
        return new BinaryMatrix(newData);
    }

    /**
     * Reduces all values in this matrix using the given integer operator.
     *
     * @param identity The identity for the given operator.
     * @param op       The operator.
     * @return The reduced value.
     */
    public int reduce(int identity, IntBinaryOperator op) {
        return Arrays.stream(data)
                .flatMapToInt(x -> IntStream.range(0, cols)
                        .map(i -> x[i] ? 1 : 0))
                .reduce(identity, op);
    }

    /**
     * Calculates and returns the sum of all entries in this matrix.
     *
     * @return The sum of all entries in this matrix.
     */
    public int sum() {
        return reduce(0, Integer::sum);
    }

    /**
     * Multiplies this matrix with another matrix.
     *
     * @param mat The matrix to multiply with.
     * @return The resulting binary matrix.
     */
    public BinaryMatrix mult(BinaryMatrix mat) {
        boolean[][] newData = new boolean[rows][mat.cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < mat.cols; j++) {
                for (int k = 0; k < cols; k++) {
                    newData[i][j] ^= data[i][k] && mat.data[k][j];
                }
            }
        }
        return new BinaryMatrix(newData);
    }

    /**
     * Returns whether this matrix is invertible, i.e. its determinant is != 0.
     *
     * @return {@code true}, iff the matrix is invertible, {@code false}
     *         otherwise.
     */
    public boolean isInvertible() {
        return det() != 0;
    }

    /**
     * Returns this matrix's determinant, calculated via gaussian algorithm.
     *
     * @return This matrix's determinant.
     */
    public int det() {
        boolean[][] newData = getData();

        // Gauss
        for (int i = 0; i < rows; i++) {
            ArrayUtils.sortDesc(newData);
            for (int j = i + 1; j < rows; j++) {
                if (gaussStep(newData, i, j)) {
                    break;
                }
            }
        }

        return Math.abs(IntStream.range(0, rows).map(i -> newData[i][i] ? 1 : 0)
                .reduce(1, (x, y) -> x * y));
    }

    /**
     * Calculates and returns the inverse matrix using Gauss-Jordan-Elimination.
     *
     * @return The inverse matrix.
     */
    public BinaryMatrix inv() {
        boolean[][] newData = getData();
        boolean[][] eye = eye(rows).data;

        // Gauss
        for (int i = 0; i < rows; i++) {
            ArrayUtils.parallelSortDesc(newData, eye);
            for (int j = i + 1; j < rows; j++) {
                if (gaussStep(newData, eye, i, j)) {
                    break;
                }
            }
        }

        // Gauss-Jordan
        for (int i = rows - 1; i > 0; i--) {
            for (int j = i - 1; j >= 0; j--) {
                gaussStep(newData, eye, i, j);
            }
        }
        return new BinaryMatrix(eye);
    }

    /**
     * Helper method performing one Gaussian step.
     *
     * @param data The matrix to perform the gaussian algorithm on.
     * @param i    The current comparing row index.
     * @param j    The current elimination row index.
     * @return {@code true} if no change was made, {@code false} otherwise.
     */
    private boolean gaussStep(boolean[][] data, int i, int j) {
        return gaussStep(data, null, i, j);
    }

    /**
     * Helper method performing one Gaussian step for 2 matrices simultaneously.
     *
     * @param data The first matrix to perform the gaussian algorithm on.
     * @param eye  The second matrix to perform the gaussian algorithm on.
     * @param i    The current comparing row index.
     * @param j    The current elimination row index.
     * @return {@code true} if no change was made, {@code false} otherwise.
     */
    private boolean gaussStep(boolean[][] data, boolean[][] eye, int i, int j) {
        boolean[] nextRow = data[j];
        if (nextRow[i]) {
            ArrayUtils.xor(nextRow, data[i]);
            if (eye != null) {
                ArrayUtils.xor(eye[j], eye[i]);
            }
            return false;
        }
        return true;
    }

    /**
     * Returns this matrix as a String consisting of 0s and 1s with spaces
     * between them.
     *
     * @return A String representing this matrix.
     */
    @Override
    public String toString() {
        String[] matStr = new String[rows];
        for (int i = 0; i < rows; i++) {
            String[] rowStr = new String[cols];
            for (int j = 0; j < cols; j++) {
                rowStr[j] = "" + (data[i][j] ? 1 : 0);
            }
            matStr[i] = String.join(" ", rowStr);
        }
        return String.join(System.lineSeparator(), matStr);
    }

}
