package de.femtopedia.ldpc;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.function.BiFunction;
import java.util.function.IntBinaryOperator;
import java.util.stream.IntStream;
import lombok.Getter;

/**
 * Represents a Matrix with binary values, i.e. 0 or 1 in each slot.
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
        if (data.length <= 0 || data[0].length <= 0) {
            throw new IllegalArgumentException("Invalid matrix size");
        }

        rows = data.length;
        cols = data[0].length;
        this.data = data;
    }

    /**
     * Constructs the zero matrix of the given size.
     *
     * @param m The row size to initialize with.
     * @param n The column size to initialize with.
     * @return The identity matrix of size m x n.
     */
    public static BinaryMatrix zero(int m, int n) {
        return fromIntFunction(m, n, (i, j) -> 0);
    }

    /**
     * Constructs the identity matrix of the given size.
     *
     * @param n The size to initialize with.
     * @return The identity matrix of size n.
     */
    public static BinaryMatrix eye(int n) {
        return fromFunction(n, n, Integer::equals);
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
    public static BinaryMatrix fromFunction(int m, int n, BiFunction<Integer,
            Integer, Boolean> f) {
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
    public static BinaryMatrix fromIntFunction(int m, int n, BiFunction<Integer,
            Integer, Integer> f) {
        int[][] data = IntStream.range(0, m)
                .mapToObj(i -> IntStream.range(0, n)
                        .map(j -> f.apply(i, j)).toArray())
                .toArray(int[][]::new);
        return new BinaryMatrix(data);
    }

    /**
     * Generates a binary matrix from the given index-to-value functions.
     *
     * @param m  Amount of the elementary matrices per column.
     * @param n  Amount of the elementary matrices per row.
     * @param mm Amount of the elementary matrices' rows.
     * @param mn Amount of the elementary matrices' columns.
     * @param f  A function from which to generate the matrix from, mapping
     *           indices to values.
     * @return The binary matrix generated from the given index-to-value
     *         function.
     */
    public static BinaryMatrix fromMatrixFunction(int m, int n, int mm, int mn,
                                                  BiFunction<Integer, Integer,
                                                          BinaryMatrix> f) {
        return fromIntFunction(m * mm, n * mn,
                (i, j) -> f.apply(i / mm, j / mn)
                        .getEntryInt(i % mm, j % mn));
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
     * Shifts the matrix circularly to the right for the given amount of steps.
     *
     * @param shift The amount of shifts to apply.
     * @return The shifted matrix.
     */
    public BinaryMatrix shiftRight(int shift) {
        return BinaryMatrix.fromFunction(rows, cols, (i, j) ->
                data[i][(j - shift + cols) % cols]);
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
     * Adds (XOR) this matrix to another matrix.
     *
     * @param mat The matrix to add this matrix to.
     * @return The resulting binary matrix.
     */
    public BinaryMatrix add(BinaryMatrix mat) {
        int[][] newData = IntStream.range(0, rows).mapToObj(i ->
                IntStream.range(0, cols)
                        .map(j -> data[i][j] ^ mat.data[i][j] ? 1 : 0)
                        .toArray())
                .toArray(int[][]::new);
        return new BinaryMatrix(newData);
    }

    /**
     * Calculates and returns the (pseudo) lower triangular matrix.
     *
     * @return The (pseudo) lower triangular matrix.
     */
    public BinaryMatrix gaussJordan() {
        int k = getCols() - getRows();
        boolean[][] a = getColumns(0, k).getData();
        boolean[][] b = getColumns(k, getCols()).getData();

        gaussJordanIntern(b, a);
        return BinaryMatrix.horizConcat(new BinaryMatrix(a),
                new BinaryMatrix(b));
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

        gaussJordanIntern(newData, eye);
        return new BinaryMatrix(eye);
    }

    /**
     * Performs the Gauss-Jordan algorithm on the given matrices simultaneously.
     *
     * @param a The first matrix to perform the gaussian algorithm on.
     * @param b The second matrix to perform the gaussian algorithm on.
     */
    private void gaussJordanIntern(boolean[][] a, boolean[][] b) {
        // Gauss
        for (int i = 0; i < rows; i++) {
            ArrayUtils.parallelSortDesc(a, b);
            for (int j = i + 1; j < rows; j++) {
                if (gaussStep(a, b, i, j)) {
                    break;
                }
            }
        }

        // Gauss-Jordan
        for (int i = rows - 1; i > 0; i--) {
            for (int j = i - 1; j >= 0; j--) {
                gaussStep(a, b, i, j);
            }
        }
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
     * Returns whether this object is equal to the given object.
     *
     * @param o The Object to compare this Object to.
     * @return {@code true} iff this instance is equal to the given instance,
     *         {@code false} otherwise.
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof BinaryMatrix)) {
            return false;
        }
        BinaryMatrix that = (BinaryMatrix) o;
        return rows == that.rows
                && cols == that.cols
                && Arrays.deepEquals(data, that.data);
    }

    /**
     * Returns a Hash code representing this Object, respecting the hashCode()-
     * equals() contract.
     *
     * @return A Hash code representing this Object.
     */
    @Override
    public int hashCode() {
        int result = Objects.hash(rows, cols);
        result = 31 * result + Arrays.hashCode(data);
        return result;
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
