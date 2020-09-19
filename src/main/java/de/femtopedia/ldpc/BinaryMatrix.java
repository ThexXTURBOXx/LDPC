package de.femtopedia.ldpc;

import java.util.Arrays;
import java.util.function.IntBinaryOperator;
import java.util.stream.IntStream;

public class BinaryMatrix {

    private final boolean[][] data;
    private final int rows;
    private final int cols;

    public BinaryMatrix(int[][] data) {
        rows = data.length;
        cols = data[0].length;
        this.data = new boolean[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                this.data[i][j] = data[i][j] == 1;
            }
        }
    }

    public BinaryMatrix(boolean[][] data) {
        rows = data.length;
        cols = data[0].length;
        this.data = data;
    }

    public static BinaryMatrix eye(int n) {
        boolean[][] data = new boolean[n][n];
        for (int i = 0; i < n; i++) {
            data[i][i] = true;
        }
        return new BinaryMatrix(data);
    }

    public static BinaryMatrix concat(BinaryMatrix... matrices) {
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

    public int getRows() {
        return rows;
    }

    public int getCols() {
        return cols;
    }

    public boolean[][] getData() {
        return Arrays.stream(data).map(boolean[]::clone)
                .toArray(boolean[][]::new);
    }

    public boolean getEntry(int row, int col) {
        return data[row][col];
    }

    public BinaryMatrix transpose() {
        boolean[][] newData = new boolean[cols][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                newData[j][i] = data[i][j];
            }
        }
        return new BinaryMatrix(newData);
    }

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

    public int reduce(int identity, IntBinaryOperator op) {
        return Arrays.stream(data)
                .flatMapToInt(x -> IntStream.range(0, cols).map(i -> x[i] ? 1 : 0))
                .reduce(identity, op);
    }

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

    public BinaryMatrix inv() {
        boolean[][] newData = getData();
        boolean[][] eye = eye(rows).data;

        // Gauss
        for (int i = 0; i < rows; i++) {
            ArrayUtils.insertionSortParallel(newData, eye);
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

    private boolean gaussStep(boolean[][] newData, boolean[][] eye, int i, int j) {
        boolean[] nextRow = newData[j];
        if (nextRow[i]) {
            ArrayUtils.xor(nextRow, newData[i]);
            ArrayUtils.xor(eye[j], eye[i]);
            return false;
        }
        return true;
    }

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
