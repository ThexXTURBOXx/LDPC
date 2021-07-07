package de.femtopedia.ldpc.util;

import de.femtopedia.ldpc.LDPC;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.bouncycastle.pqc.math.linearalgebra.GF2Matrix;
import org.bouncycastle.pqc.math.linearalgebra.GF2Vector;

/**
 * Helper class providing access to some operations on {@link GF2Vector}s and
 * {@link GF2Matrix}s.
 */
public final class MatrixUtils {

    /**
     * Don't initialize me.
     */
    private MatrixUtils() {
        throw new UnsupportedOperationException();
    }

    /**
     * Converts the given boolean array to a {@link GF2Vector}.
     *
     * @param data The boolean array to convert.
     * @return The converted {@link GF2Vector}.
     */
    public static GF2Vector toVec(boolean[] data) {
        GF2Vector vec = new GF2Vector(data.length);
        IntStream.range(0, data.length)
                .filter(i -> data[i])
                .forEach(vec::setBit);
        return vec;
    }

    /**
     * Converts the given int array to a {@link GF2Vector}.
     *
     * @param data The int array to convert.
     * @return The converted {@link GF2Vector}.
     */
    public static GF2Vector toVec(int[] data) {
        GF2Vector vec = new GF2Vector(data.length);
        IntStream.range(0, data.length)
                .filter(i -> data[i] != 0)
                .forEach(vec::setBit);
        return vec;
    }

    /**
     * Adds some noise to the given {@link GF2Vector} by flipping the bits at
     * the given indices.
     *
     * @param vec     The {@link GF2Vector} to add noise to.
     * @param indices The indices to add bitflips at.
     * @return The given {@link GF2Vector} with some noise.
     */
    public static GF2Vector addNoise(GF2Vector vec, Integer... indices) {
        return addNoise(vec, Arrays.stream(indices)
                .collect(Collectors.toSet()));
    }

    /**
     * Adds some noise to the given {@link GF2Vector} by flipping the bits at
     * the given indices.
     *
     * @param vec     The {@link GF2Vector} to add noise to.
     * @param indices The indices to add bitflips at.
     * @return The given {@link GF2Vector} with some noise.
     */
    public static GF2Vector addNoise(GF2Vector vec, Set<Integer> indices) {
        GF2Vector newVec = new GF2Vector(vec.getLength());
        for (int i = 0; i < vec.getLength(); i++) {
            int x = vec.getBit(i);
            if (indices.remove(i)) {
                x = x == 0 ? 1 : 0;
            }
            if (x != 0) {
                newVec.setBit(i);
            }
        }
        return newVec;
    }

    /**
     * Returns the boolean entry from the given {@link GF2Vector} at the given
     * position.
     *
     * @param vec The {@link GF2Vector} to return the entry from.
     * @param col The index to retrieve the entry from.
     * @return The entry at the given position in the given {@link GF2Vector}.
     */
    public static boolean getEntry(GF2Vector vec, int col) {
        return vec.getBit(col) != 0;
    }

    /**
     * Returns the boolean entry from the given {@link GF2Matrix} at the given
     * position.
     *
     * @param mat The {@link GF2Matrix} to return the entry from.
     * @param row The row index to retrieve the entry from.
     * @param col The column index to retrieve the entry from.
     * @return The entry at the given position in the given {@link GF2Matrix}.
     */
    public static boolean getEntry(GF2Matrix mat, int row, int col) {
        int q = col >> 5;
        int r = col & 0x1f;
        return ((mat.getRow(row)[q] & (1 << r)) >>> r) != 0;
    }

    /**
     * Returns the submatrix between the given columns from the given
     * {@link GF2Matrix}.
     *
     * @param mat   The {@link GF2Matrix} to retrieve the submatrix from.
     * @param start The start index for the submatrix to extract.
     * @param end   The end index for the submatrix to extract.
     * @return The extracted submatrix.
     */
    public static GF2Matrix getColumns(GF2Matrix mat, int start, int end) {
        int numRows = mat.getNumRows();
        int amount = end - start;
        int length = (amount + 31) >> 5;
        int s = start & 0x1f;
        int[][] result = new int[numRows][length];

        for (int j = start; j < end; j++) {
            int p = (j - start) >> 5;
            int q = j >> 5;
            int r = j & 0x1f;
            for (int i = 0; i < numRows; i++) {
                int[] row = mat.getRow(i);
                int[] newRow = result[i];
                newRow[p] |= Integer.rotateRight(row[q] & (1 << r), s);
            }
        }

        return new GF2Matrix(amount, result);
    }

    /**
     * Sets the given bit in the given {@link GF2Matrix}.
     *
     * @param mat The {@link GF2Matrix} to modify.
     * @param row The row index to set the bit at.
     * @param col The column index to set the bit at.
     */
    public static void setBit(GF2Matrix mat, int row, int col) {
        mat.getRow(row)[col >> 5] |= 1 << (col & 0x1f);
    }

    /**
     * Initializes a zero {@link GF2Matrix} with the given dimensions.
     *
     * @param rows The amount of rows.
     * @param cols The amount of columns.
     * @return The generated zero {@link GF2Matrix}.
     */
    public static GF2Matrix zero(int rows, int cols) {
        try {
            Constructor<GF2Matrix> constructor = GF2Matrix.class
                    .getDeclaredConstructor(int.class, int.class);
            constructor.setAccessible(true);
            return constructor.newInstance(rows, cols);
        } catch (NoSuchMethodException
                | IllegalAccessException
                | InstantiationException
                | InvocationTargetException e) {
            throw new IllegalStateException("Can't access constructor");
        }
    }

    /**
     * Prints the given {@link GF2Vector} to the console.
     *
     * @param vec The {@link GF2Vector} to print.
     */
    public static void print(GF2Vector vec) {
        print(System.out, null, vec);
    }

    /**
     * Prints the given {@link GF2Vector} to the console as a MATLAB variable.
     *
     * @param varName The name of the variable for MATLAB import.
     * @param vec     The {@link GF2Vector} to print.
     */
    public static void print(String varName, GF2Vector vec) {
        print(System.out, varName, vec);
    }

    /**
     * Prints the given {@link GF2Vector} as a MATLAB variable to the given
     * {@link PrintStream}.
     *
     * @param printStream The {@link PrintStream} to print to.
     * @param varName     The name of the variable for MATLAB import.
     * @param vec         The {@link GF2Vector} to print.
     */
    public static void print(PrintStream printStream, String varName,
                             GF2Vector vec) {
        if (varName != null) {
            printStream.print(varName + " = [");
        }
        for (int i = 0; i < vec.getLength(); i++) {
            printStream.print(vec.getBit(i));
            printStream.print(' ');
        }
        if (varName != null) {
            printStream.print("]");
        }
        printStream.println(';');
    }

    /**
     * Prints the given {@link GF2Vector} to the console.
     *
     * @param mat The {@link GF2Vector} to print.
     */
    public static void print(GF2Matrix mat) {
        print(System.out, null, mat);
    }

    /**
     * Prints the given {@link GF2Matrix} to the console as a MATLAB variable.
     *
     * @param varName The name of the variable for MATLAB import.
     * @param mat     The {@link GF2Matrix} to print.
     */
    public static void print(String varName, GF2Matrix mat) {
        print(System.out, varName, mat);
    }

    /**
     * Prints the given {@link GF2Matrix} as a MATLAB variable to the given
     * {@link PrintStream}.
     *
     * @param printStream The {@link PrintStream} to print to.
     * @param varName     The name of the variable for MATLAB import.
     * @param mat         The {@link GF2Matrix} to print.
     */
    public static void print(PrintStream printStream, String varName,
                             GF2Matrix mat) {
        if (varName != null) {
            printStream.print(varName + " = [");
        }
        for (int i = 0; i < mat.getNumRows(); i++) {
            int[] row = mat.getRow(i);
            for (int j = 0; j < mat.getNumColumns(); j++) {
                int q = j >> 5;
                int r = j & 0x1f;
                printStream.print((row[q] & (1 << r)) >>> r);
                printStream.print(' ');
            }
            if (varName != null && i == mat.getNumRows() - 1) {
                printStream.print("]");
            }
            printStream.println(";");
        }
    }

    /**
     * Provides an {@link LDPC.IterationCallback} which prints the new vector
     * into MATLAB format after each iteration.
     */
    public static class MatlabPrinter implements LDPC.IterationCallback {

        private final String filePrefix;

        private int globalIter;

        public MatlabPrinter(String filePrefix) {
            this.filePrefix = filePrefix;
            this.globalIter = 0;
        }

        @Override
        public void onIteration(int iteration, GF2Vector data, double[] llr) {
            Path path = Paths.get(filePrefix + globalIter++ + ".m");
            try (OutputStream stream = Files.newOutputStream(path);
                 PrintStream printer = new PrintStream(stream)) {
                print(printer, "TEMP" + iteration, data);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

    }

}
