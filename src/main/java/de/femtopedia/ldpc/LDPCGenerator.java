package de.femtopedia.ldpc;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.SplittableRandom;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.bouncycastle.pqc.math.linearalgebra.GF2Matrix;

import static de.femtopedia.ldpc.MatrixUtils.setBit;
import static de.femtopedia.ldpc.MatrixUtils.zero;

/**
 * Provides some methods for generating Parity Check Matrices.
 */
public final class LDPCGenerator {

    /**
     * Don't initialize me.
     */
    private LDPCGenerator() {
        throw new UnsupportedOperationException();
    }

    /**
     * Generates a Gallager Matrix from the given parameters using the Gallager
     * Ensemble Algorithm.
     *
     * @param k  Length of the messages to decode.
     * @param m  Amount of redundancy bits to use.
     * @param dv Amount of 1's per column, i.e. the amount of parity check sets.
     * @return A Gallager Matrix from the given parameters.
     */
    @Deprecated
    public static BinaryMatrix generateGallagerMatrix(int k, int m, int dv) {
        int n = k + m;
        int dc = dv * n / m;

        if ((m % dv) != 0 || ((dv * n) % m) != 0) {
            throw new IllegalArgumentException("Invalid parameters for the "
                    + "Gallager ensemble algorithm!");
        }

        BinaryMatrix[] matrices = new BinaryMatrix[dv];
        matrices[0] = BinaryMatrix.fromFunction(m / dv, n,
                (i, j) -> i * dc <= j && j < (i + 1) * dc);

        for (int i = 1; i < dv; i++) {
            List<Integer> permutation = IntStream.range(0, n).boxed()
                    .collect(Collectors.toList());
            Collections.shuffle(permutation);

            matrices[i] = matrices[0].permuteColumns(permutation);
        }

        return BinaryMatrix.vertConcat(matrices);
    }

    /**
     * Generates a matrix with girth {@literal >}= 8 (free from 4 and 6-cycles).
     *
     * @param v The general size of the matrix.
     * @param p The expansion factor.
     * @return A quasi-cyclic matrix.
     */
    @Deprecated
    public static BinaryMatrix generateCycleFreeMatrix(int v, int p) {
        BinaryMatrix[] dArr = new BinaryMatrix[v * v];
        dArr[0] = BinaryMatrix.fromFunction(v, v * v, (i, j) -> j == 0);
        for (int i = 1; i < v * v; i++) {
            dArr[i] = dArr[0].shiftRight(i);
        }
        BinaryMatrix d = BinaryMatrix.vertConcat(dArr);

        BinaryMatrix[] eArr = new BinaryMatrix[v];
        BinaryMatrix[] e0Arr = IntStream.range(0, v)
                .mapToObj(i -> BinaryMatrix.fromFunction(v, v * v,
                        Integer::equals))
                .toArray(BinaryMatrix[]::new);
        eArr[0] = BinaryMatrix.vertConcat(e0Arr);
        for (int i = 1; i < v; i++) {
            eArr[i] = eArr[0].shiftRight(i * v);
        }
        BinaryMatrix e = BinaryMatrix.vertConcat(eArr);

        BinaryMatrix[] fArr = new BinaryMatrix[v];
        fArr[0] = BinaryMatrix.fromFunction(v, v * v, (i, j) -> i * v == j);
        for (int i = 1; i < v; i++) {
            fArr[i] = fArr[0].shiftRight(i);
        }
        BinaryMatrix[] fv = IntStream.range(0, v)
                .mapToObj(i -> BinaryMatrix.vertConcat(fArr))
                .toArray(BinaryMatrix[]::new);
        BinaryMatrix f = BinaryMatrix.vertConcat(fv);

        BinaryMatrix h2 = BinaryMatrix.horizConcat(d, e, f).transpose();

        BinaryMatrix[] iArr = IntStream.range(0, p)
                .mapToObj(i -> BinaryMatrix.eye(p).shiftRight(i))
                .toArray(BinaryMatrix[]::new);
        int m = h2.getRows();
        int n = h2.getCols();

        BinaryMatrix zero = BinaryMatrix.zero(p, p);

        return BinaryMatrix.fromMatrixFunction(m, n, p, p,
                (i, j) -> h2.getEntry(i, j)
                        ? iArr[Math.abs(new SplittableRandom((long) i * n + j)
                        .nextInt() % p)]
                        : zero);
    }

    /**
     * Reads a sparse Matrix from the given alist-file.
     * Format: http://www.inference.org.uk/mackay/codes/alist.html
     *
     * @param path The path to the alist file.
     * @return The parsed {@link GF2Matrix}.
     */
    public static GF2Matrix readAList(Path path) {
        try (BufferedReader reader = Files.newBufferedReader(path)) {
            int[] dims = Arrays.stream(reader.readLine().split("\\s+"))
                    .mapToInt(Integer::parseInt).toArray();
            // Ignore weights
            reader.readLine();
            // Ignore row weights
            reader.readLine();
            // Ignore column weights
            reader.readLine();

            GF2Matrix mat = zero(dims[0], dims[1]);
            for (int row = 0; row < dims[0]; row++) {
                String line = reader.readLine();
                final int i = row;
                Arrays.stream(line.split("\\s+"))
                        .mapToInt(Integer::parseInt)
                        .filter(j -> j > 0)
                        .forEach(j -> setBit(mat, i, j - 1));
            }
            return mat;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

}
