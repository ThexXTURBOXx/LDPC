package de.femtopedia.ldpc;

import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Provides some methods for generating Parity Check Matrices.
 */
public final class LDPCGenerator {

    /**
     * Don't initialize me.
     */
    private LDPCGenerator() {
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

        return BinaryMatrix.fromMatrixFunction(m, n, p, p,
                (i, j) -> h2.getEntry(i, j)
                        ? iArr[Math.abs(new Random(i * n + j).nextInt() % p)]
                        : BinaryMatrix.zero(p, p));
    }

}
