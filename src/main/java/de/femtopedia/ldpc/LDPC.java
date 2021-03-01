package de.femtopedia.ldpc;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;
import lombok.Getter;
import lombok.Setter;
import org.bouncycastle.pqc.math.linearalgebra.GF2Matrix;
import org.bouncycastle.pqc.math.linearalgebra.GF2Vector;
import org.bouncycastle.pqc.math.linearalgebra.Matrix;
import org.bouncycastle.pqc.math.linearalgebra.Vector;

import static de.femtopedia.ldpc.MatrixUtils.getColumns;
import static de.femtopedia.ldpc.MatrixUtils.getEntry;
import static de.femtopedia.ldpc.MatrixUtils.print;

/**
 * Provides LDPC code functionality using the {@link GF2Matrix} class and the
 * sum-product-algorithm (belief propagation).
 */
public class LDPC {

    /**
     * The generator matrix to use.
     */
    @Getter
    private final GF2Matrix g;

    /**
     * The parity check matrix to use.
     */
    @Getter
    private final GF2Matrix h;

    /**
     * The chance of a bitflip occurring.
     */
    @Getter
    @Setter
    private double bitflipChance;

    /**
     * The maximum iterations for the algorithm.
     */
    @Getter
    @Setter
    private int maxIterations;

    /**
     * The file prefix to print MATLAB matrices to.
     */
    private final String filePrefix;

    /**
     * The current iteration count for printing MATLAB matrices.
     */
    private int iteration = 0;

    /**
     * Initializes the LDPC instance using the given parity check matrix.
     *
     * @param h The parity check matrix to use.
     */
    public LDPC(GF2Matrix h) {
        this(h, 0.1, 20);
    }

    /**
     * Initializes the LDPC instance using the given parity check matrix and
     * bitflip chance.
     *
     * @param h             The parity check matrix to use.
     * @param bitflipChance The chance of a bitflip occurring.
     */
    public LDPC(GF2Matrix h, double bitflipChance) {
        this(h, bitflipChance, 20);
    }

    /**
     * Initializes the LDPC instance using the given parity check matrix,
     * bitflip chance and maximum amount of iterations for the algorithm.
     *
     * @param h             The parity check matrix to use.
     * @param bitflipChance The chance of a bitflip occurring.
     * @param maxIterations The maximum iterations for the algorithm.
     */
    public LDPC(GF2Matrix h, double bitflipChance, int maxIterations) {
        this(h, bitflipChance, maxIterations, null);
    }

    /**
     * Initializes the LDPC instance using the given parity check matrix,
     * bitflip chance and maximum amount of iterations for the algorithm.
     *
     * @param h             The parity check matrix to use.
     * @param bitflipChance The chance of a bitflip occurring.
     * @param maxIterations The maximum iterations for the algorithm.
     * @param filePrefix    The file prefix to print MATLAB matrices to.
     */
    public LDPC(GF2Matrix h, double bitflipChance, int maxIterations,
                String filePrefix) {
        this.filePrefix = filePrefix;
        this.bitflipChance = bitflipChance;
        this.maxIterations = maxIterations;
        this.h = h;
        g = getGeneratorMatrix(h);
    }

    /**
     * Returns the corresponding generator matrix for the given parity check
     * matrix.
     *
     * @param h The parity check matrix to calculate its generator matrix.
     * @return The corresponding generator matrix.
     */
    public static GF2Matrix getGeneratorMatrix(GF2Matrix h) {
        int m = h.getNumRows();
        int n = h.getNumColumns();
        int k = n - m;
        Matrix a = getColumns(h, 0, k).computeTranspose();
        Matrix b = getColumns(h, k, n).computeTranspose().computeInverse();
        return ((GF2Matrix) a.rightMultiply(b)).extendRightCompactForm();
    }

    /**
     * Encodes the given binary message using this LDPC instance.
     *
     * @param msg The message to encode.
     * @return The encoded message as boolean/binary matrix.
     */
    public GF2Vector encode(GF2Vector msg) {
        return (GF2Vector) g.leftMultiply(msg);
    }

    /**
     * Calculates the LLR (log-likelihood ratio) of the given binary value.
     *
     * @param value The binary value to calculate its LLR value.
     * @return The corresponding LLR value using this LDPC instance's
     *         {@link #bitflipChance}.
     */
    public double getLLR(boolean value) {
        return getLLR(value ? 1 : 0);
    }

    /**
     * Calculates the LLR (log-likelihood ratio) of the given value.
     *
     * @param value The value to calculate its LLR value.
     * @return The corresponding LLR value using this LDPC instance's
     *         {@link #bitflipChance}.
     */
    public double getLLR(double value) {
        return Math.log((1 - bitflipChance - value) / (bitflipChance - value));
    }

    /**
     * Tries to decode the given binary/boolean message using the sum-product-
     * algorithm/belief propagation.
     *
     * @param message The message to decode.
     * @return The possibly decoded message.
     */
    public GF2Vector decode(GF2Vector message) {
        int n = h.getNumColumns();

        double[] msg = new double[n];
        for (int i = 0; i < n; i++) {
            msg[i] = getLLR(getEntry(message, i));
        }

        return decode(msg);
    }

    /**
     * Tries to decode the given LLR value message using the sum-product-
     * algorithm/belief propagation.
     *
     * @param msg The message to decode.
     * @return The possibly decoded message.
     */
    @SuppressWarnings("unchecked")
    public GF2Vector decode(double[] msg) {
        Matrix hTrans = h.computeTranspose();
        int m = h.getNumRows();
        int n = h.getNumColumns();

        List<Integer>[] colAdj = new List[n];
        List<Integer>[] rowAdj = new List[m];
        initAdjacencyMatrices(colAdj, rowAdj);

        double[][] toCheckNodes = new double[m][n];
        double[][] fromCheckNodes = new double[m][n];
        initCheckNodes(rowAdj, toCheckNodes, msg);

        int iter = 0;
        GF2Vector estimate = decisionStep(msg);
        Vector syndrome = hTrans.leftMultiply(estimate);
        while (iter++ < maxIterations && !syndrome.isZero()) {
            // Step (i)
            updateSymbolNodes(rowAdj, toCheckNodes, fromCheckNodes);

            // Step (ii)
            updateCheckNodes(colAdj, toCheckNodes, fromCheckNodes, msg);
            double[] estimateLLR = getEstimate(colAdj, fromCheckNodes, msg);

            // Step (iii)
            estimate = decisionStep(estimateLLR);
            syndrome = hTrans.leftMultiply(estimate);
        }

        return estimate;
    }

    /**
     * Initializes the adjacency matrices from the parity check matrix
     * {@link #h}.
     *
     * @param colAdj The set saving for every column n its indices of 1's.
     * @param rowAdj The set saving for every row m its indices of 1's.
     */
    private void initAdjacencyMatrices(List<Integer>[] colAdj,
                                       List<Integer>[] rowAdj) {
        int n = h.getNumColumns();

        for (int i = 0; i < n; i++) {
            colAdj[i] = new ArrayList<>();
        }
        for (int i = 0; i < h.getNumRows(); i++) {
            rowAdj[i] = new ArrayList<>();
            for (int j = 0; j < n; j++) {
                if (getEntry(h, i, j)) {
                    colAdj[j].add(i);
                    rowAdj[i].add(j);
                }
            }
        }
    }

    /**
     * Initializes the check node values.
     *
     * @param rowAdj  The set saving for every row m its indices of 1's.
     * @param ingoing The ingoing messages to the check nodes.
     * @param msg     The originally received message.
     */
    private void initCheckNodes(List<Integer>[] rowAdj, double[][] ingoing,
                                double[] msg) {
        for (int i = 0; i < h.getNumRows(); i++) {
            for (int j : rowAdj[i]) {
                ingoing[i][j] = msg[j];
            }
        }
    }

    /**
     * Updates the symbol nodes' values.
     *
     * @param rowAdj   The set saving for every row m its indices of 1's.
     * @param ingoing  The ingoing messages to the check nodes.
     * @param outgoing The outgoing messages from the check nodes.
     */
    private void updateSymbolNodes(List<Integer>[] rowAdj, double[][] ingoing,
                                   double[][] outgoing) {
        for (int i = 0; i < h.getNumRows(); i++) {
            for (int j : rowAdj[i]) {
                double prod = 1;
                for (int k : rowAdj[i]) {
                    if (k != j) {
                        prod *= Math.tanh(ingoing[i][k] / 2);
                    }
                }
                outgoing[i][j] = 2 * atanh(prod);
            }
        }
    }

    /**
     * Updates the check nodes' values.
     *
     * @param colAdj   The set saving for every column n its indices of 1's.
     * @param ingoing  The ingoing messages to the check nodes.
     * @param outgoing The outgoing messages from the check nodes.
     * @param msg      The originally received message.
     */
    private void updateCheckNodes(List<Integer>[] colAdj, double[][] ingoing,
                                  double[][] outgoing, double[] msg) {
        for (int j = 0; j < h.getNumColumns(); j++) {
            for (int i : colAdj[j]) {
                double sum = 0;
                for (int k : colAdj[j]) {
                    if (k != i) {
                        sum += outgoing[k][j];
                    }
                }
                ingoing[i][j] = msg[j] + sum;
            }
        }
    }

    /**
     * Returns an estimate of the decoded matrix based on the current nodes'
     * values.
     *
     * @param colAdj   The set saving for every column n its indices of 1's.
     * @param outgoing The outgoing messages from the check nodes.
     * @param msg      The originally received message.
     * @return The current LLR value estimate.
     */
    private double[] getEstimate(List<Integer>[] colAdj, double[][] outgoing,
                                 double[] msg) {
        int n = h.getNumColumns();

        double[] l = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int k : colAdj[i]) {
                sum += outgoing[k][i];
            }
            l[i] = msg[i] + sum;
        }
        return l;
    }

    /**
     * Decides for each LLR value whether it is more likely that it is a 0 or 1.
     * This return value represents the current estimate for the decoded message
     * in binary form.
     *
     * @param llrValues The LLR values to decide on.
     * @return The current estimate for the decoded message in binary form.
     */
    private GF2Vector decisionStep(double[] llrValues) {
        GF2Vector vec = new GF2Vector(llrValues.length);
        IntStream.range(0, llrValues.length)
                .filter(i -> llrValues[i] < 0)
                .forEach(vec::setBit);

        // Print as MATLAB matrix if wanted.
        if (filePrefix != null) {
            Path path = Paths.get(filePrefix + iteration++ + ".m");
            try (OutputStream stream = Files.newOutputStream(path);
                 PrintStream printer = new PrintStream(stream)) {
                print(printer, "TEMP", vec);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        return vec;
    }

    /**
     * Calculates the arctanh(x) of a given value x.
     *
     * @param x The value to calculate the arctanh(x) of.
     * @return The arctanh(x) of the given value x.
     */
    private static double atanh(double x) {
        return Math.log((1 + x) / (1 - x)) / 2;
    }

}
