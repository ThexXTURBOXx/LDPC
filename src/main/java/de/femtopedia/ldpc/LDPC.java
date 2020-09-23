package de.femtopedia.ldpc;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import lombok.Getter;
import lombok.Setter;

/**
 * Provides LDPC code functionality using the {@link BinaryMatrix} class and the
 * sum-product-algorithm (belief propagation).
 */
public class LDPC {

    /**
     * The generator matrix to use.
     */
    @Getter
    private final BinaryMatrix g;

    /**
     * The parity check matrix to use.
     */
    @Getter
    private final BinaryMatrix h;

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
     * Initializes the LDPC instance using the given parity check matrix.
     *
     * @param h The parity check matrix to use.
     */
    public LDPC(BinaryMatrix h) {
        this(h, 0.1, 20);
    }

    /**
     * Initializes the LDPC instance using the given parity check matrix and
     * bitflip chance.
     *
     * @param h             The parity check matrix to use.
     * @param bitflipChance The chance of a bitflip occurring.
     */
    public LDPC(BinaryMatrix h, double bitflipChance) {
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
    public LDPC(BinaryMatrix h, double bitflipChance, int maxIterations) {
        this.bitflipChance = bitflipChance;
        this.maxIterations = maxIterations;
        this.h = h;
        g = getGeneratorMatrix(h);
    }

    /**
     * Generates a Gallager Matrix from the given parameters using the Gallager
     * Ensemble Algorithm.
     *
     * @param k  Length of the messages to decode.
     * @param l  Amount of redundancy bits to use.
     * @param dv Amount of 1's per column, i.e. the amount of parity check sets.
     * @return A Gallager Matrix from the given parameters.
     */
    public static BinaryMatrix generateGallagerMatrix(int k, int l, int dv) {
        int n = k + l;
        int m = n - k;
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
     * Returns the corresponding generator matrix for the given parity check
     * matrix.
     *
     * @param h The parity check matrix to calculate its generator matrix.
     * @return The corresponding generator matrix.
     */
    public static BinaryMatrix getGeneratorMatrix(BinaryMatrix h) {
        int k = h.getCols() - h.getRows();
        BinaryMatrix a = h.getColumns(0, k);
        BinaryMatrix b = h.getColumns(k, h.getCols());
        if (b.isInvertible()) {
            return BinaryMatrix.horizConcat(BinaryMatrix.eye(k),
                    a.transpose().mult(b.transpose().inv()));
        } else {
            throw new IllegalArgumentException("Matrix B is not invertible!");
        }
    }

    /**
     * Encodes the given binary message using this LDPC instance.
     *
     * @param msg The message to encode.
     * @return The encoded message as boolean/binary matrix.
     */
    public BinaryMatrix encode(BinaryMatrix msg) {
        return msg.mult(g);
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
    public BinaryMatrix decode(BinaryMatrix message) {
        int n = h.getCols();

        double[] msg = new double[n];
        for (int i = 0; i < n; i++) {
            msg[i] = getLLR(message.getEntry(0, i));
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
    public BinaryMatrix decode(double[] msg) {
        int m = h.getRows();
        int n = h.getCols();

        List<Integer>[] colAdj = new List[n];
        List<Integer>[] rowAdj = new List[m];
        initAdjacencyMatrices(colAdj, rowAdj);

        double[][] toCheckNodes = new double[m][n];
        double[][] fromCheckNodes = new double[m][n];
        initCheckNodes(rowAdj, toCheckNodes, msg);

        int iter = 0;
        BinaryMatrix estimate = decisionStep(msg);
        BinaryMatrix syndrome = estimate.mult(h.transpose());
        while (iter++ < maxIterations && syndrome.sum() != 0) {
            // Step (i)
            updateSymbolNodes(rowAdj, toCheckNodes, fromCheckNodes);

            // Step (ii)
            updateCheckNodes(colAdj, toCheckNodes, fromCheckNodes, msg);
            double[] estimateLLR = getEstimate(colAdj, fromCheckNodes, msg);

            // Step (iii)
            estimate = decisionStep(estimateLLR);
            syndrome = estimate.mult(h.transpose());
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
        int n = h.getCols();

        for (int i = 0; i < n; i++) {
            colAdj[i] = new ArrayList<>();
        }
        for (int i = 0; i < h.getRows(); i++) {
            rowAdj[i] = new ArrayList<>();
            for (int j = 0; j < n; j++) {
                if (h.getEntry(i, j)) {
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
        for (int i = 0; i < h.getRows(); i++) {
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
        for (int i = 0; i < h.getRows(); i++) {
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
        for (int j = 0; j < h.getCols(); j++) {
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
        int n = h.getCols();

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
    private BinaryMatrix decisionStep(double[] llrValues) {
        int n = llrValues.length;
        boolean[][] estimateData = new boolean[1][n];
        for (int i = 0; i < n; i++) {
            estimateData[0][i] = llrValues[i] < 0;
        }
        return new BinaryMatrix(estimateData);
    }

    /**
     * Calculates the arctanh(x) of a given value x.
     *
     * @param x The value to calculate the arctanh(x) of.
     * @return The arctanh(x) of the given value x.
     */
    private double atanh(double x) {
        return Math.log((1 + x) / (1 - x)) / 2;
    }

}
