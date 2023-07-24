package de.femtopedia.ldpc;

import de.femtopedia.ldpc.util.MatrixUtils;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;
import lombok.Getter;
import lombok.Setter;
import org.bouncycastle.pqc.legacy.math.linearalgebra.GF2Matrix;
import org.bouncycastle.pqc.legacy.math.linearalgebra.GF2Vector;
import org.bouncycastle.pqc.legacy.math.linearalgebra.Matrix;
import org.bouncycastle.pqc.legacy.math.linearalgebra.Vector;

import static de.femtopedia.ldpc.util.MatrixUtils.getColumns;
import static de.femtopedia.ldpc.util.MatrixUtils.getEntry;

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
     * The consumer being called after each iteration of the decoder.
     */
    @Getter
    private final IterationCallback consumer;

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
     * The column-indexed adjacency matrix.
     */
    private final List<Integer>[] colAdj;

    /**
     * The row-indexed adjacency matrix.
     */
    private final List<Integer>[] rowAdj;

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
     * bitflip chance and maximum amount of iterations for the algorithm and
     * calculates the generator matrix.
     *
     * @param h             The parity check matrix to use.
     * @param bitflipChance The chance of a bitflip occurring.
     * @param maxIterations The maximum iterations for the algorithm.
     * @param consumer      The consumer being called after each iteration of
     *                      the decoder.
     */
    public LDPC(GF2Matrix h, double bitflipChance, int maxIterations,
                IterationCallback consumer) {
        this(getGeneratorMatrix(h), h,
                bitflipChance, maxIterations, consumer);
    }

    /**
     * Initializes the LDPC instance using the given parity check matrix,
     * bitflip chance and maximum amount of iterations for the algorithm using
     * the given pre-computed generator matrix without checking its validity.
     *
     * @param g             The pre-computed generator matrix to use.
     * @param h             The parity check matrix to use.
     * @param bitflipChance The chance of a bitflip occurring.
     * @param maxIterations The maximum iterations for the algorithm.
     */
    public LDPC(GF2Matrix g, GF2Matrix h, double bitflipChance,
                int maxIterations) {
        this(g, h, bitflipChance, maxIterations, null);
    }

    /**
     * Initializes the LDPC instance using the given parity check matrix,
     * bitflip chance and maximum amount of iterations for the algorithm using
     * the given pre-computed generator matrix without checking its validity.
     *
     * @param g             The pre-computed generator matrix to use.
     * @param h             The parity check matrix to use.
     * @param bitflipChance The chance of a bitflip occurring.
     * @param maxIterations The maximum iterations for the algorithm.
     * @param consumer      The consumer being called after each iteration of
     *                      the decoder.
     */
    @SuppressWarnings("unchecked")
    public LDPC(GF2Matrix g, GF2Matrix h, double bitflipChance,
                int maxIterations, IterationCallback consumer) {
        this.g = g;
        this.h = h;
        this.bitflipChance = bitflipChance;
        this.maxIterations = maxIterations;
        this.consumer = consumer;

        colAdj = new List[h.getNumColumns()];
        rowAdj = new List[h.getNumRows()];
        initAdjacencyMatrices(colAdj, rowAdj);
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
        for (int i = 0; i < getParityBits(); i++) {
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
     * Returns the number of bits that an encoded message consists of (n).
     *
     * @return The number of bits in an encoded message (n).
     */
    public int getEncodedBits() {
        return g.getNumColumns();
    }

    /**
     * Returns the number of message bits in an encoded message (k).
     *
     * @return The number of message bits in an encoded message (k).
     */
    public int getMessageBits() {
        return g.getNumRows();
    }

    /**
     * Returns the number of parity bits in an encoded message (n-k).
     *
     * @return The number of parity bits in an encoded message (n-k).
     */
    public int getParityBits() {
        return h.getNumRows();
    }

    /**
     * Extracts the data from the given {@link GF2Vector} (may be padded).
     *
     * @param vec The code vector.
     * @return The extracted data vector.
     */
    public GF2Vector extractData(GF2Vector vec) {
        int dataBits = getMessageBits();
        int bits = dataBits + getParityBits();
        int amount = vec.getLength() / bits;
        GF2Vector data = new GF2Vector(amount * dataBits);

        IntStream.range(0, amount)
                .forEach(i -> IntStream.range(0, dataBits)
                        .filter(j -> vec.getBit(i * bits + j) != 0)
                        .forEach(j -> data.setBit(i * dataBits + j)));

        return data;
    }

    /**
     * Extracts the data from the given {@link GF2Vector} cut off at the given
     * index.
     *
     * @param vec     The code vector.
     * @param dataLen The end index to cut off at.
     * @return The extracted data vector.
     */
    public GF2Vector extractData(GF2Vector vec, int dataLen) {
        int dataBits = getMessageBits();
        int bits = dataBits + getParityBits();
        int amount = vec.getLength() / bits;
        GF2Vector data = new GF2Vector(dataLen);

        IntStream.range(0, amount)
                .forEach(i -> IntStream.range(0, dataBits)
                        .filter(j -> i * dataBits + j < dataLen)
                        .filter(j -> vec.getBit(i * bits + j) != 0)
                        .forEach(j -> data.setBit(i * dataBits + j)));

        return data;
    }

    /**
     * Preprocesses the given {@link GF2Vector} by splitting it into some parts
     * that fit the current parity check matrix.
     *
     * @param len The length to pad to.
     * @param vec The {@link GF2Vector} to preprocess.
     * @return Some {@link GF2Vector}s fitting the current scheme.
     */
    private GF2Vector[] preprocess(int len, GF2Vector vec) {
        if (len == vec.getLength()) {
            return new GF2Vector[]{vec};
        }

        int parts = (int) Math.ceil((double) vec.getLength() / len);
        GF2Vector[] vectors = new GF2Vector[parts];
        for (int i = 0; i < parts; i++) {
            int offset = i * len;
            GF2Vector part = new GF2Vector(len);
            for (int j = 0; j < len; j++) {
                int index = offset + j;
                if (index >= vec.getLength()) {
                    break;
                }
                if (MatrixUtils.getEntry(vec, index)) {
                    part.setBit(j);
                }
            }
            vectors[i] = part;
        }
        return vectors;
    }

    /**
     * Postprocesses the given {@link GF2Vector}s by extracting their data and
     * concatenating them to a data {@link GF2Vector}.
     *
     * @param vectors The {@link GF2Vector}s to process.
     * @return The extracted and concatenated data.
     */
    private GF2Vector postprocess(GF2Vector[] vectors) {
        int length = Arrays.stream(vectors).mapToInt(Vector::getLength).sum();
        if (vectors.length == 1 && length == h.getNumColumns()) {
            return vectors[0];
        }

        GF2Vector result = new GF2Vector(length);
        for (int i = 0; i < vectors.length; i++) {
            GF2Vector vec = vectors[i];
            int len = vec.getLength();
            int offset = i * len;
            IntStream.range(0, len)
                    .filter(j -> getEntry(vec, j))
                    .forEach(j -> result.setBit(offset + j));
        }
        return result;
    }

    /**
     * Encodes the given binary message using this LDPC instance (pads if
     * necessary).
     *
     * @param msg The message to encode.
     * @return The encoded message as boolean/binary vector.
     */
    public GF2Vector encode(GF2Vector msg) {
        return postprocess(Arrays.stream(preprocess(getMessageBits(), msg))
                .map(v -> (GF2Vector) g.leftMultiply(v))
                .toArray(GF2Vector[]::new));
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
        return Math.log(Math.abs((1 - bitflipChance - value)
                / (bitflipChance - value)));
    }

    /**
     * Tries to decode the given binary/boolean message using the sum-product-
     * algorithm/belief propagation.
     *
     * @param message The message to decode.
     * @return The possibly decoded message.
     */
    public GF2Vector decode(GF2Vector message) {
        int n = getEncodedBits();
        return postprocess(Arrays.stream(preprocess(n, message))
                .map(v -> IntStream.range(0, n)
                        .mapToDouble(i -> getLLR(getEntry(v, i)))
                        .toArray())
                .map(this::decode)
                .toArray(GF2Vector[]::new));
    }

    /**
     * Tries to decode the given LLR value message using the sum-product-
     * algorithm/belief propagation.
     *
     * @param msg The message to decode.
     * @return The possibly decoded message.
     */
    public GF2Vector decode(double[] msg) {
        int m = getParityBits();
        int n = h.getNumColumns();

        double[][] toCheckNodes = new double[m][n];
        double[][] fromCheckNodes = new double[m][n];
        initCheckNodes(toCheckNodes, msg);

        int iter = 0;
        GF2Vector estimate = decisionStep(iter, msg);
        Vector syndrome = h.rightMultiply(estimate);
        while (iter++ < maxIterations && !syndrome.isZero()) {
            // Step (i)
            updateSymbolNodes(toCheckNodes, fromCheckNodes);

            // Step (ii)
            updateCheckNodes(toCheckNodes, fromCheckNodes, msg);
            double[] estimateLLR = getEstimate(fromCheckNodes, msg);

            // Step (iii)
            estimate = decisionStep(iter, estimateLLR);
            syndrome = h.rightMultiply(estimate);
        }

        return estimate;
    }

    /**
     * Initializes the check node values.
     *
     * @param ingoing The ingoing messages to the check nodes.
     * @param msg     The originally received message.
     */
    private void initCheckNodes(double[][] ingoing, double[] msg) {
        for (int i = 0; i < getParityBits(); i++) {
            for (int j : rowAdj[i]) {
                ingoing[i][j] = msg[j];
            }
        }
    }

    /**
     * Updates the symbol nodes' values.
     *
     * @param ingoing  The ingoing messages to the check nodes.
     * @param outgoing The outgoing messages from the check nodes.
     */
    private void updateSymbolNodes(double[][] ingoing, double[][] outgoing) {
        for (int i = 0; i < getParityBits(); i++) {
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
     * @param ingoing  The ingoing messages to the check nodes.
     * @param outgoing The outgoing messages from the check nodes.
     * @param msg      The originally received message.
     */
    private void updateCheckNodes(double[][] ingoing, double[][] outgoing,
                                  double[] msg) {
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
     * @param outgoing The outgoing messages from the check nodes.
     * @param msg      The originally received message.
     * @return The current LLR value estimate.
     */
    private double[] getEstimate(double[][] outgoing, double[] msg) {
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
     * @param iter      The current iteration count.
     * @param llrValues The LLR values to decide on.
     * @return The current estimate for the decoded message in binary form.
     */
    private GF2Vector decisionStep(int iter, double[] llrValues) {
        GF2Vector vec = new GF2Vector(llrValues.length);
        IntStream.range(0, llrValues.length)
                .filter(i -> llrValues[i] < 0)
                .forEach(vec::setBit);

        // Call consumer if available
        if (consumer != null) {
            consumer.onIteration(iter, vec, llrValues);
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

    /**
     * Provides a consumer which is being called after each iteration of the
     * decoding algorithm.
     */
    public interface IterationCallback {

        /**
         * Consumes the given parameters and executes an action.
         *
         * @param iteration The current iteration count.
         * @param data      The new data after the decision step.
         * @param llr       The current LLR Values before the decision step.
         */
        void onIteration(int iteration, GF2Vector data, double[] llr);

    }

}
