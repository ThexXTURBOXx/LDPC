package de.femtopedia.ldpc;

import java.util.ArrayList;
import java.util.List;

public class LDPC {

    private final BinaryMatrix g;
    private final BinaryMatrix h;
    private final double bitflipChance;
    private final int maxIterations;

    public LDPC(BinaryMatrix h, double bitflipChance, int maxIterations) {
        this.bitflipChance = bitflipChance;
        this.maxIterations = maxIterations;
        this.h = h;
        g = getGeneratorMatrix(h);
    }

    public static BinaryMatrix getGeneratorMatrix(BinaryMatrix h) {
        int k = h.getCols() - h.getRows();
        BinaryMatrix A = h.getColumns(0, k);
        BinaryMatrix B = h.getColumns(k, h.getCols());
        return BinaryMatrix.concat(BinaryMatrix.eye(k),
                A.transpose().mult(B.transpose().inv()));
    }

    public BinaryMatrix getGeneratorMatrix() {
        return g;
    }

    public BinaryMatrix getParityCheckMatrix() {
        return h;
    }

    public BinaryMatrix encode(BinaryMatrix msg) {
        return msg.mult(g);
    }

    @SuppressWarnings("unchecked")
    public BinaryMatrix decode(BinaryMatrix msg) {
        int m = h.getRows();
        int n = h.getCols();
        double[] L = new double[n];
        for (int i = 0; i < n; i++) {
            L[i] = Math.log(msg.getEntry(0, i)
                    ? (bitflipChance / (1 - bitflipChance))
                    : ((1 - bitflipChance) / bitflipChance));
        }
        List<Integer>[] M = new List[n];
        List<Integer>[] N = new List[m];
        for (int i = 0; i < n; i++) {
            M[i] = new ArrayList<>();

        }
        for (int i = 0; i < m; i++) {
            N[i] = new ArrayList<>();
            for (int j = 0; j < n; j++) {
                if (h.getEntry(i, j)) {
                    M[j].add(i);
                    N[i].add(j);
                }
            }
        }
        double[][] lambda = new double[m][n];
        double[][] Lambda = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j : N[i]) {
                lambda[i][j] = L[j];
            }
        }

        int iter = 0;
        BinaryMatrix estimate = msg;
        BinaryMatrix syndrome = msg.mult(h.transpose());
        while (iter++ < maxIterations
                && syndrome.reduce(0, Integer::sum) != 0) {
            // Step (i)
            for (int i = 0; i < m; i++) {
                for (int j : N[i]) {
                    double prod = 1;
                    for (int k : N[i]) {
                        if (k != j) {
                            prod *= Math.tanh(lambda[i][k] / 2);
                        }
                    }
                    Lambda[i][j] = 2 * atanh(prod);
                }
            }

            // Step (ii)
            for (int j = 0; j < n; j++) {
                for (int i : M[j]) {
                    double sum = 0;
                    for (int k : M[j]) {
                        if (k != i) {
                            sum += Lambda[k][j];
                        }
                    }
                    lambda[i][j] = L[j] + sum;
                }
            }

            double[] l = new double[n];
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for (int k : M[i]) {
                    sum += Lambda[k][i];
                }
                l[i] = L[i] + sum;
            }

            // Step (iii)
            boolean[][] estimateData = new boolean[1][n];
            for (int i = 0; i < n; i++) {
                estimateData[0][i] = l[i] < 0;
            }
            estimate = new BinaryMatrix(estimateData);
            syndrome = estimate.mult(h.transpose());
        }

        return estimate;
    }

    private double atanh(double x) {
        return Math.log((1 + x) / (1 - x)) / 2;
    }

}
