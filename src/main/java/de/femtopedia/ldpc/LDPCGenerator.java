package de.femtopedia.ldpc;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import org.bouncycastle.pqc.math.linearalgebra.GF2Matrix;

import static de.femtopedia.ldpc.util.MatrixUtils.setBit;
import static de.femtopedia.ldpc.util.MatrixUtils.zero;

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

    /**
     * Reads the {@link GF2Matrix} from the given {@link Path}.
     *
     * @param path The {@link Path} to read from.
     * @return The parsed {@link GF2Matrix}.
     */
    public static GF2Matrix readBinaryMatrix(Path path) {
        GF2Matrix matrix = null;
        try {
            matrix = new GF2Matrix(Files.readAllBytes(path));
        } catch (IOException e) {
            e.printStackTrace();
        }
        return matrix;
    }

    /**
     * Writes the {@link GF2Matrix} to the given {@link Path}.
     *
     * @param path The {@link Path} to write to.
     * @param mat  The {@link GF2Matrix} to be written.
     */
    public static void writeBinaryMatrix(Path path, GF2Matrix mat) {
        try {
            Files.write(path, mat.getEncoded());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
