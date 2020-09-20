package de.femtopedia.ldpc;

import java.util.Arrays;

/**
 * Helper class providing access to some operations on boolean/binary arrays.
 */
public final class ArrayUtils {

    /**
     * Don't initialize me.
     */
    private ArrayUtils() {
    }

    /**
     * Performs a sorting algorithm descending on one 2d array's rows.
     *
     * @param a The array to sort.
     */
    public static void sortDesc(boolean[][] a) {
        Arrays.sort(a, (b1, b2) -> -ArrayUtils.compare(b1, b2));
    }

    /**
     * Performs insertion sort descending on one 2d array's rows and performs
     * the same swaps on another array.
     *
     * @param a The array to sort.
     * @param b The array to also perform the steps on.
     */
    public static void parallelSortDesc(boolean[][] a, boolean[][] b) {
        for (int i = 0; i < a.length; i++) {
            for (int j = i; j > 0; j--) {
                if (ArrayUtils.compare(a[j - 1], a[j]) < 0) {
                    ArrayUtils.swap(a, j, j - 1);
                    ArrayUtils.swap(b, j, j - 1);
                } else {
                    break;
                }
            }
        }
    }

    /**
     * Compares two boolean/binary arrays.
     *
     * @param a The array to compare to.
     * @param b The array to compare with.
     * @return the value {@code 0} if {@code a == b};
     *         a value less than {@code 0} if {@code a < b}; and
     *         a value greater than {@code 0} if {@code a > b}.
     */
    public static int compare(boolean[] a, boolean[] b) {
        for (int i = 0; i < a.length; i++) {
            int comp = Boolean.compare(a[i], b[i]);
            if (comp != 0) {
                return comp;
            }
        }
        return 0;
    }

    /**
     * Swaps the given indices in a given 2d boolean/binary array.
     *
     * @param a  The array to swap entries in.
     * @param i1 The first index to swap.
     * @param i2 The second index to swap.
     */
    public static void swap(boolean[][] a, int i1, int i2) {
        boolean[] temp = Arrays.copyOf(a[i1], a[i1].length);
        a[i1] = Arrays.copyOf(a[i2], a[i2].length);
        a[i2] = temp;
    }

    /**
     * Performs the XOR operation (binary addition) on two boolean/binary
     * arrays.
     *
     * @param a The first array (gets changed).
     * @param b The second array (doesn't get changed).
     */
    public static void xor(boolean[] a, boolean[] b) {
        for (int i = 0; i < a.length; i++) {
            a[i] ^= b[i];
        }
    }

}
