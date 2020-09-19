package de.femtopedia.ldpc;

import java.util.Arrays;

public class ArrayUtils {

    public static void insertionSortParallel(boolean[][] a, boolean[][] b) {
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

    public static int compare(boolean[] a, boolean[] b) {
        for (int i = 0; i < a.length; i++) {
            int comp = Boolean.compare(a[i], b[i]);
            if (comp != 0) {
                return comp;
            }
        }
        return 0;
    }

    public static void swap(boolean[][] a, int i1, int i2) {
        boolean[] temp = Arrays.copyOf(a[i1], a[i1].length);
        a[i1] = Arrays.copyOf(a[i2], a[i2].length);
        a[i2] = temp;
    }

    public static void xor(boolean[] a, boolean[] b) {
        for (int i = 0; i < a.length; i++) {
            a[i] ^= b[i];
        }
    }

}
