package de.femtopedia.ldpc;

/**
 * Provides a simple example as test case.
 */
public final class Test {

    /**
     * Don't initialize me.
     */
    private Test() {
    }

    /**
     * Starting point of the program.
     *
     * @param args The starting arguments.
     */
    public static void main(String[] args) {
        int[][] dataParityCheck = new int[][]{
                {0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1},
                {1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0},
                {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0},
                {0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0},
                {1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0},
                {1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0}};
        int[][] dataMsg = new int[][]{{1, 1, 1, 0, 0, 1}};

        BinaryMatrix h = new BinaryMatrix(dataParityCheck);
        LDPC ldpc = new LDPC(h, 0.1, 20);
        BinaryMatrix msg = new BinaryMatrix(dataMsg);
        BinaryMatrix encoded = ldpc.encode(msg);

        boolean[][] recvData = encoded.getData();
        recvData[0][6] = !recvData[0][6];
        BinaryMatrix recvMsg = new BinaryMatrix(recvData);

        System.out.println("Generator Matrix:");
        System.out.println(ldpc.getG());
        System.out.println();
        System.out.println("Parity Check Matrix:");
        System.out.println(ldpc.getH());
        System.out.println();
        System.out.println("Test message:");
        System.out.println(msg);
        System.out.println();
        System.out.println("Test message encoded:");
        System.out.println(encoded);
        System.out.println();
        System.out.println("Received message decoded:");
        System.out.println(ldpc.decode(recvMsg));
        System.out.println();
    }

}
