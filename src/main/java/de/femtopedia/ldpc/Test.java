package de.femtopedia.ldpc;

public class Test {

    public static void main(String[] args) {
        int[][] dataParityCheck = new int[][]{
                {0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1},
                {1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0},
                {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0},
                {0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0},
                {1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0},
                {1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0}};
        int[][] dataMsg = new int[][]{{1, 1, 1, 0, 0, 1}};

        BinaryMatrix H = new BinaryMatrix(dataParityCheck);
        LDPC ldpc = new LDPC(H, 0.1, 20);
        BinaryMatrix MSG = new BinaryMatrix(dataMsg);
        BinaryMatrix ENC = ldpc.encode(MSG);

        boolean[][] recvData = ENC.getData();
        recvData[0][6] = !recvData[0][6];
        BinaryMatrix RMSG = new BinaryMatrix(recvData);

        System.out.println("Generator Matrix:");
        System.out.println(ldpc.getGeneratorMatrix());
        System.out.println();
        System.out.println("Parity Check Matrix:");
        System.out.println(ldpc.getParityCheckMatrix());
        System.out.println();
        System.out.println("Test message:");
        System.out.println(MSG);
        System.out.println();
        System.out.println("Test message encoded:");
        System.out.println(ENC);
        System.out.println();
        System.out.println("Received message decoded:");
        System.out.println(ldpc.decode(RMSG));
        System.out.println();
    }

}
