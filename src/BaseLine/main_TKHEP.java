package BaseLine;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import Stock_Package.Stock;
import Tool_Package.MemoryLogger;

public class main_TKHEP {

	public static void main(String[] args) throws IOException {

		String dataset = "accidents";
		String input = dataset + ".txt";
		String input2 = dataset + "Invest.txt";

		Stock stock = new Stock();
		stock.loadFile(input2);

		int kStart = 20000;
		int decreaseStep = 2000;
		int repeat = 5;
		int dbSize = Integer.MAX_VALUE;

		int currentK = kStart;

		for (int run = 1; run <= repeat; run++) {

			System.out.println("\n---------------- RUN " + run + " ----------------");
			System.out.println("Current k = " + currentK);

			// === FILE OUTPUT NAME ===
			String outFileName = "results/"+dataset + "_k" + currentK + ".txt";

			// === Prepare writer ===
			PrintWriter writer = new PrintWriter(new FileWriter(outFileName));

			writer.println("=========== EXPERIMENT RUN " + run + " ===========");
			writer.println("Dataset: " + dataset);
			writer.println("k = " + currentK);
			writer.println("------------------------------------------");

			AlgoTKHEP algo = new AlgoTKHEP();

			long start = System.currentTimeMillis();
			algo.runAlgorithm(currentK, input, null, stock, true, dbSize, true);
			long end = System.currentTimeMillis();
			long runtime = end - start;

			// ===== WRITE RUNTIME =====
			writer.println("Runtime: " + runtime + " ms (" + (runtime / 1000.0) + " s)");

			// ===== WRITE STATS =====
			writer.println("\n========== TKHEP - STATS ============");
			writer.println("minEfficiency = " + algo.minEfficiency);
			writer.println("High Efficiency itemsets count = " + algo.globalK);
			writer.println("Total time = " + ((algo.endTimestamp - algo.startTimestamp) / 1000.0) + " s");
			writer.println("Transaction merge count = " + algo.mergeCount);
			writer.println("Transaction read count = " + algo.transactionReadingCount);
			writer.println("Max memory = " + MemoryLogger.getInstance().getMaxMemory());
			writer.println("Candidate count = " + algo.candidateCount);
			writer.println("=====================================");

			// ===== WRITE RESULTS (nếu muốn) =====
			// algo.printResult(writer);  // <-- nếu bạn có phiên bản printResult(PrintWriter)

			writer.close();

			System.out.println(">> Output saved to: " + outFileName);

			algo = null;
			System.gc();

			currentK -= decreaseStep;
			if (currentK <= 0) break;
		}

		System.out.println("\n============== FINISHED ALL RUNS ==============");
	}
}
