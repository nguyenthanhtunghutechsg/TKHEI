package BaseLine;


import java.io.IOException;

import Stock_Package.Stock;

// dEFIM TESTER, OUTPUT TO SCREEN
public class main_TKHEP {

	public static void main(String[] arg) throws IOException {
		String input = "DB_Utility.txt";
		String output = "output.txt";
		String input2 = "Stock.txt";
		Stock stock = new Stock();
		stock.loadFile(input2);
		int k = 50;
		int dbSize =  67557*25/100;
		AlgoTKHEP algo = new AlgoTKHEP(); // Create the dEFIM algorithm object

		// execute the algorithm
		Itemsets itemsets = algo.runAlgorithm(k, input, output, stock, true, dbSize, true);

		algo.printStats(); // Print statistics
		algo.printResult();
		// itemsets.printItemsets(); // Print the itemsets
	}
}
