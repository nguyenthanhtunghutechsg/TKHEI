package Stock_Package;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class GenInvest {
    Map<Integer,Integer> itemToInvest=new HashMap<>();
    Set<Integer> itemSet=new HashSet<>();
    public static void main(String[] args) throws IOException {
        String input = "SUSY.txt";
        String output = "SUSYInvest.txt";
        GenInvest genInvest=new GenInvest();
        genInvest.read(input);

        genInvest.testRandom();
        genInvest.writeFile(output);
    }

    private void writeFile(String output) throws IOException {
        FileWriter fw = new FileWriter(output);
        for (Integer item:itemToInvest.keySet()) {
            fw.write(item+" "+itemToInvest.get(item));
            fw.write("\n");
            fw.flush();
        }
    }

    public void testRandom() {
        Random random = new Random();
        for (Integer item:itemSet) {
            int num = (int)Math.round(10000+random.nextGaussian()*10000);
            while (num<=0){
                num = (int)Math.round(100+random.nextGaussian()*25);
            }
            itemToInvest.put(item,num);
        }

    }

    public void read(String input) throws IOException {

        String thisline0;
        BufferedReader br0 = new BufferedReader(new FileReader(input));
        while ((thisline0 = br0.readLine()) != null) {
            String[] Trans = thisline0.split(":");
            String[] items = Trans[0].split(" ");
            for (String item:items) {
                Integer itemInt=Integer.valueOf(item);
                itemSet.add(itemInt);
            }
        }
    }
}
