package Stock_Package;
import java.io.*;
import java.util.HashMap;
import java.util.Map;

public class Stock {
   // Map<Integer,Integer> profitMap=new HashMap<>();
    public Map<Integer,Integer> investMap=new HashMap<>();
    public void loadFile(String path) throws IOException {
        String thisLine;
        BufferedReader myInput = null;
        try {
            FileInputStream fin = new FileInputStream(new File(path));
            myInput = new BufferedReader(new InputStreamReader(fin));
            // for each transaction (line) in the input file
            while ((thisLine = myInput.readLine()) != null) {
                // if the line is  a comment, is  empty or is a
                // kind of metadata
                if (thisLine.isEmpty() == true ||
                        thisLine.charAt(0) == '#' || thisLine.charAt(0) == '%'
                        || thisLine.charAt(0) == '@') {
                    continue;
                }

                // process the transaction
              String[] str=thisLine.split(" ");
              Integer item=Integer.valueOf(str[0]);
              Integer invest=Integer.valueOf(str[1]);
              //profitMap.put(item,profit);
              investMap.put(item,invest);
            }
        } catch (Exception e) {
            // catch exceptions
            e.printStackTrace();
        }finally {
            if(myInput != null){
                // close the file
                myInput.close();
            }
        }
    }

}
