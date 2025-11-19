package Tool_Package;
public class Tool {
    public static int compare(double double1, double double2)  {
        double diff=double1-double2;
        if (diff>0.00001){
            return 1;
        } else if (diff<-0.00001) {
            return -1;
        }else {
            return 0;
        }
        //BigDecimal data1 = new BigDecimal(double1);
        //BigDecimal data2 = new BigDecimal(double2);
        //return data1.compareTo(data2);
    }
}
