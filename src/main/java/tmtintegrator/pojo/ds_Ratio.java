package tmtintegrator.pojo;

public class ds_Ratio implements Comparable<ds_Ratio>
{
    public double preInt = 0;
    public double rt = 0;
    public double ratio = 0;
    public double weight = 0;

    @Override
    public int compareTo(ds_Ratio r) {
        return r.toString().compareTo(r.toString());
    }
}
