package tmtintegrator;

public class ds_Ratio implements Comparable<ds_Ratio>
{
    double preInt = 0;
    double rt = 0;
    double ratio = 0;
    double weight = 0;

    @Override
    public int compareTo(ds_Ratio r) {
        return r.toString().compareTo(r.toString());
    }
}
