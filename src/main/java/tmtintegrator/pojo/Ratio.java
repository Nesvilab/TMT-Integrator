package tmtintegrator.pojo;

public class Ratio implements Comparable<Ratio>
{
    public double preInt = 0;
    public double rt = 0;
    public double ratio = 0;
    public double weight = 0;

    @Override
    public int compareTo(Ratio r) {
        return r.toString().compareTo(r.toString());
    }
}
