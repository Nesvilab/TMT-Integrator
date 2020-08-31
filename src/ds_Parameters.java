import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

public class ds_Parameters
{
    File fastaF;
    String reportPath;
    File combinedF;
    int channelNum;
    String refTag;
    int groupBy;
    boolean psmNorm;
    boolean outlierRemoval;
    int protNorm;
    float minPepProb;
    float minPurity;
    float minPercent;
    boolean uniquePep;
    int uniqueGene;
    boolean bestPsm;
    String[] protExcludeAry;
    boolean allow_overlabel;
    boolean allow_unlabeled;
    List<String> modTagLi = new ArrayList<String>();
    float minSiteProb;
    boolean ms1Int;
    boolean  top3Pep;
    boolean print_RefInt;
    int add_Ref=-1;
    double max_pep_prob_thres=0;
    int min_ntt=0;
    int abn_type=0; //0: ratio-based; 1: raw-based

    List<File> FileLi = new ArrayList<File>();
    List<String> fNameLi = new ArrayList<String>();
    TreeMap<String, String> fastaMap = new TreeMap<String, String>(); //Key: z; Value: ds_Peak
    ds_ColumnIndex ci = new ds_ColumnIndex();
    List<String> AllGeneLi = new ArrayList<String>();
    TreeMap<String, String> TitleMap = new TreeMap<String, String>();
}
