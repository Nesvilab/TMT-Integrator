package tmtintegrator;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

public class ds_Parameters
{
    File fastaF; // FIXME 03: no usage
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
    float minSiteProb;
    boolean ms1Int;
    boolean  top3Pep;
    boolean print_RefInt;
    int add_Ref=-1;
    double max_pep_prob_thres=0;
    int min_ntt=0;
    int abn_type=0; //0: ratio-based; 1: raw-based
    int aggregation_method=0; //psm aggregation, 0: median; 1: weighted-ratio
    String modAA="";
    String columntag="";
    List<String> modTagLi = new ArrayList<String>();

    boolean geneflag = false;
    boolean glycoflag = false;
    boolean useGlycoComposition = false;    // whether to use glycan composition string or glycan mass for indexing
    float glycoQval = -1;

    String prefix = "rev_";
    boolean log2transformed = true;

    List<File> FileLi = new ArrayList<File>();
    List<String> fNameLi = new ArrayList<String>();
    //TreeMap<String, String> fastaMap = new TreeMap<String, String>(); //Key: z; Value: ds_Peak
    List<String> AllGeneLi = new ArrayList<String>();
    TreeMap<String, String> TitleMap = new TreeMap<String, String>();
    TreeMap<String, ds_Index> indMap = new TreeMap<String, ds_Index>();
    //ds_ColumnIndex ci = new ds_ColumnIndex();

    TreeMap<String, String> ppMap = new TreeMap<String, String>(); //proteinID-protein map
    //TreeMap<String, String> phMap = new TreeMap<String, String>(); //protein-fasta header map
}
