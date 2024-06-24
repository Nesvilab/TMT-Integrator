package tmtintegrator.pojo;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

public class Parameters
{
    public File fastaF; // FIXME 03: no usage
    public String reportPath;
    public File combinedF;
    public int channelNum;
    public String refTag;
    public int groupBy;
    public boolean psmNorm;
    public boolean outlierRemoval;
    public int protNorm;
    public float minPepProb;
    public float minPurity;
    public float minPercent;
    public boolean uniquePep;
    public int uniqueGene;
    public boolean bestPsm;
    public String[] protExcludeAry;
    public boolean allow_overlabel;
    public boolean allow_unlabeled;
    public float minSiteProb;
    public boolean ms1Int;
    public boolean  top3Pep;
    public boolean print_RefInt;
    public int add_Ref=-1;
    public double max_pep_prob_thres=0;
    public int min_ntt=0;
    public int abn_type=0; //0: ratio-based; 1: raw-based
    public int aggregation_method=0; //psm aggregation, 0: median; 1: weighted-ratio
    public String modAA="";
    public String columntag="";
    public List<String> modTagLi = new ArrayList<>();

    public boolean geneflag = false;
    public boolean glycoflag = false;
    public boolean useGlycoComposition = false;    // whether to use glycan composition string or glycan mass for indexing
    public float glycoQval = -1;

    public boolean log2transformed = true;

    public List<File> FileLi = new ArrayList<>();
    public List<String> fNameLi = new ArrayList<>();
    //TreeMap<String, String> fastaMap = new TreeMap<String, String>(); //Key: z; Value: ds_Peak
    public List<String> AllGeneLi = new ArrayList<>();
    public TreeMap<String, String> TitleMap = new TreeMap<>();
    public TreeMap<String, Index> indMap = new TreeMap<>();
    //ds_ColumnIndex ci = new ds_ColumnIndex();

    public TreeMap<String, String> ppMap = new TreeMap<>(); //proteinID-protein map
    //TreeMap<String, String> phMap = new TreeMap<String, String>(); //protein-fasta header map
}
