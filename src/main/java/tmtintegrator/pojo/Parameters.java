package tmtintegrator.pojo;

import java.io.File;
import java.util.*;

public class Parameters {
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
    public boolean top3Pep;
    public boolean print_RefInt;
    public int add_Ref = -1;
    public double max_pep_prob_thres = 0;
    public int min_ntt = 0;
    public int abn_type = 0; //0: ratio-based; 1: raw-based
    public int aggregation_method = 0; //psm aggregation, 0: median; 1: weighted-ratio
    public String modAA = "";
    public String columntag = "";
    public Set<String> modTagSet = new HashSet<>();
    public boolean geneflag = false;
    public boolean glycoflag = false;
    public boolean useGlycoComposition = false;    // whether to use glycan composition string or glycan mass for indexing
    public float glycoQval = -1;
    public boolean log2transformed = true;

    public List<File> fileList = new ArrayList<>();
    public List<File> proteinFileList = new ArrayList<>();
    public List<String> fNameLi = new ArrayList<>();
    public Set<String> allGeneSet = new HashSet<>();
    public Map<String, String> titleMap = new HashMap<>();
    public Map<String, Index> indMap = new HashMap<>();
    public Map<String, ProteinIndex> proteinIndexMap = new HashMap<>();
}
