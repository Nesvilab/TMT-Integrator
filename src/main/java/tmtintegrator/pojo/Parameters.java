/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package tmtintegrator.pojo;

import tmtintegrator.constants.Constants;

import java.io.File;
import java.util.*;

public class Parameters {
    public String reportPath;
    public File combinedF;
    public int channelNum;
    public List<Subplex> subplexes = new ArrayList<>();
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
    public boolean isNglyco = false;
    public boolean useGlycoComposition = false;
    public float glycoQval = -1;
    public boolean log2transformed = true;
    public float[] labels;
    public int minResolution = Constants.DEFAULT_MIN_RESOLUTION;
    public float minSNR = Constants.DEFAULT_MIN_SNR;
    public boolean addIsobaricFilter = false;

    public List<File> fileList = new ArrayList<>();
    public List<File> proteinFileList = new ArrayList<>();
    public List<String> fNameLi = new ArrayList<>();
    public Set<String> allGeneSet = new HashSet<>();
    public Map<String, String> titleMap = new HashMap<>();
    public Map<String, Index> indMap = new HashMap<>();
    public Map<String, ProteinIndex> proteinIndexMap = new HashMap<>();
}
