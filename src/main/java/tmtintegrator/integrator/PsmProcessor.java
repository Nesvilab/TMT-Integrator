package tmtintegrator.integrator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import tmtintegrator.constants.Constants;
import tmtintegrator.constants.GroupBy;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.PsmInfo;
import tmtintegrator.pojo.Ratio;
import tmtintegrator.utils.Utils;

/**
 * Utility class to process PSMs.
 *
 * @author rogerli on 05/2024
 */

public class PsmProcessor {

    private static final Pattern p = Pattern.compile("[A-Z](\\d+)");
    private final Parameters parameters;
    private final GroupBy groupBy;
    private final Map<String, Map<String, PsmInfo>> groupPsmMap; // <groupKey(proteinId), <fileName, psmInfo>>
    private final Map<String, Map<String, double[]>> groupAbundanceMap; // <groupKey(proteinId), <fileName, abundance>>

    public PsmProcessor(Parameters parameters, GroupBy groupBy) {
        this.parameters = parameters;
        this.groupPsmMap = new HashMap<>();
        this.groupAbundanceMap = new TreeMap<>();
        this.groupBy = groupBy;
    }

    public Map<String, Map<String, double[]>> getGroupAbundanceMap() {
        return groupAbundanceMap;
    }

    /**
     * Group PSMs by protein ID into a map.
     *
     * @param fileMap map of file path to list of PSM entries
     */
    public void groupPsm(Map<String, List<String>> fileMap) {
        for (Map.Entry<String, List<String>> entry : fileMap.entrySet()) {
            String filePath = entry.getKey();
            List<String> psmList = entry.getValue();
            Index index = parameters.indMap.get(filePath);

            for (String psm : psmList) {
                groupPsmEntry(filePath, psm, index);
            }
        }

        // compute the total reference intensity
        computeTotalRefInt();
    }

    /**
     * Remove outliers from PSMs.
     */
    public void removeOutlier() {
        for (Map.Entry<String, Map<String, PsmInfo>> groupEntry : groupPsmMap.entrySet()) {
            Map<String, PsmInfo> fileMap = groupEntry.getValue();
            for (Map.Entry<String, PsmInfo> entry : fileMap.entrySet()) {
                String filePath = entry.getKey();
                PsmInfo psmInfo = entry.getValue();
                Index index = parameters.indMap.get(filePath);
                // Only process psmInfo with psmList over threshold for speed
                if (psmInfo.psmList.size() >= Constants.PSM_NUM_THRESHOLD) {
                    List<String> normPsmList = processPsmListByIQR(psmInfo.psmList, index);
                    psmInfo.psmList.clear();
                    psmInfo.psmList.addAll(normPsmList);
                }
            }
        }
    }

    /**
     * Collapse PSMs into abundance map.
     */
    public void collapse() {
        for (Map.Entry<String, Map<String, PsmInfo>> groupEntry : groupPsmMap.entrySet()) {
            Map<String, PsmInfo> fileMap = groupEntry.getValue();
            Map<String, double[]> fileAbundanceMap = new HashMap<>();
            Map<String, Set<String>> proteinMap = new TreeMap<>();
            int numPsm = 0;
            double maxPeptideProb = 0;
            Set<String> geneSet = new HashSet<>();

            for (Map.Entry<String, PsmInfo> entry : fileMap.entrySet()) {
                String filePath = entry.getKey();
                PsmInfo psmInfo = entry.getValue();
                Index index = parameters.indMap.get(filePath);
                double[] medianValues = new double[index.totLen];
                double[][] ratio2DValues = new double[psmInfo.psmList.size()][index.plexNum];
                Ratio[][] ratioObj2DValues = new Ratio[psmInfo.psmList.size()][index.plexNum];

                // get max peptide probability and all peptide sequences
                maxPeptideProb = updateMaxPepProbAndProteinMap(psmInfo, index, maxPeptideProb, geneSet, proteinMap);

                // get ratio values
                if (parameters.aggregation_method == 0) {
                    populateRatio2DValues(psmInfo, index, ratio2DValues);
                } else if (parameters.aggregation_method == 1) {
                    populateRatioObj2DValues(psmInfo, index, ratioObj2DValues);
                }

                // take median ratios and update abundance map
                takeMedianRatios(index, medianValues, ratio2DValues, ratioObj2DValues);
                medianValues[index.plexNum] = psmInfo.totalRefInt; // add total reference intensity
                fileAbundanceMap.put(filePath, medianValues);
                if (groupBy == GroupBy.GENE || groupBy == GroupBy.PROTEIN_ID) {
                    numPsm += psmInfo.psmList.size();
                }
            }

            String groupKey = groupEntry.getKey();
            int l = 0, r = 0;
            if (!parameters.modTagSet.stream().allMatch(e -> e.equals("none"))) {
                if (groupBy == GroupBy.PEPTIDE || groupBy == GroupBy.MULTI_MASS_GLYCO) {
                    String[] ss = groupKey.split("%");
                    l = Integer.parseInt(ss[1]);
                    r = Integer.parseInt(ss[2]);
                } else if (groupBy == GroupBy.SINGLE_PHOSPHO_SITE || groupBy == GroupBy.MULTI_PHOSPHO_SITE) {
                    l = Integer.MAX_VALUE;
                    r = Integer.MIN_VALUE;
                    String[] ss = groupKey.split("%");
                    Matcher m = p.matcher(ss[ss.length - 1]);
                    while (m.find()) {
                        int t = Integer.parseInt(m.group(1));
                        l = Math.min(l, t);
                        r = Math.max(r, t);
                    }
                }
            }

            String globalGenePepSeq = createGlobalGenePepSeq(geneSet, proteinMap, numPsm, l, r);
            groupKey = updateGroupKey(groupKey, globalGenePepSeq, maxPeptideProb);
            groupAbundanceMap.put(groupKey, fileAbundanceMap);
        }
    }

    /**
     * Generate single site, remove multiple sites if single site exists.
     * Calculate the median abundance and update the group abundance map.
     */
    public void generateSingleSite() {
        // <groupKey, List<groupKey>>
        Map<String, List<String>> keyMap = new HashMap<>();
        int location = 5;

        // Cluster keys based on the index
        clusterKeys(keyMap, location);
        // Remove multiple sites if single site exists
        removeMultiSites(keyMap);
        // Calculate the median abundance
        Map<String, Map<String, double[]>> newGroupAbundanceMap = calculateMedianAbundance(keyMap);
        // Update the group abundance map
        groupAbundanceMap.clear();
        groupAbundanceMap.putAll(newGroupAbundanceMap);
    }

    // region helper methods
    private void groupPsmEntry(String filePath, String psm, Index index) {
        String[] fields = psm.split("\t");
        String[] psmKeyParts = fields[0].split("#");
        String groupKey = psmKeyParts[0]; // generated by groupBy option while loading
        String newPepSequence = psmKeyParts[2];
        String[] lr = Utils.getLRStrings(fields[index.extpepIndex]);
        String gene = fields[index.genecIndex];
        int pepsIndex = Integer.parseInt(fields[index.protsIndex]) - 1;

        // Ensure the group key exists, and create a new PSM info if necessary
        Map<String, PsmInfo> fileMap = groupPsmMap.computeIfAbsent(groupKey, k -> new HashMap<>());
        PsmInfo psmInfo = fileMap.computeIfAbsent(filePath, k -> new PsmInfo());
        psmInfo.gene = gene;
        psmInfo.addP(newPepSequence, lr[0], lr[1], pepsIndex);
        psmInfo.psmList.add(psm);
    }

    private void computeTotalRefInt() {
        for (Map.Entry<String, Map<String, PsmInfo>> groupEntry : groupPsmMap.entrySet()) {
            Map<String, PsmInfo> fileMap = groupEntry.getValue();
            for (PsmInfo psmInfo : fileMap.values()) {
                computePsmInfoTotalRefInt(psmInfo);
            }
        }
    }

    private void computePsmInfoTotalRefInt(PsmInfo psmInfo) {
        List<Double> refIntensityList = new ArrayList<>();
        for (String psm : psmInfo.psmList) {
            double refIntensity = extractRefIntensity(psm);
            refIntensityList.add(refIntensity);
            psmInfo.totalRefInt += refIntensity; // sum all peptide intensities
        }
        Collections.sort(refIntensityList);

        // compute top 3 intensive peptides TODO: why overwriting totalRefInt?
        if (parameters.top3Pep && psmInfo.psmList.size() >= 3) {
            psmInfo.totalRefInt = 0;
            for (int i = refIntensityList.size() - 1; i >= refIntensityList.size() - 3; i--) {
                psmInfo.totalRefInt += refIntensityList.get(i);
            }
        }
    }

    private double extractRefIntensity(String psm) {
        String[] fields = psm.split("\t");
        String[] psmKeyParts = fields[0].split("#");
        return Double.parseDouble(psmKeyParts[1]);
    }

    private List<String> processPsmListByIQR(List<String> psmList, Index index) {
        double[][] ratio2DValues = Utils.convertTo2DArray(psmList, index);
        // remove outlier from each channel
        removeOutlierByChannel(ratio2DValues, index);
        return Utils.updatePsmRatios(psmList, ratio2DValues, index);
    }

    private void removeOutlierByChannel(double[][] ratio2DValues, Index index) {
        for (int j = 0; j < index.plexNum; j++) {
            List<Double> ratios = new ArrayList<>();
            for (double[] row : ratio2DValues) {
                if (!Double.isNaN(row[j])) {
                    ratios.add(row[j]);
                }
            }

            // calculate IQR for outlier removal
            if (ratios.size() >= Constants.PSM_NUM_THRESHOLD) {
                double[] iqrBounds = Utils.computeIQR(ratios); // lower and upper bounds
                for (double[] row : ratio2DValues) {
                    // remove outliers
                    if (!Double.isNaN(row[j]) && (row[j] < iqrBounds[0] || row[j] > iqrBounds[1])) {
                        row[j] = Double.NaN;
                    }
                }
            }
        }
    }

    private double updateMaxPepProbAndProteinMap(PsmInfo psmInfo, Index index, double maxPeptideProb,
                                                 Set<String> geneSet, Map<String, Set<String>> proteinMap) {
        for (String psm : psmInfo.psmList) {
            String[] fields = psm.split("\t");
            double peptideProb = Double.parseDouble(fields[index.pepProbcIndex]);
            String proteinId = fields[index.proteinIDcIndex];
            // update max peptide probability
            maxPeptideProb = Math.max(maxPeptideProb, peptideProb);
            // add gene to gene list
            geneSet.add(psmInfo.gene);

            // add peptide index to protein map
            Set<String> peptideIndices = proteinMap.computeIfAbsent(proteinId, k -> new HashSet<>());
            peptideIndices.addAll(psmInfo.getPeptideIndex());
        }
        return maxPeptideProb;
    }

    private void populateRatio2DValues(PsmInfo psmInfo, Index index, double[][] ratio2DValues) {
        for (int i = 0; i < psmInfo.psmList.size(); i++) {
            String[] fields = psmInfo.psmList.get(i).split("\t");
            for (int j = index.abnIndex; j < fields.length; j++) {
                try {
                    ratio2DValues[i][j - index.abnIndex] = Double.parseDouble(fields[j]);
                } catch (NumberFormatException e) {
                    ratio2DValues[i][j - index.abnIndex] = Double.NaN;
                }
            }
        }
    }

    private void populateRatioObj2DValues(PsmInfo psmInfo, Index index, Ratio[][] ratioObj2DValues) {
        for (int i = 0; i < psmInfo.psmList.size(); i++) {
            String[] fields = psmInfo.psmList.get(i).split("\t");
            for (int j = index.abnIndex; j < fields.length; j++) {
                Ratio ratio = new Ratio();
                ratio.preInt = Double.parseDouble(fields[index.ms1IntIndex]);
                ratio.rt = Double.parseDouble(fields[index.rtIndex]);
                ratio.ratio = Double.parseDouble(fields[j]);
                ratioObj2DValues[i][j - index.abnIndex] = ratio;
            }
        }
    }

    private void takeMedianRatios(Index index, double[] medianValues, double[][] ratio2DValues,
                                  Ratio[][] ratioObj2DValues) {
        for (int j = 0; j < index.plexNum; j++) {
            if (parameters.aggregation_method == 0) {
                List<Double> channelValues = new ArrayList<>();
                for (double[] ratios : ratio2DValues) {
                    if (!Double.isNaN(ratios[j])) {
                        channelValues.add(ratios[j]);
                    }
                }
                medianValues[j] = Utils.takeMedian(channelValues);
            } else if (parameters.aggregation_method == 1) {
                List<Ratio> channelValues = new ArrayList<>();
                for (Ratio[] ratios : ratioObj2DValues) {
                    if (!Double.isNaN(ratios[j].ratio)) {
                        channelValues.add(ratios[j]);
                    }
                }
                medianValues[j] = Utils.takeWeightedMedian(channelValues);
            }
        }
    }

    private String createGlobalGenePepSeq(Set<String> geneSet, Map<String, Set<String>> proteinMap, int numPsm, int indexStart, int indexEnd) {
        String globalGenePepSeq = String.join(";", geneSet);
        String proteinIdSeq = String.join(";", proteinMap.keySet());
        switch (groupBy) {
            case GENE:
                globalGenePepSeq = numPsm + "\t" + proteinIdSeq;
                break;
            case PROTEIN_ID:
                globalGenePepSeq = numPsm + "\t" + globalGenePepSeq;
                break;
            case MULTI_PHOSPHO_SITE:
            case MULTI_MASS_GLYCO:
                globalGenePepSeq = handleMultiSiteAndMass(proteinMap, globalGenePepSeq, proteinIdSeq, indexStart, indexEnd);
                break;
            default:
                globalGenePepSeq = createDefaultGenePepSeq(proteinMap, globalGenePepSeq, proteinIdSeq, indexStart, indexEnd);
        }
        return globalGenePepSeq;
    }

    /**
     * Create global gene peptide sequence for multi-site and multi-mass.
     *
     * @param proteinMap       map of protein ID to list of peptide indices
     * @param globalGenePepSeq global gene peptide sequence
     * @param proteinIdSeq     protein ID sequence
     * @return global gene peptide sequence
     */
    private String handleMultiSiteAndMass(Map<String, Set<String>> proteinMap, String globalGenePepSeq,
                                          String proteinIdSeq, int indexStart, int indexEnd) {
        Set<String> t = new TreeSet<>();
        int start, i, j, len;
        String[] abc = null;

        for (Map.Entry<String, Set<String>> entry : proteinMap.entrySet()) {
            Set<String> peptideList = entry.getValue();
            for (String peptide : peptideList) { // peptide format: peptide@pepIndex@extendedPeptide
                String[] parts = peptide.split("@");
                if (abc == null && (indexStart != Integer.MAX_VALUE || indexEnd != Integer.MIN_VALUE)) {
                    len = parts[0].length();
                    i = parts[2].indexOf('.');
                    j = parts[2].indexOf('.', i + 1);
                    start = Integer.parseInt(parts[1]);
                    abc = ff(parts[2].substring(0, i), parts[2].substring(i + 1, j), parts[2].substring(j + 1), indexStart, indexEnd, start, start + len + 1);
                }
                t.add(parts[0]);
            }
        }

        if (abc == null) {
            return String.format("%s\t%s\t%s\t\t\t", globalGenePepSeq, proteinIdSeq, String.join(";", t));
        } else {
            return String.format("%s\t%s\t%s\t%s\t%d\t%d", globalGenePepSeq, proteinIdSeq, String.join(";", t), abc[0] + "." + abc[1] + "." + abc[2], indexStart, indexEnd);
        }
    }

    private String createDefaultGenePepSeq(Map<String, Set<String>> proteinMap, String globalGenePepSeq, String proteinIdSeq, int indexStart, int indexEnd) {
        Set<String> peptideSet = new TreeSet<>();
        int start, i, j, len;
        String[] abc = null;

        for (Map.Entry<String, Set<String>> entry : proteinMap.entrySet()) {
            Set<String> peptideList = entry.getValue();
            for (String peptide : peptideList) {
                String[] parts = peptide.split("@");
                peptideSet.add(parts[0]);
                if (abc == null) {
                    len = parts[0].length();
                    start = Integer.parseInt(parts[1]);
                    if (parameters.modTagSet.stream().allMatch(e -> e.contentEquals("none"))) {
                        String[] ss = ("________" + parts[2] + "________").split("\\.");
                        abc = new String[]{ss[0].substring(ss[0].length() - 7), ss[1], ss[2].substring(0, 7)};
                        indexStart = start + 1;
                        indexEnd = start + len;
                    } else {
                        i = parts[2].indexOf('.');
                        j = parts[2].indexOf('.', i + 1);
                        abc = ff(parts[2].substring(0, i), parts[2].substring(i + 1, j), parts[2].substring(j + 1), indexStart, indexEnd, start, start + len + 1);
                    }
                }
            }
        }

        if (abc == null) {
            return String.format("%s\t%s\t%s\t%s\t%d\t%d", globalGenePepSeq, proteinIdSeq, String.join(";", peptideSet), "", indexStart, indexEnd);
        } else {
            return String.format("%s\t%s\t%s\t%s\t%d\t%d", globalGenePepSeq, proteinIdSeq, String.join(";", peptideSet), abc[0] + "." + abc[1] + "." + abc[2], indexStart, indexEnd);
        }
    }

    static String[] ff(String s1, String s2, String s3, int indexStart, int indexEnd, int indexFirstDot, int indexSecondDot) {
        if (indexStart - indexFirstDot - 1 < 0 || indexSecondDot - indexEnd < 0) {
            System.err.printf("Error in ff: %s %s %s %d %d %d %d\n", s1, s2, s3, indexStart, indexEnd, indexFirstDot, indexSecondDot);
        }

        String a = "________" + s1;
        String b = s3 + "________";

        String s11, s22;
        if (8 - indexStart + indexFirstDot > 0) {
            s11 = a.substring(a.length() - 8 + indexStart - indexFirstDot) + s2.substring(0, indexStart - indexFirstDot - 1);
        } else {
            s11 = s2.substring(indexStart - indexFirstDot - 8, indexStart - indexFirstDot - 1);
        }

        if (8 - indexSecondDot + indexEnd > 0) {
            s22 = s2.substring(indexEnd - indexFirstDot) + b.substring(0, 8 - indexSecondDot + indexEnd);
        } else {
            s22 = s2.substring(indexEnd - indexFirstDot, indexEnd - indexFirstDot + 7);
        }

        return new String[]{s11, s2.substring(indexStart - indexFirstDot - 1, indexEnd - indexFirstDot), s22};
    }

    private String updateGroupKey(String groupKey, String globalGenePepSeq, double maxPeptideProb) {
        if (globalGenePepSeq.isEmpty()) {
            return groupKey + "\t" + maxPeptideProb;
        }
        return groupKey + "\t" + globalGenePepSeq + "\t" + maxPeptideProb;
    }

    private void clusterKeys(Map<String, List<String>> keyMap, int location) {
        for (String groupKey : groupAbundanceMap.keySet()) {
            String[] keyParts = groupKey.split("\t");
            String[] indexParts = keyParts[0].split("%");

            if (Utils.tryParseInt(indexParts[indexParts.length - 1]) < 0) {
                String[] parts = indexParts[location].split(parameters.modAA);

                // find the peptide start position in the protein sequence
                String pepseq = keyParts[3];
                int firstIdx = findPepStartIndex(pepseq);
                int pepStartIdx = Integer.parseInt(parts[1]) - firstIdx;

                for (int i = 1; i < parts.length; i++) {
                    String newGroupKey = indexParts[0] + "%" +
                            indexParts[location].charAt(indexParts[location].indexOf(parts[i]) - 1) + parts[i];

                    // update keyMap
                    List<String> keyList = keyMap.computeIfAbsent(newGroupKey, k -> new ArrayList<>());
                    if (!keyList.contains(groupKey)) {
                        keyList.add(groupKey);
                    }
                    // update keyPepMap
                    int pepIndex = Integer.parseInt(parts[i]) - pepStartIdx;
                    Map<String, Integer> indexMap = new HashMap<>();
                    indexMap.put(pepseq, pepIndex);
                }
            }
        }
    }

    private int findPepStartIndex(String peptide) {
        for (int i = 0; i < peptide.length(); i++) {
            if (Character.isLowerCase(peptide.charAt(i))) {
                return i;
            }
        }
        return -1;
    }

    private void removeMultiSites(Map<String, List<String>> keyMap) {
        for (Map.Entry<String, List<String>> entry : keyMap.entrySet()) {
            List<String> keyList = entry.getValue();
            if (keyList.size() > 1) {
                // filter out single site list
                List<String> singleSiteList = keyList.stream()
                        .filter(key -> key.split("[\t%]")[4].equalsIgnoreCase("1"))
                        .collect(Collectors.toList());
                if (!singleSiteList.isEmpty()) {
                    entry.setValue(singleSiteList);
                }
            }
        }
    }

    private Map<String, Map<String, double[]>> calculateMedianAbundance(Map<String, List<String>> keyMap) {
        Map<String, Map<String, double[]>> newGroupAbundanceMap = new HashMap<>();
        for (String key : keyMap.keySet()) {
            List<String> keyList = keyMap.get(key);
            if (keyList.size() > 1) {
                // <fileName, abundances>
                Map<String, List<double[]>> groupFileAbnMap = new HashMap<>();
                aggregateAbundance(keyList, groupFileAbnMap);
                Map<String, double[]> updatedAbnMap = calculateUpdatedAbundance(groupFileAbnMap);

                String newKey = createNewGroupKey(key, keyList);
                newGroupAbundanceMap.put(newKey, updatedAbnMap);
            } else {
                String[] keyParts = keyList.get(0).split("\t");
                String newKey = key + "\t" + keyParts[1] + "\t" + keyParts[2] + "\t" + keyParts[3] + "\t"
                        + keyParts[4] + "\t" + keyParts[5] + "\t" + keyParts[6] + "\t" + keyParts[7];
                newGroupAbundanceMap.put(newKey, groupAbundanceMap.get(keyList.get(0)));
            }
        }
        return newGroupAbundanceMap;
    }

    private void aggregateAbundance(List<String> keyList, Map<String, List<double[]>> groupFileAbnMap) {
        for (String groupKey : keyList) {
            Map<String, double[]> fileAbundanceMap = groupAbundanceMap.get(groupKey);
            for (String fileName : parameters.fNameLi) {
                double[] medians = fileAbundanceMap.get(fileName);
                if (medians != null) {
                    List<double[]> abundanceList = groupFileAbnMap.computeIfAbsent(fileName, k -> new ArrayList<>());
                    abundanceList.add(medians);
                }
            }
        }
    }

    private Map<String, double[]> calculateUpdatedAbundance(Map<String, List<double[]>> groupFileAbnMap) {
        Map<String, double[]> updatedAbnMap = new HashMap<>();
        for (String fileName : groupFileAbnMap.keySet()) {
            Index index = parameters.indMap.get(fileName);
            List<double[]> abundanceList = groupFileAbnMap.get(fileName);

            if (abundanceList.size() > 1) {
                double[] fileMedians = new double[index.totLen];
                for (int i = 0; i < index.totLen; i++) {
                    List<Double> channelValues = new ArrayList<>();
                    for (double[] medians : abundanceList) {
                        if (!Double.isNaN(medians[i])) {
                            channelValues.add(medians[i]);
                        }
                    }
                    fileMedians[i] = Utils.takeMedian(channelValues);
                }
                updatedAbnMap.put(fileName, fileMedians);
            } else {
                updatedAbnMap.put(fileName, abundanceList.get(0));
            }
        }
        return updatedAbnMap;
    }

    private String createNewGroupKey(String groupKey, List<String> keyList) {
        String[] keyParts = keyList.get(0).split("\t"); // FIXME 09: buggy logic here, all parts after the first part should be identical, in fact not.
        String gene = keyParts[1];
        String proteinId = keyParts[2];
        String extPepList = keyParts[4];
        String start = keyParts[5];
        String end = keyParts[6];

        double maxPeptideProb = 0;
        Set<String> peptideSet = new TreeSet<>();
        for (String key : keyList) {
            String[] parts = key.split("\t");
            peptideSet.add(parts[3]);
            maxPeptideProb = Math.max(maxPeptideProb, Float.parseFloat(parts[7]));
        }

        return String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", groupKey, gene, proteinId,
                String.join(";", peptideSet), extPepList, start, end, maxPeptideProb);
    }
    // endregion
}
