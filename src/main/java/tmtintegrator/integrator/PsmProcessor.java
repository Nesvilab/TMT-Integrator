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
import tmtintegrator.pojo.psm.Psm;
import tmtintegrator.pojo.psm.PsmRecord;
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
    private final Map<String, Map<String, PsmInfo>> groupPsmMap; // <groupKey(report index), <fileName, psmInfo>>
    private final Map<String, Map<String, double[]>> groupAbundanceMap; // <groupKey(report index), <fileName, abundance>>

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
     * Group PSMs by report level into a map.
     *
     * @param psmList list of PSMs
     */
    public void groupPsm(List<Psm> psmList) {
        for (Psm psm : psmList) {
            String filePath = psm.getPsmFile().getAbsolutePath();
            for (PsmRecord psmRecord : psm.getPsmRecords()) {
                groupPsmEntry(filePath, psmRecord);
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
                if (psmInfo.psmRecords.size() >= Constants.PSM_NUM_THRESHOLD) {
                    removeOutlierByChannel(psmInfo.psmRecords, index);
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
                double[] medianValues = new double[index.usedChannelNum + 1]; // one additional for total reference intensity

                // get max peptide probability and all peptide sequences
                maxPeptideProb = updateMaxPepProbAndProteinMap(psmInfo, maxPeptideProb, geneSet, proteinMap);

                // take median ratios and update abundance map
                takeMedianRatios(medianValues, psmInfo, index);
                medianValues[medianValues.length - 1] = psmInfo.totalRefInt; // add total reference intensity
                fileAbundanceMap.put(filePath, medianValues);
                if (groupBy == GroupBy.GENE || groupBy == GroupBy.PROTEIN_ID) {
                    numPsm += psmInfo.psmRecords.size();
                }
            }

            String groupKey = groupEntry.getKey();
            int[] lr = extractLRIndex(groupKey);

            String globalGenePepSeq = createGlobalGenePepSeq(geneSet, proteinMap, numPsm, lr[0], lr[1]);
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
        int location = 5; // last index of group key parts

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
    private void groupPsmEntry(String filePath, PsmRecord psmRecord) {
        String[] lr = Utils.getLRStrings(psmRecord.getExtendedPeptide());
        Map<String, PsmInfo> fileMap = groupPsmMap.computeIfAbsent(psmRecord.getGroupKey(), k -> new HashMap<>());
        PsmInfo psmInfo = fileMap.computeIfAbsent(filePath, k -> new PsmInfo());
        psmInfo.gene = psmRecord.getGene();
        psmInfo.addP(psmRecord.getPeptide(), lr[0], lr[1], psmRecord.getPepsIndex());
        psmInfo.psmRecords.add(psmRecord);
    }

    private void computeTotalRefInt() {
        for (Map.Entry<String, Map<String, PsmInfo>> groupEntry : groupPsmMap.entrySet()) {
            Map<String, PsmInfo> fileMap = groupEntry.getValue();
            for (PsmInfo psmInfo : fileMap.values()) {
                computeTotalRefInt(psmInfo);
            }
        }
    }

    private void computeTotalRefInt(PsmInfo psmInfo) {
        List<Double> refIntensityList = new ArrayList<>();
        for (PsmRecord psmRecord : psmInfo.psmRecords) {
            double refIntensity = psmRecord.getRefIntensity();
            refIntensityList.add(refIntensity);
            psmInfo.totalRefInt += refIntensity; // sum all peptide intensities
        }
        Collections.sort(refIntensityList);

        // compute top 3 intensive peptides
        if (parameters.top3Pep && psmInfo.psmRecords.size() >= 3) {
            psmInfo.totalRefInt = 0;
            for (int i = refIntensityList.size() - 1; i >= refIntensityList.size() - 3; i--) {
                psmInfo.totalRefInt += refIntensityList.get(i);
            }
        }
    }

    private void removeOutlierByChannel(List<PsmRecord> psmRecords, Index index) {
        for (int j = 0; j < index.usedChannelNum; j++) {
            List<Double> ratios = new ArrayList<>();
            for (PsmRecord psmRecord : psmRecords) {
                double ratio = psmRecord.getChannels().get(j);
                if (!Double.isNaN(ratio)) {
                    ratios.add(ratio);
                }
            }

            // calculate IQR for outlier removal
            if (ratios.size() >= Constants.PSM_NUM_THRESHOLD) {
                double[] iqrBounds = Utils.computeIQR(ratios); // lower and upper bounds
                for (PsmRecord psmRecord : psmRecords) {
                    // remove outliers
                    double ratio = psmRecord.getChannels().get(j);
                    if (!Double.isNaN(ratio) && (ratio < iqrBounds[0] || ratio > iqrBounds[1])) {
                        psmRecord.getChannels().set(j, Double.NaN);
                    }
                }
            }
        }
    }

    private double updateMaxPepProbAndProteinMap(PsmInfo psmInfo, double maxPeptideProb,
                                                 Set<String> geneSet, Map<String, Set<String>> proteinMap) {
        for (PsmRecord psmRecord : psmInfo.psmRecords) {
            double probability = psmRecord.getProbability();
            String proteinId = psmRecord.getProteinId();
            // update max peptide probability
            maxPeptideProb = Math.max(maxPeptideProb, probability);
            // add gene to gene list
            geneSet.add(psmInfo.gene);

            // add peptide index to protein map
            Set<String> peptideIndices = proteinMap.computeIfAbsent(proteinId, k -> new HashSet<>());
            peptideIndices.addAll(psmInfo.getPeptideIndex());
        }
        return maxPeptideProb;
    }

    private void takeMedianRatios(double[] medianValues, PsmInfo psmInfo, Index index) {
        for (int j = 0; j < index.usedChannelNum; j++) {
            if (parameters.aggregation_method == 0) {
                List<Double> channelsValues = new ArrayList<>();
                for (PsmRecord psmRecord : psmInfo.psmRecords) {
                    double ratio = psmRecord.getChannels().get(j);
                    if (!Double.isNaN(ratio)) {
                        channelsValues.add(ratio);
                    }
                }
                medianValues[j] = Utils.takeMedian(channelsValues);
            } else if (parameters.aggregation_method == 1) {
                List<Ratio> channelsValues = new ArrayList<>();
                for (PsmRecord psmRecord : psmInfo.psmRecords) {
                    double ratio = psmRecord.getChannels().get(j);
                    if (!Double.isNaN(ratio)) {
                        Ratio ratioObj = new Ratio();
                        ratioObj.preInt = psmRecord.getMs1Intensity();
                        ratioObj.rt = psmRecord.getRetention();
                        ratioObj.ratio = ratio;
                        channelsValues.add(ratioObj);
                    }
                }
                medianValues[j] = Utils.takeWeightedMedian(channelsValues);
            }
        }
    }

    private int[] extractLRIndex(String groupKey) {
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
        return new int[]{l, r};
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
                if (abc == null) { // FIXME 02: this only works for the first one
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

                for (int i = 1; i < parts.length; i++) {
                    String newGroupKey = indexParts[0] + "%" +
                            indexParts[location].charAt(indexParts[location].indexOf(parts[i]) - 1) + parts[i];

                    // update keyMap
                    List<String> keyList = keyMap.computeIfAbsent(newGroupKey, k -> new ArrayList<>());
                    if (!keyList.contains(groupKey)) {
                        keyList.add(groupKey);
                    }
                }
            }
        }
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
                double[] fileMedians = new double[index.usedChannelNum + 1];
                for (int i = 0; i < index.usedChannelNum + 1; i++) {
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
