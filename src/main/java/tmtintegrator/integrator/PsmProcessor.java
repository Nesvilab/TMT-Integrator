package tmtintegrator.integrator;

import tmtintegrator.constants.GroupBy;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.PsmInfo;
import tmtintegrator.pojo.Ratio;
import tmtintegrator.utils.Utils;

import java.util.*;

/**
 * Utility class to process PSMs.
 *
 * @author rogerli on 05/2024
 */

public class PsmProcessor {
    private final Parameters parameters;
    private final GroupBy groupBy;
    private final Map<String, Map<String, PsmInfo>> groupPsmMap; // <groupKey(proteinId), <fileName, psmInfo>>
    private final Map<String, Map<String, double[]>> groupAbundanceMap; // <groupKey(proteinId), <fileName, abundance>>
    public static final int PSM_NUM_THRESHOLD = 4; // threshold for outlier removal

    public PsmProcessor(Parameters parameters, GroupBy groupBy) {
        this.parameters = parameters;
        this.groupPsmMap = new TreeMap<>();
        this.groupAbundanceMap = new TreeMap<>();
        this.groupBy = groupBy;
    }

    public Map<String, Map<String, PsmInfo>> getGroupPsmMap() {
        return groupPsmMap;
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
                if (psmInfo.psmList.size() >= PSM_NUM_THRESHOLD) {
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
            Map<String, double[]> fileAbundanceMap = new TreeMap<>(); // TODO: HashMap?
            Map<String, List<String>> proteinMap = new TreeMap<>(); // <proteinId, peptideIndices> TODO: HashMap?
            int numPsm = 0;
            double maxPeptideProb = 0;
            List<String> geneList = new ArrayList<>();

            for (Map.Entry<String, PsmInfo> entry : fileMap.entrySet()) {
                String filePath = entry.getKey();
                PsmInfo psmInfo = entry.getValue();
                Index index = parameters.indMap.get(filePath);
                double[] medianValues = new double[index.totLen];
                double[][] ratio2DValues = new double[psmInfo.psmList.size()][index.plexNum];
                Ratio[][] ratioObj2DValues = new Ratio[psmInfo.psmList.size()][index.plexNum];

                // get max peptide probability and all peptide sequences
                maxPeptideProb = updateMaxPepProbAndProteinMap(psmInfo, index, maxPeptideProb, geneList, proteinMap);

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

            String globalGenePepSeq = createGlobalGenePepSeq(geneList, proteinMap, numPsm);
            String groupKey = groupEntry.getKey();
            groupKey = updateGroupKey(groupKey, globalGenePepSeq, maxPeptideProb);
            groupAbundanceMap.put(groupKey, fileAbundanceMap);
        }
    }

    // region helper methods
    private void groupPsmEntry(String filePath, String psm, Index index) {
        String[] fields = psm.split("\t");
        String[] psmKeyParts = fields[0].split("#");
        String groupKey = psmKeyParts[0]; // generated by groupBy option while loading
        String newPepSequence = psmKeyParts[2];
        String extPepSequence = Utils.refineExtendedSequence(fields[index.extpepIndex]);
        String gene = fields[index.genecIndex];
        int pepsIndex = Integer.parseInt(fields[index.protsIndex]) - 1;

        // Ensure the group key exists, and create a new PSM info if necessary
        Map<String, PsmInfo> fileMap = groupPsmMap.computeIfAbsent(groupKey, k -> new TreeMap<>());
        PsmInfo psmInfo = fileMap.computeIfAbsent(filePath, k -> new PsmInfo());
        psmInfo.gene = gene;
        if (psmInfo.peptide.isEmpty()) {
            // for new PsmInfo, set peptide and pepsIndex
            psmInfo.peptide = newPepSequence;
            psmInfo.pepsIndex = pepsIndex;
        } else if (psmInfo.peptide.length() < newPepSequence.length()) {
            // update peptide and pepsIndex if new peptide is longer
            psmInfo.peptide = newPepSequence;
//            psmInfo.pepsIndex = pepsIndex; // FIXME: should be updated, but it will always be the first pepsIndex
            // FIXME: Here is the original code, where seems the intention is to update pepsIndex as well
            //   however, with the current implementation, it will always be the first pepsIndex
            //   psmInfo.peptide = (psmInfo.peptide.length() < newPepSequence.length()) ? newPepSequence : psmInfo.peptide;
            //   psmInfo.pepsIndex = (psmInfo.peptide.length() < newPepSequence.length()) ? pepsIndex: psmInfo.pepsIndex; // pepIndex will never be updated
            // FIXME END
        }
        psmInfo.extpep = extPepSequence;
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
                if (row[j] != -9999) { // FIXME: !Double.isNaN(row[j]) is better
                    ratios.add(row[j]);
                }
            }

            // calculate IQR for outlier removal
            if (ratios.size() >= PSM_NUM_THRESHOLD) {
                double[] iqrBounds = Utils.computeIQR(ratios); // lower and upper bounds
                for (double[] row : ratio2DValues) {
                    // remove outliers
                    if (!Double.isNaN(row[j]) && (row[j] < iqrBounds[0] || row[j] > iqrBounds[1])) {
                        row[j] = -9999; // FIXME: Double.NaN is better
                    }
                }
            }
        }
    }

    private double updateMaxPepProbAndProteinMap(PsmInfo psmInfo, Index index, double maxPeptideProb,
                                                 List<String> geneList, Map<String, List<String>> proteinMap) {
        for (String psm : psmInfo.psmList) {
            String[] fields = psm.split("\t");
            double peptideProb = Double.parseDouble(fields[index.pepProbcIndex]);
            String proteinId = fields[index.proteinIDcIndex];
            // update max peptide probability
            maxPeptideProb = Math.max(maxPeptideProb, peptideProb);
            // add gene to gene list
            if (!geneList.contains(psmInfo.gene)) {
                geneList.add(psmInfo.gene);
            }

            // add peptide index to protein map
            String peptideIndex = psmInfo.peptide + "@" + psmInfo.pepsIndex + "@" + psmInfo.extpep;
            List<String> peptideIndices = proteinMap.computeIfAbsent(proteinId, k -> new ArrayList<>());
            if (!peptideIndices.contains(peptideIndex)) {
                peptideIndices.add(peptideIndex); // FIXME: potential bug, remove this fixme once tested.
            }
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
                    ratio2DValues[i][j - index.abnIndex] = -9999; // FIXME: Double.NaN is better
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
                    if (ratios[j] != -9999) { // FIXME: !Double.isNaN(ratios[j]) is better
                        channelValues.add(ratios[j]);
                    }
                }
                medianValues[j] = Utils.takeMedian(channelValues);
            } else if (parameters.aggregation_method == 1) {
                List<Ratio> channelValues = new ArrayList<>();
                for (Ratio[] ratios : ratioObj2DValues) {
                    if (ratios[j].ratio != -9999) { // FIXME: !Double.isNaN(ratios[j].ratio) is better
                        channelValues.add(ratios[j]);
                    }
                }
                medianValues[j] = Utils.takeWeightedMedian(channelValues);
            }
        }
    }

    private String createGlobalGenePepSeq(List<String> geneList, Map<String, List<String>> proteinMap, int numPsm) {
        String globalGenePepSeq = String.join(";", geneList);
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
                globalGenePepSeq = handleMultiSiteAndMass(proteinMap, globalGenePepSeq, proteinIdSeq);
                break;
            default:
                globalGenePepSeq = createDefaultGenePepSeq(proteinMap, globalGenePepSeq, proteinIdSeq);
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
    private String handleMultiSiteAndMass(Map<String, List<String>> proteinMap, String globalGenePepSeq,
                                          String proteinIdSeq) {
        Map<String, String> map = new TreeMap<>(); // <peptide, extendedPeptide> // TODO: HashMap?
        String peptideIdxSeq = "";

        for (Map.Entry<String, List<String>> entry : proteinMap.entrySet()) {
            int lowestStart = Integer.MAX_VALUE;
            int highestEnd = Integer.MIN_VALUE;

            List<String> peptideList = entry.getValue();
            for (String peptide : peptideList) {
                String[] parts = peptide.split("@");
                int pepIndex = Integer.parseInt(parts[1]);
                int endIdx = pepIndex + parts[0].length();
                // update start/end if longer peptide is found
                if ((pepIndex < lowestStart) || (endIdx > highestEnd)) {
                    peptideIdxSeq = (pepIndex + 1) + "\t" + (endIdx + 1);
                    lowestStart = Math.min(pepIndex, lowestStart);
                    highestEnd = Math.max(endIdx, highestEnd);
                }

                // process extended peptide and update map
                int firstDotIdx = parts[2].indexOf('.');
                if (firstDotIdx < 0) {
                    throw new IllegalArgumentException("Invalid extended peptide(no '.'): " + parts[2]);
                }
                int lastDotIdx = parts[2].indexOf('.', firstDotIdx + 1);
                map.put(parts[0], parts[2].substring(0, firstDotIdx + 1) + parts[0] + parts[2].substring(lastDotIdx));
            }
        }

        return String.format("%s\t%s\t%s\t%s\t%s", globalGenePepSeq, proteinIdSeq, String.join(";", map.keySet()),
                String.join(";", map.values()), peptideIdxSeq);
    }

    private String createDefaultGenePepSeq(Map<String, List<String>> proteinMap, String globalGenePepSeq,
                                           String proteinIdSeq) {
        Set<String> peptideSet = new TreeSet<>(); // TODO: HashSet?
        String extendedPeptideSeq = "";
        String peptideIdxSeq = "";

        for (Map.Entry<String, List<String>> entry : proteinMap.entrySet()) {
            List<String> peptideList = entry.getValue();
            int lowestStart = Integer.MAX_VALUE;
            int highestEnd = Integer.MIN_VALUE;

            for (String peptide : peptideList) {
                String[] parts = peptide.split("@");
                peptideSet.add(parts[0]);

                int pepIndex = Integer.parseInt(parts[1]);
                int endIdx = pepIndex + parts[0].length();
                // update start/end if longer peptide is found
                if ((pepIndex < lowestStart) || (endIdx > highestEnd)) {
                    peptideIdxSeq = (pepIndex + 1) + "\t" + (endIdx + 1);
                    lowestStart = Math.min(pepIndex, lowestStart);
                    highestEnd = Math.max(endIdx, highestEnd);
                }

                extendedPeptideSeq = parts[2];
            }
        }
        return String.format("%s\t%s\t%s\t%s\t%s", globalGenePepSeq, proteinIdSeq, String.join(";", peptideSet),
                extendedPeptideSeq, peptideIdxSeq);
    }

    private String updateGroupKey(String groupKey, String globalGenePepSeq, double maxPeptideProb) {
        if (globalGenePepSeq.isEmpty()) {
            return groupKey + "\t" + maxPeptideProb;
        }
        return groupKey + "\t" + globalGenePepSeq + "\t" + maxPeptideProb;
    }
    // endregion
}
