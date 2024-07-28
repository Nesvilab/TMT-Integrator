package tmtintegrator.pojo;

import tmtintegrator.constants.GroupBy;
import tmtintegrator.utils.Utils;

import java.util.*;

/**
 * Processed PSM entry class.
 */

public class ProcessedPsmEntry {

    // region PSM file information
    private final Parameters parameters;
    private final Index index;
    private final String psmLine;
    private final String[] fields;
    // endregion

    // region PSM entry Fields
    private boolean isUsed;
    private double referenceIntensity;
    private String peptide;
    private String assignedModification;
    private String ptmLocal;
    private String proteinId;
    private String gene;
    private double ms1Intensity;
    private int pepIndex;
    // endregion

    private static class PhosphoSiteData {
        private int siteLocalCount;
        private final StringBuilder siteLocalPos; // site localized position
        private final Map<Integer, List<String>> siteLocalMassMap; // Key: position; Value: site localized mass
        private final List<String> siteLocalPosList;
        private int siteCount;
        private final Map<Integer, Double> probMap; // Key: position; Value: site probability
        private int startIndex;
        private int endIndex;

        private PhosphoSiteData() {
            siteLocalCount = 0;
            siteLocalPos = new StringBuilder();
            siteLocalMassMap = new TreeMap<>(); // TODO: HashMap?
            siteLocalPosList = new ArrayList<>();
            siteCount = 0;
            probMap = new TreeMap<>(); // TODO: HashMap?
        }
    }

    public ProcessedPsmEntry(Parameters parameters, Index index, String psm) {
        this.parameters = parameters;
        this.index = index;
        this.psmLine = psm;
        this.fields = psm.split("\t");
    }

    // region Getters
    public boolean isUsed() {
        return isUsed;
    }

    public double getReferenceIntensity() {
        return referenceIntensity;
    }

    public String getPeptide() {
        return peptide;
    }

    public String getAssignedModification() {
        return assignedModification;
    }

    public String getPtmLocal() {
        return ptmLocal;
    }

    public String getProteinId() {
        return proteinId;
    }

    public String getGene() {
        return gene;
    }

    public double getMs1Intensity() {
        return ms1Intensity;
    }
    // endregion

    /**
     * Parse PSM entry fields.
     */
    public void parsePsmEntry() {
        isUsed = Boolean.parseBoolean(fields[index.isUsedIndex]);
        referenceIntensity = Double.parseDouble(fields[index.refIndex]);
        peptide = fields[index.pepcIndex];
        assignedModification = fields[index.assignedModcIndex];
        ptmLocal = index.ptmLocalcIndex >= 0 ? fields[index.ptmLocalcIndex] : "";
        proteinId = fields[index.proteinIDcIndex];
        gene = !fields[index.genecIndex].isEmpty() ? fields[index.genecIndex] : proteinId;
        ms1Intensity = Double.parseDouble(fields[index.ms1IntIndex]);
        pepIndex = Integer.parseInt(fields[index.protsIndex]) - 1;
    }

    /**
     * handle the calculation of the number of phospho site modifications and their positions.
     * @param psmList list of PSM entries
     * @param groupBy group by option
     */
    public void analyzePhosphoSites(List<String> psmList, GroupBy groupBy) {
        PhosphoSiteData data = new PhosphoSiteData();
        // get num of phospho sites with probability > threshold and their positions
        data.startIndex = pepIndex;
        data.endIndex = Integer.parseInt(fields[index.proteIndex]) - 1;

        if (parameters.minSiteProb > 0) {
            analyzeWithThreshold(data);
        } else if (parameters.minSiteProb == 0) {
            analyzeWithoutThreshold(data);
        }

        // use MS1 intensity as reference intensity
        if (parameters.ms1Int) {
            referenceIntensity = useMs1Intensity();
        }

        // set up key
        String groupKey = handleByGroup(psmList, data, groupBy);

        // add to PSM list
        if (!groupKey.isEmpty()) {
            psmList.add(groupKey + "#" + referenceIntensity + "#" + peptide + "#" + psmLine);
        }
    }

    // region Helper methods
    private void analyzeWithThreshold(PhosphoSiteData data) {
        countModifiedSites(data);
        populateProbMap(data);
        extractSitePositions(data);
        // update start and end indices
        data.startIndex = Math.max(0, data.startIndex + pepIndex);
        data.endIndex = Math.max(0, data.endIndex + pepIndex);
    }

    private void analyzeWithoutThreshold(PhosphoSiteData data) {
        List<Integer> positions = extractSitePositionsInMod(data);
        Collections.sort(positions);
        for (int position : positions) {
            data.probMap.put(position, 0d);
            data.siteLocalPos.append(peptide.charAt(position - 1)).append(pepIndex + position);
            data.siteLocalPosList.add(peptide.charAt(position - 1) + String.valueOf(pepIndex + position));
        }
        data.siteLocalCount = positions.size();
        data.siteCount = positions.size();
        data.startIndex = positions.get(0) + pepIndex;
        data.endIndex = positions.get(positions.size() - 1) + pepIndex;
    }

    private void countModifiedSites(PhosphoSiteData data) {
        String massTag = "(" + parameters.columntag.substring(parameters.columntag.indexOf(":") + 1);
        int index = assignedModification.indexOf(massTag);
        while (index >= 0) {
            data.siteCount++;
            index = assignedModification.indexOf(massTag, index + 1);
        }
    }

    private void populateProbMap(PhosphoSiteData data) {
        int leftIndex = ptmLocal.indexOf("(");
        int rightIndex = ptmLocal.indexOf(")");
        int gap = 0;
        boolean isFirst = true;

        while (leftIndex >= 0) {
            // extract probability
            double probability = Double.parseDouble(ptmLocal.substring(leftIndex + 1, rightIndex));
            // TODO: key updating logic needs to be reviewed
            int key;
            if (isFirst) {
                key = leftIndex;
                isFirst = false;
            } else {
                gap += (rightIndex - leftIndex + 1);
                key = leftIndex - gap;
            }
            data.probMap.put(key, probability);
            // update indices
            leftIndex = ptmLocal.indexOf("(", leftIndex + 1);
            rightIndex = ptmLocal.indexOf(")", rightIndex + 1);
        }
    }

    private void extractSitePositions(PhosphoSiteData data) {
        boolean isStart = false;
        for (Map.Entry<Integer, Double> entry : data.probMap.entrySet()) {
            int key = entry.getKey();
            double probability = entry.getValue();
            if (probability > parameters.minSiteProb) {
                data.siteLocalCount++;
                // TODO: need add twice?
                data.siteLocalPos.append(peptide.charAt(key - 1)).append(pepIndex + key);
                data.siteLocalPosList.add(peptide.charAt(key - 1) + String.valueOf(pepIndex + key));
            }
            if (!isStart) {
                data.startIndex = key;
                isStart = true;
            }
            data.endIndex = key;
        }
    }

    private List<Integer> extractSitePositionsInMod(PhosphoSiteData data) {
        List<Integer> positions = new ArrayList<>();
        String[] assignedMods = assignedModification.split(",");
        for (String assignedMod : assignedMods) {
            String glycoComp = index.glycoCompositionIndex == -1 ? "" : fields[index.glycoCompositionIndex];
            String mod = Utils.getAssignedModIndex(assignedMod, glycoComp, parameters.useGlycoComposition);

            if (Utils.isModTagMatch(mod, parameters.modTagSet)) {
                // extract position number, ex: 11N(2569.9045) -> 11
                int position = Integer.parseInt(assignedMod.substring(0, assignedMod.indexOf("(") - 1).trim());
                positions.add(position);
                updateSiteLocalMassMap(data, position, mod);
            }
        }
        return positions;
    }

    private void updateSiteLocalMassMap(PhosphoSiteData data, int position, String mod) {
        // simplify the mod string for index: remove parentheses and first AA code
        mod = mod.replace("(", "").replace(")", "").substring(1);
        List<String> modMassList = data.siteLocalMassMap.getOrDefault(position, new ArrayList<>());
        modMassList.add(mod);
        data.siteLocalMassMap.put(position, modMassList);
    }

    private double useMs1Intensity() {
        double sum = 0;
        for (int i = index.abnIndex; i < fields.length; i++) {
            sum += Double.parseDouble(fields[i]);
        }
        return ms1Intensity * (referenceIntensity / sum);
    }

    private String handleByGroup(List<String> psmList, PhosphoSiteData data, GroupBy groupBy) {
        String groupKey = "";
        switch (groupBy) {
            case GENE:
                groupKey = gene;
                break;
            case PROTEIN_ID:
                groupKey = proteinId;
                break;
            case PEPTIDE:
                groupKey = generatePeptideGroupKey(data);
                break;
            case MULTI_PHOSPHO_SITE:
                groupKey = handleMultiPhosphoSites(data);
                break;
            case SINGLE_PHOSPHO_SITE:
                // no group key generated for single phospho site, directly add to psm list
                handleSinglePhosphoSite(data, psmList);
                break;
            case MULTI_MASS_GLYCO:
                groupKey = handleMultiMassGlyco(data);
                break;
        }
        return groupKey;
    }

    private String generatePeptideGroupKey(PhosphoSiteData data) {
        if (parameters.minSiteProb >= 0 && data.siteCount > 0) {
            return proteinId + "%" + data.startIndex + "%" + data.endIndex;
        } else if (parameters.minSiteProb < 0) {
            return proteinId + "%" + peptide;
        }
        return "";
    }

    private String handleMultiPhosphoSites(PhosphoSiteData data) {
        if (data.siteCount <= 0) {
            return "";
        }
        // generate group key
        StringBuilder key = new StringBuilder(proteinId + "%" + data.startIndex + "%" + data.endIndex +
                "%" + data.siteCount + "%" + data.siteLocalCount);
        if (data.siteLocalPos.length() > 0) {
            key.append("%").append(data.siteLocalPos);
        }
        String[] peptideKeys = peptide.split("-");
        if (peptideKeys.length > 1) {
            key.append("%").append(peptideKeys[1]);
        }
        // mark localized sites in peptide sequence
        markLocalizedSites(data);

        return key.toString();
    }

    private void handleSinglePhosphoSite(PhosphoSiteData data, List<String> psmList) {
        if (data.siteLocalPosList.isEmpty()) {
            return;
        }
        markLocalizedSites(data);
        // handle single phospho site without generate group key
        for (String site : data.siteLocalPosList) {
            psmList.add(proteinId + "%" + site + "#" + referenceIntensity + "#" + peptide + "#" + psmLine);
        }
    }

    private String handleMultiMassGlyco(PhosphoSiteData data) {
        if (data.siteCount <= 0) {
            return "";
        }
        // generate site localized mass string
        StringBuilder siteLocalMass = new StringBuilder();
        for (List<String> massList : data.siteLocalMassMap.values()) {
            siteLocalMass.append(String.join("_", massList)).append("_");
        }
        siteLocalMass.deleteCharAt(siteLocalMass.length() - 1); // remove trailing '_'
        // generate group key
        StringBuilder key = new StringBuilder(proteinId + "%" + data.startIndex + "%" + data.endIndex +
                "%" + data.siteCount + "%" + data.siteLocalCount);
        if (data.siteLocalPos.length() > 0) {
            key.append("%").append(data.siteLocalPos);
        }
        key.append("%").append(siteLocalMass);
        // mark localized sites in peptide sequence
        multiMassMarkLocalizedSites(data);

        return key.toString();
    }

    /**
     * Mark localized sites in peptide sequence.
     *
     * @param data PhosphoSiteData object
     */
    private void markLocalizedSites(PhosphoSiteData data) {
        StringBuilder newPeptideSeq = new StringBuilder();
        for (int i = 0; i < peptide.length(); i++) {
            char aa = peptide.charAt(i);
            if (data.probMap.containsKey(i + 1) && data.probMap.get(i + 1) >= parameters.minSiteProb) {
                newPeptideSeq.append(Character.toLowerCase(aa));
            } else {
                newPeptideSeq.append(aa);
            }
        }
        peptide = newPeptideSeq.toString();
    }

    private void multiMassMarkLocalizedSites(PhosphoSiteData data) {
        StringBuilder newPeptideSeq = new StringBuilder();
        for (int i = 0; i < peptide.length(); i++) {
            char aa = peptide.charAt(i);
            if (data.probMap.containsKey(i + 1) && data.probMap.get(i + 1) > parameters.minSiteProb) {
                newPeptideSeq.append(Character.toLowerCase(aa));
            } else {
                newPeptideSeq.append(aa);
            }
        }
        peptide = newPeptideSeq.toString();
    }
    // endregion
}
