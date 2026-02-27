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

package tmtintegrator.pojo.psm;

import static tmtintegrator.utils.Utils.matchLabels;

import java.util.*;
import java.util.regex.Matcher;
import java.util.stream.Collectors;

import tmtintegrator.constants.Constants;
import tmtintegrator.constants.GroupBy;
import tmtintegrator.constants.ReferenceType;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.Subplex;
import tmtintegrator.utils.Utils;

/**
 * Represents a PSM row in a PSM.tsv file
 */
public class PsmRecord {

    // region PSM file info
    private final Parameters parameters;
    private final Index index;
    // endregion

    // region psm record fields, following the order in the PSM.tsv file
    private String spectrum;
    private String peptide;
    private String modifiedPeptide;
    private String extendedPeptide;
    private int charge;
    private double retention;
    private double calculatedPepMass;
    private double probability;
    private int numEnzymaticTermini;
    private int proteinStart;
    private int proteinEnd;
    private double ms1Intensity;
    private String assignedModifications;
    private String ptmLocalization;
    private double purity;
    private boolean isUnique;
    private String protein;
    private String proteinId;
    private String gene;
    private String mappedGenes;
    // endregion

    // region per-plex channel data
    private List<List<Double>> subplexChannels = new ArrayList<>();
    private List<Double> subplexRefIntensities = new ArrayList<>();
    // endregion

    // region active plex fields (set by setActivePlex)
    private List<Double> channels;
    private double refIntensity;
    // endregion

    // region extra fields
    private double sumTmtIntensity;
    private boolean isUsed;
    private boolean isExcluded;
    private int pepsIndex;
    private String groupKey;
    private PhosphoSiteData phosphoSiteData;
    // endregion

    // region glycan fields
    private double glycanQValue;
    private String glycanComposition;
    // endregion

    // region fields for filtering
    private boolean passResolutionSnr = true;
    private boolean labelFlag;
    private boolean modficationFlag;
    private boolean isPeptideRetained;
    private boolean proteinExclusionFlag;
    private int geneCategory;
    // endregion

    // region backup fields (for groupBy re-iteration)
    private String copyPeptide;
    private List<List<Double>> copySubplexChannels;
    // endregion


    public PsmRecord(Parameters parameters, Index index) {
        this.parameters = parameters;
        this.index = index;
    }

    // region getters

    public String getModifiedPeptide() {
        return modifiedPeptide;
    }

    public String getAssignedModifications() {
        return assignedModifications;
    }

    public String getExtendedPeptide() {
        return extendedPeptide;
    }

    public double getRetention() {
        return retention;
    }

    public double getSumTmtIntensity() {
        return sumTmtIntensity;
    }

    public boolean isUsed() {
        return isUsed;
    }

    public void setUsed(boolean used) {
        isUsed = used;
    }

    public boolean isExcluded() {
        return isExcluded;
    }

    public double getRefIntensity() {
        return refIntensity;
    }

    public List<Double> getChannels() {
        return channels;
    }

    public String getGroupKey() {
        return groupKey;
    }

    public String getGene() {
        return gene;
    }

    public String getPeptide() {
        return peptide;
    }

    public int getPepsIndex() {
        return pepsIndex;
    }

    public double getProbability() {
        return probability;
    }

    public String getProteinId() {
        return proteinId;
    }

    public double getMs1Intensity() {
        return ms1Intensity;
    }

    // endregion

    /**
     * Set the active plex. channels and refIntensity become references to the
     * specified plex's data, so normalization modifies them in-place.
     */
    public void setActivePlex(int plexIdx) {
        channels = subplexChannels.get(plexIdx);
        refIntensity = subplexRefIntensities.get(plexIdx);
    }

    /**
     * Parse a PSM record from a line in a PSM.tsv file
     */
    public void parsePsmRecord(String line, List<Double> tmtIntensities, List<Integer> naChannels) {
        String[] fields = line.split("\t");
        parseRegularFields(fields);
        pepsIndex = proteinStart - 1;
        passResolutionSnr = passResolutionSnr(fields, naChannels);
        parseChannels(fields, tmtIntensities, naChannels);
        updateRefIntensity(fields);
        extractGlycanInfo(fields);
    }

    public void backup() {
        copyPeptide = peptide;
        copySubplexChannels = new ArrayList<>();
        for (List<Double> ch : subplexChannels) {
            copySubplexChannels.add(new ArrayList<>(ch));
        }
    }

    public void reset() {
        isExcluded = false;
        peptide = copyPeptide;
        for (int i = 0; i < subplexChannels.size(); i++) {
            subplexChannels.set(i, new ArrayList<>(copySubplexChannels.get(i)));
        }
    }

    /**
     * Update flags for filtering
     */
    public void updateFlags(Set<String> modTags) {
        labelFlag = updateLabelFlag();
        modficationFlag = updateModificationFlag(modTags);
        isPeptideRetained = !parameters.uniquePep || isUnique;
        proteinExclusionFlag = checkProteinExclusion();
        geneCategory = determineGeneCategory();
    }

    /**
     * Check if the PSM entry passes the criteria
     */
    public boolean isPassCriteria(double tmtThreshold) {
        boolean passAllRefs = true;
        for (double ref : subplexRefIntensities) {
            if (ref <= 0) {
                passAllRefs = false;
                break;
            }
        }
        return passResolutionSnr &&
            purity >= parameters.minPurity &&
            probability >= parameters.minPepProb &&
            sumTmtIntensity >= tmtThreshold &&
            labelFlag &&
            modficationFlag &&
            isPeptideRetained &&
            proteinExclusionFlag &&
            geneCategory >= parameters.uniqueGene &&
            numEnzymaticTermini >= parameters.min_ntt &&
            passAllRefs;
    }

    public String generatePsmKey() {
        String fileName = spectrum.substring(spectrum.lastIndexOf("_") + 1, spectrum.indexOf("."));
        return fileName + "_" + peptide + "_" + charge + "_" + calculatedPepMass;
    }

    public void analyzeByGroup(GroupBy groupBy) {
        phosphoSiteData = new PhosphoSiteData();
        groupKey = "";

        phosphoSiteData.startIndex = pepsIndex;
        phosphoSiteData.endIndex = proteinEnd - 1;

        if (parameters.minSiteProb > 0) {
            processSitesWithProbT();
        } else if (parameters.minSiteProb == 0) {
            processSitesWithoutProbT();
        }

        generateGroupKey(groupBy);

        if (groupBy != GroupBy.SINGLE_PHOSPHO_SITE && groupKey.isEmpty()) {
            isExcluded = true;
        }
    }

    public void useMS1Intensity() {
        for (int i = 0; i < subplexRefIntensities.size(); i++) {
            double ref = subplexRefIntensities.get(i);
            subplexRefIntensities.set(i, ms1Intensity * (ref / sumTmtIntensity));
        }
    }

    // region helper methods
    private void updateRefIntensity(String[] fields) {
        ReferenceType refType = ReferenceType.fromValue(parameters.add_Ref);
        if (refType == ReferenceType.REAL || refType == ReferenceType.RAW_ABUNDANCE) {
            for (Index.SubplexIndex si : index.subplexIndices) {
                subplexRefIntensities.add(Double.parseDouble(fields[si.refIndex]));
            }
        } else {
            for (Subplex s : parameters.subplexes) {
                s.refTag = "Virtual_Reference";
            }
            computeVirtualReference(refType);
        }
    }

    private void parseRegularFields(String[] fields) {
        spectrum = fields[index.spectrumIndex];
        peptide = fields[index.peptideIndex];
        extendedPeptide = fields[index.extpepIndex];
        charge = Integer.parseInt(fields[index.chargeIndex]);
        retention = Double.parseDouble(fields[index.rtIndex]);
        calculatedPepMass = Double.parseDouble(fields[index.pepMassIndex]);
        probability = Double.parseDouble(fields[index.pepProbcIndex]);
        numEnzymaticTermini = Integer.parseInt(fields[index.numEnzyTermi]);
        proteinStart = Integer.parseInt(fields[index.protsIndex]);
        proteinEnd = Integer.parseInt(fields[index.proteIndex]);
        ms1Intensity = Double.parseDouble(fields[index.ms1IntIndex]);
        assignedModifications = fields[index.assignedModcIndex];
        ptmLocalization = index.ptmLocalcIndex == -1 ? "" : fields[index.ptmLocalcIndex];
        purity = Double.parseDouble(fields[index.purityIndex]);
        isUnique = Boolean.parseBoolean(fields[index.isUniquecIndex]);
        protein = fields[index.proteincIndex];
        proteinId = fields[index.proteinIDcIndex];
        gene = fields[index.genecIndex];
        mappedGenes = index.mapGeneIndex == -1 ? "" : fields[index.mapGeneIndex];

        modifiedPeptide = generateModifiedSequence(peptide, assignedModifications);
    }

    private boolean passResolutionSnr(String[] fields, List<Integer> naChannels) {
        if (parameters.addIsobaricFilter) {
            float snrSum = 0;
            for (int i = 0; i < parameters.channelNum; ++i) {
                if (!naChannels.contains(index.subplexIndices.get(0).abnIndex + i)) {
                    snrSum += Float.parseFloat(fields[index.snrOffset + i]);
                }
            }

            if (snrSum < parameters.minSNR) {
                return false;
            } else {
                for (int i = 0; i < parameters.channelNum; ++i) {
                    if (!naChannels.contains(index.subplexIndices.get(0).abnIndex + i)
                            && Integer.parseInt(fields[index.resOffset + i]) < parameters.minResolution
                            && Float.parseFloat(fields[index.snrOffset + i]) >= 1) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    private void parseChannels(String[] fields, List<Double> tmtIntensities, List<Integer> naChannels) {
        double sum = 0;
        for (int p = 0; p < parameters.subplexes.size(); p++) {
            Index.SubplexIndex si = index.subplexIndices.get(p);
            int channelCount = parameters.subplexes.get(p).channelCount;
            List<Double> plexChannels = new ArrayList<>();
            for (int i = 0; i < channelCount; i++) {
                try {
                    if (naChannels.contains(si.abnIndex + i)) {
                        continue;
                    }
                    double intensity = Double.parseDouble(fields[si.abnIndex + i]);
                    plexChannels.add(intensity);
                    sum += intensity;
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException("Invalid TMT intensity value: " + fields[si.abnIndex + i]);
                }
            }
            subplexChannels.add(plexChannels);
        }
        tmtIntensities.add(sum);
        sumTmtIntensity = sum;
    }

    private void extractGlycanInfo(String[] fields) {
        try {
            glycanQValue = parameters.glycoQval >= 0 ? Utils.tryParseDouble(fields[index.glycoQvalIndex]) : -1;
        } catch (IndexOutOfBoundsException ex) {
            System.out.println("Glycan FDR control requested but no such column found. No Glycan FDR applied.");
            parameters.glycoQval = -1;
            glycanQValue = -1;
        }
        glycanComposition = index.glycoCompositionIndex >= 0 ? fields[index.glycoCompositionIndex] : "";
    }

    private boolean updateLabelFlag() {
        if (parameters.allow_unlabeled) {
            return true;
        }

        if (matchLabels(assignedModifications, parameters.labels, 0.1f)) {
            return parameters.allow_overlabel || (!assignedModifications.contains("S(229.") && !assignedModifications.contains("S(304."));
        } else {
            return false;
        }
    }

    private boolean updateModificationFlag(Set<String> modTags) {
        if (parameters.modTagSet.contains("none")) {
            return true;
        }
        if (assignedModifications.isEmpty()) {
            return false;
        }
        for (String term : parameters.modTagSet) {
            if (isSkipModification(term)) {
                continue;
            }
            term = modifyTermForGlyco(term);
            if (assignedModifications.contains(term)) {
                return true;
            }

            if (term.equalsIgnoreCase("n-glyco") || term.equalsIgnoreCase("o-glyco")) {
                updateModAAParameter(term);
                if (processGlycoModification(term, modTags)) {
                    return true;
                }
            }
        }
        return false;
    }

    private boolean isSkipModification(String term) {
        if (term.equalsIgnoreCase("n-glyco") || term.equalsIgnoreCase("o-glyco")) {
            return parameters.glycoQval >= 0 && glycanQValue > parameters.glycoQval;
        }
        return false;
    }

    private String modifyTermForGlyco(String term) {
        if (!term.equalsIgnoreCase("N-glyco") && !term.equalsIgnoreCase("O-glyco")) {
            return term.substring(0, term.indexOf(".") + 3);
        }
        return term;
    }

    private void updateModAAParameter(String term) {
        if (term.equalsIgnoreCase("n-glyco")) {
            parameters.modAA = "N";
        } else if (term.equalsIgnoreCase("o-glyco")) {
            parameters.modAA = "S|T";
        }
    }

    private boolean processGlycoModification(String term, Set<String> modTags) {
        String[] modifications = assignedModifications.split(",");
        boolean foundGlycan = false;
        for (String modification : modifications) {
            if ((term.equalsIgnoreCase("N-glyco") && modification.contains("N(")) ||
                    (term.equalsIgnoreCase("O-glyco") &&
                            (modification.contains("S(") || modification.contains("T(")))) {
                double mass = extractMass(modification);
                if (mass >= 100) {
                    String modIndex = Utils.getAssignedModIndex(modification, glycanComposition, parameters.useGlycoComposition);
                    modTags.add(modIndex);
                    foundGlycan = true;
                }
            }
        }
        return foundGlycan;
    }

    private double extractMass(String modification) {
        Matcher matcher = Constants.GLYCO_MOD_PATTERN.matcher(modification);
        try {
            if (matcher.find()) {
                return Double.parseDouble(matcher.group(Constants.MASS_GROUP));
            }
            return Double.NaN;
        } catch (NumberFormatException e) {
            return Double.NaN;
        }
    }

    private boolean checkProteinExclusion() {
        if (protein.contains("contam_")) {
            return false;
        }
        if (parameters.protExcludeAry != null && !parameters.protExcludeAry[0].contains("none")) {
            for (String term : parameters.protExcludeAry) {
                if (protein.contains(term)) {
                    return false;
                }
            }
        }
        return true;
    }

    private int determineGeneCategory() {
        if (mappedGenes == null || mappedGenes.trim().isEmpty()) {
            return 0;
        }

        String[] genes = mappedGenes.trim().split(",");
        if (genes.length <= 1) {
            return handleSingleGene(genes[0]);
        } else {
            return handleMultipleGenes(genes);
        }
    }

    private int handleSingleGene(String mappedGene) {
        mappedGene = mappedGene.trim();
        if (mappedGene.isEmpty() || mappedGene.equalsIgnoreCase(gene.trim())) {
            return 2;
        } else {
            return parameters.allGeneSet.contains(mappedGene) ? 0 : 1;
        }
    }

    private int handleMultipleGenes(String[] mappedGenes) {
        for (String mappedGene : mappedGenes) {
            mappedGene = mappedGene.trim();
            if (!mappedGene.equalsIgnoreCase(gene.trim()) && parameters.allGeneSet.contains(mappedGene)) {
                return 0;
            }
        }
        return 1;
    }

    private void computeVirtualReference(ReferenceType refType) {
        switch (refType) {
            case SUMMATION:
                for (List<Double> ch : subplexChannels) {
                    subplexRefIntensities.add(ch.stream().mapToDouble(Double::doubleValue).sum());
                }
                break;
            case AVERAGE:
                for (List<Double> ch : subplexChannels) {
                    subplexRefIntensities.add(ch.stream().filter(v -> v > 0).mapToDouble(Double::doubleValue).average().orElse(0));
                }
                break;
            case MEDIAN:
                for (List<Double> ch : subplexChannels) {
                    List<Double> nonZero = ch.stream().filter(v -> v > 0).collect(Collectors.toList());
                    subplexRefIntensities.add(nonZero.isEmpty() ? 0 : Utils.takeMedian(nonZero));
                }
                break;
            default:
                for (int i = 0; i < subplexChannels.size(); i++) {
                    subplexRefIntensities.add(0.0);
                }
                break;
        }
    }

    private void processSitesWithProbT() {
        countModifiedSites();
        populateProbMap();
        extractSitePositions();
        phosphoSiteData.startIndex = Math.max(0, phosphoSiteData.startIndex + pepsIndex);
        phosphoSiteData.endIndex = Math.max(0, phosphoSiteData.endIndex + pepsIndex);
    }

    private void processSitesWithoutProbT() {
        List<Integer> positions = extractSitePositionsInMod();
        Collections.sort(positions);
        for (int position : positions) {
            phosphoSiteData.probMap.put(position, 0d);
            phosphoSiteData.siteLocalPos.append(peptide.charAt(position - 1)).append(pepsIndex + position);
        }
        phosphoSiteData.siteLocalCount = positions.size();
        phosphoSiteData.siteCount = positions.size();
        phosphoSiteData.startIndex = positions.get(0) + pepsIndex;
        phosphoSiteData.endIndex = positions.get(positions.size() - 1) + pepsIndex;
    }

    private void countModifiedSites() {
        String massTag = "(" + parameters.columntag.substring(parameters.columntag.indexOf(":") + 1);
        int index = assignedModifications.indexOf(massTag);
        while (index >= 0) {
            phosphoSiteData.siteCount++;
            index = assignedModifications.indexOf(massTag, index + 1);
        }
    }

    private void populateProbMap() {
        int leftIndex = ptmLocalization.indexOf("(");
        int rightIndex = ptmLocalization.indexOf(")");
        int gap = 0;
        boolean isFirst = true;

        while (leftIndex >= 0) {
            double probability = Double.parseDouble(ptmLocalization.substring(leftIndex + 1, rightIndex));
            int key;
            if (isFirst) {
                key = leftIndex;
                isFirst = false;
            } else {
                key = leftIndex - gap;
            }
            phosphoSiteData.probMap.put(key, probability);
            gap += (rightIndex - leftIndex + 1);

            leftIndex = ptmLocalization.indexOf("(", leftIndex + 1);
            rightIndex = ptmLocalization.indexOf(")", rightIndex + 1);
        }
    }

    private void extractSitePositions() {
        boolean notStart = false;
        for (Map.Entry<Integer, Double> entry : phosphoSiteData.probMap.entrySet()) {
            int key = entry.getKey();
            double probability = entry.getValue();
            if (probability > parameters.minSiteProb) {
                phosphoSiteData.siteLocalCount++;
                phosphoSiteData.siteLocalPos.append(peptide.charAt(key - 1)).append(pepsIndex + key);
            }
            if (!notStart) {
                phosphoSiteData.startIndex = key;
                notStart = true;
            }
            phosphoSiteData.endIndex = key;
        }
    }

    private List<Integer> extractSitePositionsInMod() {
        List<Integer> positions = new ArrayList<>();
        String[] assignedMods = assignedModifications.split(",");
        for (String assignedMod : assignedMods) {
            String mod = Utils.getAssignedModIndex(assignedMod, glycanComposition, parameters.useGlycoComposition);
            if (Utils.isModTagMatch(mod, parameters.modTagSet)) {
                int position = Integer.parseInt(assignedMod.substring(0, assignedMod.indexOf("(") - 1).trim());
                positions.add(position);
                updateSiteLocalMassMap(position, mod);
            }
        }
        return positions;
    }

    private void updateSiteLocalMassMap(int position, String mod) {
        mod = mod.replace("(", "").replace(")", "").substring(1);
        List<String> modMassList = phosphoSiteData.siteLocalMassMap.computeIfAbsent(position, k -> new ArrayList<>());
        modMassList.add(mod);
    }

    private void generateGroupKey(GroupBy groupBy) {
        switch (groupBy) {
            case GENE:
                groupKey = gene.isEmpty() ? proteinId : gene;
                break;
            case PROTEIN_ID:
                groupKey = proteinId;
                break;
            case PEPTIDE:
                generatePeptideGroupKey();
                break;
            case MULTI_PHOSPHO_SITE:
                handleMultiPhosphoSites();
                break;
            case SINGLE_PHOSPHO_SITE:
                break;
            case MULTI_MASS_GLYCO:
                handleMultiMassGlyco();
                break;
            case MODIFIED_PEPTIDE:
                generateModifiedPeptideGroupKey();
                break;
        }
    }

    private void generateModifiedPeptideGroupKey() {
        groupKey = proteinId + "%" + modifiedPeptide;
    }

    private void generatePeptideGroupKey() {
        if (parameters.minSiteProb >= 0 && phosphoSiteData.siteCount > 0) {
            groupKey = proteinId + "%" + phosphoSiteData.startIndex + "%" + phosphoSiteData.endIndex;
        } else if (parameters.minSiteProb < 0) {
            groupKey = proteinId + "%" + peptide;
        }
    }

    private void handleMultiPhosphoSites() {
        if (phosphoSiteData.siteCount <= 0) {
            return;
        }
        StringBuilder key = new StringBuilder(proteinId + "%" + phosphoSiteData.startIndex + "%" + phosphoSiteData.endIndex +
                "%" + phosphoSiteData.siteCount + "%" + phosphoSiteData.siteLocalCount);
        if (phosphoSiteData.siteLocalPos.length() > 0) {
            key.append("%").append(phosphoSiteData.siteLocalPos);
        }
        String[] peptideKeys = peptide.split("-");
        if (peptideKeys.length > 1) {
            key.append("%").append(peptideKeys[1]);
        }
        markLocalizedSites();
        groupKey = key.toString();
    }

    private void handleMultiMassGlyco() {
        if (phosphoSiteData.siteCount <= 0) {
            return;
        }
        StringBuilder siteLocalMass = new StringBuilder();
        for (List<String> massList : phosphoSiteData.siteLocalMassMap.values()) {
            siteLocalMass.append(String.join("_", massList)).append("_");
        }
        siteLocalMass.deleteCharAt(siteLocalMass.length() - 1);
        StringBuilder key = new StringBuilder(proteinId + "%" + phosphoSiteData.startIndex + "%" + phosphoSiteData.endIndex +
                "%" + phosphoSiteData.siteCount + "%" + phosphoSiteData.siteLocalCount);
        if (phosphoSiteData.siteLocalPos.length() > 0) {
            key.append("%").append(phosphoSiteData.siteLocalPos);
        }
        key.append("%").append(siteLocalMass);
        markLocalizedSites();
        groupKey = key.toString();
    }

    private void markLocalizedSites() {
        StringBuilder newPeptideSeq = new StringBuilder();
        for (int i = 0; i < peptide.length(); i++) {
            char aa = peptide.charAt(i);
            if (phosphoSiteData.probMap.containsKey(i + 1) && phosphoSiteData.probMap.get(i + 1) >= parameters.minSiteProb) {
                newPeptideSeq.append(Character.toLowerCase(aa));
            } else {
                newPeptideSeq.append(aa);
            }
        }
        peptide = newPeptideSeq.toString();
    }

    private String generateModifiedSequence(String sequence, String assignMods) {
        if (assignMods.isEmpty()) return sequence;

        int seqLen = sequence.length();
        String[] modifications = assignMods.split(",");
        Map<Integer, String> modMap = new LinkedHashMap<>();

        for (String mod : modifications) {
            mod = mod.trim();
            int openIdx = mod.indexOf('(');
            int closeIdx = mod.indexOf(')');

            if (openIdx == -1 || closeIdx == -1 || closeIdx < openIdx) {
                System.err.println("Warning: Could not parse assigned modification: " + mod + " for spectrum " + spectrum);
                System.exit(1);
            }

            String site = mod.substring(0, openIdx);
            String mass = "[" + mod.substring(openIdx + 1, closeIdx) + "]";

            if (site.equalsIgnoreCase("N-term")) {
                modMap.put(0, mass);
            } else if (site.equalsIgnoreCase("C-term")) {
                modMap.put(seqLen + 1, mass);
            } else {
                int aaIdx = -1;
                for (int i = 0; i < site.length(); i++) {
                    if (Character.isLetter(site.charAt(i))) {
                        aaIdx = i;
                        break;
                    }
                }
                if (aaIdx != -1) {
                    int index = Integer.parseInt(site.substring(0, aaIdx));
                    modMap.put(index, mass);
                }
            }
        }

        StringBuilder modifiedSequence = new StringBuilder();
        int sequenceLength = sequence.length();

        if (modMap.containsKey(0)) {
            modifiedSequence.append("n").append(modMap.get(0));
        }

        for (int i = 1; i <= sequenceLength; i++) {
            modifiedSequence.append(sequence.charAt(i - 1));
            if (modMap.containsKey(i)) {
                modifiedSequence.append(modMap.get(i));
            }
        }

        if (modMap.containsKey(seqLen + 1)) {
            modifiedSequence.append("c").append(modMap.get(0));
        }

        return modifiedSequence.toString();
    }
    // endregion
}
