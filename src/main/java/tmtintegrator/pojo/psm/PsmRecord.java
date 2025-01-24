package tmtintegrator.pojo.psm;

import tmtintegrator.constants.Constants;
import tmtintegrator.constants.GroupBy;
import tmtintegrator.constants.ReferenceType;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.utils.Utils;

import java.util.*;
import java.util.regex.Matcher;
import java.util.stream.Collectors;

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
    private double observedMz;
    private double calculatedPepMass;
    private double probability;
    private int numEnzymaticTermini;
    private int proteinStart;
    private int proteinEnd;
    private double ms1Intensity;
    private String assignedModifications;
    private String observedModifications;
    private String ptmLocalization;
    private double purity;
    private boolean isUnique;
    private String protein;
    private String proteinId;
    private String entryName;
    private String gene;
    private String proteinDescription;
    private String mappedGenes;
    private String mappedProteins;
    private List<Double> channels;
    // endregion

    // region extra fields
    private double sumTmtIntensity;
    private double refIntensity;
    private boolean isVirtualReference;
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
    private boolean allowOverLabel;
    private boolean labelFlag;
    private boolean modficationFlag;
    private boolean isPeptideRetained;
    private boolean proteinExclusionFlag;
    private int geneCategory;
    // endregion

    // region backup fields
    private String copyPeptide;
    private double copyRefIntensity;
    private List<Double> copyChannels;
    // endregion


    public PsmRecord(Parameters parameters, Index index) {
        this.parameters = parameters;
        this.index = index;
        channels = new ArrayList<>();
    }

    // region getters and setters
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

    public double getCopyRefIntensity() {
        return copyRefIntensity;
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
     * Parse a PSM record from a line in a PSM.tsv file
     *
     * @param line         line in the PSM.tsv file
     * @param tmtIntensities list of TMT intensities
     * @param naChannels list of NA channels indices
     */
    public void parsePsmRecord(String line, List<Double> tmtIntensities, List<Integer> naChannels) {
        String[] fields = line.split("\t");
        parseRegularFields(fields);
        pepsIndex = proteinStart - 1;
        parseChannels(fields, tmtIntensities, naChannels);
        updateRefIntensity(fields);
        extractGlycanInfo(fields);
    }

    public void backup() {
        // backup the original values for fields that will be modified
        copyRefIntensity = refIntensity;
        copyPeptide = peptide;
        copyChannels = new ArrayList<>(channels);
    }

    public void reset() {
        // reset the modified fields to the original values
        isExcluded = false;
        refIntensity = copyRefIntensity;
        peptide = copyPeptide;
        channels = new ArrayList<>(copyChannels);
    }

    /**
     * Update flags for filtering
     *
     * @param modTags set of mod tags
     */
    public void updateFlags(Set<String> modTags) {
        allowOverLabel = parameters.allow_overlabel || (!assignedModifications.contains("S(229.")); // FIXME
        labelFlag = updateLabelFlag();
        modficationFlag = updateModificationFlag(modTags);
        isPeptideRetained = !parameters.uniquePep || isUnique;
        proteinExclusionFlag = checkProteinExclusion();
        geneCategory = determineGeneCategory();
    }

    /**
     * Check if the PSM entry passes the criteria
     *
     * @param tmtThreshold TMT threshold
     * @return true if the PSM entry passes the criteria
     */
    public boolean isPassCriteria(double tmtThreshold) {
        return purity >= parameters.minPurity &&
                probability >= parameters.minPepProb &&
                sumTmtIntensity >= tmtThreshold &&
                allowOverLabel &&
                labelFlag &&
                modficationFlag &&
                isPeptideRetained &&
                proteinExclusionFlag &&
                geneCategory >= parameters.uniqueGene &&
                numEnzymaticTermini >= parameters.min_ntt &&
                refIntensity > 0;
    }

    public String generatePsmKey() {
        // extract the filename from Spectrum
        String fileName = spectrum.substring(spectrum.lastIndexOf("_") + 1, spectrum.indexOf("."));
        return fileName + "_" + peptide + "_" + charge + "_" + calculatedPepMass;
    }

    /**
     * Analyze phospho sites in the PSM entry
     *
     * @param groupBy group by option
     */
    public void analyzePhosphoSites(GroupBy groupBy) {
        // FIXME: analyze phospho sites only need once
        // initialize phospho site data and group key
        phosphoSiteData = new PhosphoSiteData();
        groupKey = "";

        // get num of phospho sites with probability > threshold and their positions
        phosphoSiteData.startIndex = pepsIndex;
        phosphoSiteData.endIndex = proteinEnd - 1;

        if (parameters.minSiteProb > 0) {
            analyzeWithThreshold();
        } else if (parameters.minSiteProb == 0) {
            analyzeWithoutThreshold();
        }

        // use MS1 intensity as reference intensity
        if (parameters.ms1Int) {
            refIntensity = useMs1Intensity();
        }

        // generate group key
        generateGroupKey(groupBy);

        // filter out psm without phospho sites
        if (groupBy != GroupBy.SINGLE_PHOSPHO_SITE && groupKey.isEmpty()) {
            isExcluded = true;
        }
    }

    // region helper methods
    private void updateRefIntensity(String[] fields) {
        ReferenceType refType = ReferenceType.fromValue(parameters.add_Ref);
        if (refType == ReferenceType.NONE || refType == ReferenceType.RAW_ABUNDANCE) {
            refIntensity = Double.parseDouble(fields[index.refIndex]);
        } else {
            isVirtualReference = true;
            parameters.refTag = "Virtual_Reference";
            refIntensity = computeVirtualReference(refType);
        }
    }

    private void parseRegularFields(String[] fields) {
        spectrum = fields[index.spectrumIndex];
        peptide = fields[index.peptideIndex];
        modifiedPeptide = fields[index.modifiedPeptideIndex];
        extendedPeptide = fields[index.extpepIndex];
        charge = Integer.parseInt(fields[index.chargeIndex]);
        retention = Double.parseDouble(fields[index.rtIndex]);
        observedMz = Double.parseDouble(fields[index.observedMzIndex]);
        calculatedPepMass = Double.parseDouble(fields[index.pepMassIndex]);
        probability = Double.parseDouble(fields[index.pepProbcIndex]);
        numEnzymaticTermini = Integer.parseInt(fields[index.numEnzyTermi]);
        proteinStart = Integer.parseInt(fields[index.protsIndex]);
        proteinEnd = Integer.parseInt(fields[index.proteIndex]);
        ms1Intensity = Double.parseDouble(fields[index.ms1IntIndex]);
        assignedModifications = fields[index.assignedModcIndex];
        observedModifications = fields[index.observedModIndex];
        ptmLocalization = index.ptmLocalcIndex == -1 ? "" : fields[index.ptmLocalcIndex];
        purity = Double.parseDouble(fields[index.purityIndex]);
        isUnique = Boolean.parseBoolean(fields[index.isUniquecIndex]);
        protein = fields[index.proteincIndex];
        proteinId = fields[index.proteinIDcIndex];
        entryName = fields[index.entryNameIndex];
        gene = fields[index.genecIndex];
        proteinDescription = fields[index.proteinDescIndex];
        mappedGenes = index.mapGeneIndex == -1 ? "" : fields[index.mapGeneIndex];
        mappedProteins = index.mappedProteinsIndex == -1 ? "" : fields[index.mappedProteinsIndex];
    }

    private void parseChannels(String[] fields, List<Double> tmtIntensities, List<Integer> naChannels) {
        double sum = 0;
        for (int i = index.abnIndex; i < fields.length; i++) {
            try {
                if (naChannels.contains(i)) {
                    continue;
                }
                double intensity = Double.parseDouble(fields[i]);
                channels.add(intensity);
                sum += intensity;
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Invalid TMT intensity value: " + fields[i]);
            }
        }
        tmtIntensities.add(sum);
        sumTmtIntensity = sum;
    }

    private void extractGlycanInfo(String[] fields) {
        try {
            glycanQValue = parameters.glycoQval >= 0 ? Utils.tryParseDouble(fields[index.glycoQvalIndex]) : -1;
        } catch (IndexOutOfBoundsException ex) {
            // non-glyco search
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
        return assignedModifications.contains("n(42.") ||
                assignedModifications.contains("n(229.") ||
                assignedModifications.contains(", 1S(229.") ||
                assignedModifications.indexOf("1S(229.") == 0 ||
                assignedModifications.contains("N-term(42.") ||
                assignedModifications.contains("N-term(229.") ||
                assignedModifications.contains("N-term(144.1") ||
                assignedModifications.contains("9K(304.2") ||
                assignedModifications.contains("N-term(304.2");
    }

    private boolean updateModificationFlag(Set<String> modTags) {
        if (parameters.modTagSet.contains("none")) {
            // all PSMs are valid in terms of modification
            return true;
        }
        if (assignedModifications.isEmpty()) {
            // skip PSMs with nothing in assigned mods column if mods requested
            // (shouldn't happen, but can if two mods are reported on same site and Philosopher can't separate them)
            return false;
        }
        for (String term : parameters.modTagSet) {
            // skip this mod tag if psm with glycan q-value above threshold
            if (isSkipModification(term)) {
                continue;
            }
            term = modifyTermForGlyco(term);
            if (assignedModifications.contains(term)) {
                return true;
            }

            // special handling for glycan
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
            return term.substring(0, term.indexOf(".") + 3); // TODO: use floor to cut float
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

    /**
     * Process glycan modification
     *
     * @param term    N-glyco or O-glyco
     * @param modTags list of new mod tags
     * @return true if glycan modification is processed
     */
    private boolean processGlycoModification(String term, Set<String> modTags) {
        String[] modifications = assignedModifications.split(",");
        for (String modification : modifications) {
            if ((term.equalsIgnoreCase("N-glyco") && modification.contains("N(")) ||
                    (term.equalsIgnoreCase("O-glyco") &&
                            (modification.contains("S(") || modification.contains("T(")))) {
                double mass = extractMass(modification); // TODO: this only works for modification with one "(mass)"
                if (mass >= 100) {
                    String modIndex = Utils.getAssignedModIndex(modification, glycanComposition, parameters.useGlycoComposition);
                    modTags.add(modIndex);
                    return true;
                }
            }
        }
        return false;
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
        // exclude contaminants
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
                // found a gene in the list that is not the target gene
                return 0;
            }
        }
        return 1;
    }

    private double computeVirtualReference(ReferenceType refType) {
        switch (refType) {
            case SUMMATION:
                return summationReference();
            case AVERAGE:
                return averageReference();
            case MEDIAN:
                return medianReference();
        }
        return 0;
    }

    private double summationReference() {
        return sumTmtIntensity;
    }

    private double averageReference() {
        // take average of all non zero values
        return channels.stream().filter(value -> value > 0).mapToDouble(Double::doubleValue).average().orElse(0);
    }

    private double medianReference() {
        // take median of all non zero values
        List<Double> nonZeroValues = channels.stream().filter(value -> value > 0).collect(Collectors.toList());
        return nonZeroValues.isEmpty() ? 0 : Utils.takeMedian(nonZeroValues);
    }

    private void analyzeWithThreshold() {
        countModifiedSites();
        populateProbMap();
        extractSitePositions();
        // update start and end index
        phosphoSiteData.startIndex = Math.max(0, phosphoSiteData.startIndex + pepsIndex);
        phosphoSiteData.endIndex = Math.max(0, phosphoSiteData.endIndex + pepsIndex);
    }

    private void analyzeWithoutThreshold() {
        List<Integer> positions = extractSitePositionsInMod();
        Collections.sort(positions);
        for (int position : positions) {
            phosphoSiteData.probMap.put(position, 0d);
            // extracted position is 1-based, need to convert to 0-based for peptide
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
            // extract probability
            double probability = Double.parseDouble(ptmLocalization.substring(leftIndex + 1, rightIndex));
            int key;
            if (isFirst) {
                key = leftIndex; // 1-based AA local position
                isFirst = false;
            } else {
                key = leftIndex - gap;
            }
            phosphoSiteData.probMap.put(key, probability);
            // record gap to locate AA in peptide sequence without probability
            gap += (rightIndex - leftIndex + 1); // size of a probability string "(0.1234)" including parentheses

            // update indices
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
                // extract position number, ex: 11N(2569.9045) -> 11
                int position = Integer.parseInt(assignedMod.substring(0, assignedMod.indexOf("(") - 1).trim());
                positions.add(position);
                updateSiteLocalMassMap(position, mod); // for multi-mass glycan
            }
        }
        return positions;
    }

    private void updateSiteLocalMassMap(int position, String mod) {
        // simplify the mod string for index: remove parentheses and first AA code
        mod = mod.replace("(", "").replace(")", "").substring(1);
        List<String> modMassList = phosphoSiteData.siteLocalMassMap.computeIfAbsent(position, k -> new ArrayList<>());
        modMassList.add(mod);
    }

    private double useMs1Intensity() {
        if (isVirtualReference) {
            return ms1Intensity * (refIntensity / (sumTmtIntensity + refIntensity));
        }
        return ms1Intensity * (refIntensity / sumTmtIntensity);
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
                //  single phospho site processing will be handled in the second round
                break;
            case MULTI_MASS_GLYCO:
                handleMultiMassGlyco();
                break;
        }
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
        // generate group key
        StringBuilder key = new StringBuilder(proteinId + "%" + phosphoSiteData.startIndex + "%" + phosphoSiteData.endIndex +
                "%" + phosphoSiteData.siteCount + "%" + phosphoSiteData.siteLocalCount);
        if (phosphoSiteData.siteLocalPos.length() > 0) {
            key.append("%").append(phosphoSiteData.siteLocalPos);
        }
        String[] peptideKeys = peptide.split("-"); // FIXME: what's the format
        if (peptideKeys.length > 1) {
            key.append("%").append(peptideKeys[1]);
        }
        // mark localized sites in peptide sequence
        markLocalizedSites();

        groupKey = key.toString();
    }

    private void handleMultiMassGlyco() {
        if (phosphoSiteData.siteCount <= 0) {
            return;
        }
        // generate site localized mass string
        StringBuilder siteLocalMass = new StringBuilder();
        for (List<String> massList : phosphoSiteData.siteLocalMassMap.values()) {
            siteLocalMass.append(String.join("_", massList)).append("_");
        }
        siteLocalMass.deleteCharAt(siteLocalMass.length() - 1); // remove trailing '_'
        // generate group key
        StringBuilder key = new StringBuilder(proteinId + "%" + phosphoSiteData.startIndex + "%" + phosphoSiteData.endIndex +
                "%" + phosphoSiteData.siteCount + "%" + phosphoSiteData.siteLocalCount);
        if (phosphoSiteData.siteLocalPos.length() > 0) {
            key.append("%").append(phosphoSiteData.siteLocalPos);
        }
        key.append("%").append(siteLocalMass);
        // mark localized sites in peptide sequence
        markLocalizedSites();

        groupKey = key.toString();
    }

    /**
     * Mark localized sites in peptide sequence.
     */
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
    // endregion
}
