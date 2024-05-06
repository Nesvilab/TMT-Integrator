package tmtintegrator.pojo;

import tmtintegrator.utils.Utils;

import java.util.List;

/**
 * PSM entry class representing a single line in the PSM file
 *
 * @author rogerli on 04/2024
 */
public class PsmEntry {

    // region PSM file information
    private final ds_Parameters parameters;
    private final ds_Index index;
    private final String[] fields;
    // endregion

    // region Psm Entry Fields
    private String filename;
    private String peptideSequence;
    private double purity;
    private double tmtIntensity;
    private double peptideProbability;
    private boolean isUnique;
    private String assignedModification;
    private String proteinName;
    private String proteinId;
    private String gene;
    private String mappedGenes;
    private double referenceIntensity;
    private int numEnzymaticTermini;
    private double glycanQValue;
    private String glycanComposition;
    // endregion

    // region Configuration Fields
    private boolean allowOverLabel;
    private boolean labelFlag;
    private boolean modficationFlag;
    private boolean uniquePeptideFlag;
    private boolean proteinExclusionFlag;
    // endregion

    private int geneCategory;

    public PsmEntry(String psm, ds_Parameters parameters, ds_Index index) {
        this.parameters = parameters;
        this.index = index;
        this.fields = psm.split("\t");
    }

    // region Getters
    public String getFilename() {
        return filename;
    }

    public String getPeptideSequence() {
        return peptideSequence;
    }

    public double getPurity() {
        return purity;
    }

    public double getTmtIntensity() {
        return tmtIntensity;
    }

    public double getPeptideProbability() {
        return peptideProbability;
    }

    public boolean isUnique() {
        return isUnique;
    }

    public String getAssignedModification() {
        return assignedModification;
    }

    public String getProteinName() {
        return proteinName;
    }

    public String getProteinId() {
        return proteinId;
    }

    public String getGene() {
        return gene;
    }

    public String getMappedGenes() {
        return mappedGenes;
    }

    public double getReferenceIntensity() {
        return referenceIntensity;
    }

    public int getNumEnzymaticTermini() {
        return numEnzymaticTermini;
    }

    public double getGlycanQValue() {
        return glycanQValue;
    }

    public String getGlycanComposition() {
        return glycanComposition;
    }

    public boolean isAllowOverLabel() {
        return allowOverLabel;
    }

    public boolean isLabelFlag() {
        return labelFlag;
    }

    public boolean isModficationFlag() {
        return modficationFlag;
    }

    public boolean isUniquePeptideFlag() {
        return uniquePeptideFlag;
    }

    public boolean isProteinExclusionFlag() {
        return proteinExclusionFlag;
    }

    public int getGeneCategory() {
        return geneCategory;
    }

    // endregion

    /**
     * Parse the PSM entry and extract the relevant information
     */
    public void parsePsmEntry() {
        extractFileName();
        peptideSequence = fields[index.pepcIndex];
        purity = Double.parseDouble(fields[index.purityIndex]);
        tmtIntensity = Double.parseDouble(fields[fields.length - 1]);
        peptideProbability = Double.parseDouble(fields[index.pepProbcIndex]);
        isUnique = Boolean.parseBoolean(fields[index.isUniquecIndex]);
        assignedModification = fields[index.assignedModcIndex];
        proteinName = fields[index.proteincIndex];
        proteinId = fields[index.proteinIDcIndex];
        gene = !fields[index.genecIndex].isEmpty() ? fields[index.genecIndex] : proteinId;
        mappedGenes = index.mapGeneIndex >= 0 ? fields[index.mapGeneIndex] : "";
        referenceIntensity = parameters.add_Ref < 0 ? Double.parseDouble(fields[index.refIndex]) : 10000; // TODO: set random value to pass the criteria
        numEnzymaticTermini = Integer.parseInt(fields[index.numEnzyTermi]);
        extractGlycoQValue();
        glycanComposition = index.glycoCompositionIndex == -1 ? "" : fields[index.glycoCompositionIndex]; // required for glycan processing
    }

    /**
     * Check the configurations for the PSM entry and update the flags <br>
     * Update modTagList if new mod tags are found
     *
     * @param newModTagList list of new mod tags
     */
    public void checkConfigurations(List<String> newModTagList) {
        allowOverLabel = parameters.allow_overlabel || (!assignedModification.contains("S(229."));
        labelFlag = updateLabelFlag();
        modficationFlag = checkModifications(newModTagList);
        uniquePeptideFlag = !parameters.uniquePep || isUnique; // TODO: seems like a bug
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
                peptideProbability >= parameters.minPepProb &&
                tmtIntensity >= tmtThreshold &&
                allowOverLabel &&
                labelFlag &&
                modficationFlag &&
                uniquePeptideFlag &&
                proteinExclusionFlag &&
                geneCategory >= parameters.uniqueGene &&
                numEnzymaticTermini >= parameters.min_ntt &&
                referenceIntensity > 0;
    }

    /**
     * Get the processed PSM entry
     *
     * @return processed PSM entry string line
     */
    public String getProcessedPsm() {
        fields[index.isUsedIndex] = "true";
        StringBuilder newPsm = new StringBuilder();
        for (String field : fields) {
            newPsm.append(field).append("\t");
        }
        newPsm.append(geneCategory);
        return newPsm.toString();
    }

    /**
     * Generate the PSM key
     *
     * @return PSM key
     */
    public String generatePsmKey() {
        return filename + "_" + fields[index.peptideIndex] + "_" +
                fields[index.chargeIndex] + "_" + fields[index.pepMassIndex];
    }

    // region Helper Methods============================================================================================
    private void extractFileName() {
        filename = fields[0].substring(fields[0].lastIndexOf('_') + 1, fields[0].indexOf('.'));
    }

    private void extractGlycoQValue() {
        try {
            glycanQValue = parameters.glycoQval >= 0 ? Utils.tryParseDouble(fields[index.glycoQvalIndex]) : -1;
        } catch (IndexOutOfBoundsException ex) {
            // non-glyco search
            System.out.println("Glycan FDR control requested but no such column found. No Glycan FDR applied.");
            parameters.glycoQval = -1;
            glycanQValue = -1;
        }
    }

    private boolean updateLabelFlag() {
        if (parameters.allow_unlabeled) {
            return true;
        }
        String firstAAtag = peptideSequence.charAt(0) + "(229."; // TODO: No usage, don't know what for
        return assignedModification.contains("1S(229.") ||
                assignedModification.contains("N-term(42.") ||
                assignedModification.contains("N-term(229.") ||
                assignedModification.contains("N-term(144.1") ||
                assignedModification.contains("9K(304.2") ||
                assignedModification.contains("N-term(304.2");
        // TODO: seems incorrect below
        // assignedModification.contains("n(42.") ||
        // assignedModification.contains("n(229.");
    }

    /**
     * Check if the assigned modification meets the criteria
     *
     * @param newModTagList list of new mod tags
     * @return true if the modification is valid
     */
    private boolean checkModifications(List<String> newModTagList) {
        if (parameters.modTagLi.get(0).trim().equalsIgnoreCase("none")) {
            // If 'none' is the specified modification, all PSMs are valid in terms of modification.
            return true;
        }
        if (assignedModification.isEmpty()) {
            // TODO: may need skip current PSM
            // skip PSMs with nothing in assigned mods column if mods requested
            // (shouldn't happen, but can if two mods are reported on same site and Philosopher can't separate them)
            return false;
        }
        for (String term : parameters.modTagLi) {
            // skip if psm glycan q-value is above threshold
            if (isSkipModification(term)) {
                continue;
            }
            term = modifyTermForGlyco(term);
            if (assignedModification.contains(term)) {
                return true;
            }

            // Special handling for glycans
            if (term.equalsIgnoreCase("N-glyco") || term.equalsIgnoreCase("O-glyco")) {
                updateModAAParameter(term);
                if (processGlycoModification(term, newModTagList)) {
                    return true;
                }
            }
        }
        return false;
    }

    private boolean isSkipModification(String term) {
        if (term.equalsIgnoreCase("N-glyco") || term.equalsIgnoreCase("O-glyco")) {
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
        if (term.equalsIgnoreCase("N-glyco")) {
            parameters.modAA = "N";
        } else if (term.equalsIgnoreCase("O-glyco")) {
            parameters.modAA = "S|T";
        }
    }

    /**
     * Process glycan modification
     *
     * @param term          N-glyco or O-glyco
     * @param newModTagList list of new mod tags
     * @return true if glycan modification is processed
     */
    private boolean processGlycoModification(String term, List<String> newModTagList) {
        String[] modifications = assignedModification.split(",");
        for (String modification : modifications) {
            if ((term.equalsIgnoreCase("N-glyco") && modification.contains("N(")) ||
                    (term.equalsIgnoreCase("O-glyco") &&
                            (modification.contains("S(") || modification.contains("T(")))) {
                double mass = extractMass(modification); // TODO: this only works for modification with one "(mass)"
                if (mass >= 100) {
                    String modIndex = getAssignedModIndex(modification, glycanComposition, parameters.useGlycoComposition);
                    if (!newModTagList.contains(modIndex)) {
                        newModTagList.add(modIndex);
                    }
                    return true;
                }
            }
        }
        return false;
    }

    private double extractMass(String modification) {
        try {
            int start = modification.indexOf("(") + 1;
            int end = modification.indexOf(")");
            return Double.parseDouble(modification.substring(start, end));
        } catch (NumberFormatException e) {
            return Double.NaN; // TODO: need test
        }
    }

    /**
     * Helper method for getting the part of a modification to use as the index. Supports using the glycan composition
     * instead of mass if specified.
     * Note: if observed mods is empty, default to using the mass value instead
     *
     * @param inputMod             single Assigned modification string
     * @param glycanComposition    contents of the corresponding glycan composition column (only needed for glyco mode)
     * @param useGlycanComposition whether to use the glycan composition or mass as the index
     * @return mod string to use as index
     */
    public String getAssignedModIndex(String inputMod, String glycanComposition, boolean useGlycanComposition) {
        String mod;
        if (useGlycanComposition && !glycanComposition.isEmpty()) {
            // if using composition for index, read it from glycan composition column. Still get AA site from assigned mods
            mod = inputMod.substring(inputMod.indexOf("(") - 1, inputMod.indexOf("("));
            mod = String.format("%s(%s)", mod, glycanComposition);
        } else {
            // read mass from assigned mod, as would do for any other mod
            mod = inputMod.substring(inputMod.indexOf("(") - 1);
        }
        return mod;
    }

    private boolean checkProteinExclusion() {
        // exclude contaminants
        if (proteinName.contains("contam")) {
            return false;
        }
        if (parameters.protExcludeAry != null && !parameters.protExcludeAry[0].contains("none")) {
            for (String term : parameters.protExcludeAry) {
                if (proteinName.contains(term)) {
                    return false;
                }
            }
        }
        return true;
    }

    private int determineGeneCategory() {
        if (mappedGenes == null || mappedGenes.trim().isEmpty()) {
            // Assuming 0 for no mapping genes
            return 0; // TODO: magic number, use parameter
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
            return 2; // TODO: magic number, use parameter
        } else {
            return parameters.AllGeneLi.contains(mappedGene) ? 0 : 1; // TODO: magic number, use parameter
        }
    }

    private int handleMultipleGenes(String[] mappedGenes) {
        for (String mappedGene : mappedGenes) {
            mappedGene = mappedGene.trim();
            if (!mappedGene.equalsIgnoreCase(gene.trim()) && parameters.AllGeneLi.contains(mappedGene)) {
                // found a gene in the list that is not the target gene
                return 0; // TODO: magic number, use parameter
            }
        }
        return 1; // TODO: magic number, use parameter
    }
    // endregion========================================================================================================
}
