package tmtintegrator.integrator;

import tmtintegrator.pojo.ProteinIndex;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.psm.Psm;
import tmtintegrator.pojo.psm.PsmRecord;
import tmtintegrator.utils.ReportData;

import java.io.*;
import java.util.*;

/**
 * Preprocess PSM files, check for missing values, create index for each PSM file
 *
 * @author rogerli on 05/2024
 */
public class PsmPreProcessor {

    private final Parameters parameters;
    private final ReportData reportData;

    public PsmPreProcessor(Parameters parameters, ReportData reportData) {
        this.parameters = parameters;
        this.reportData = reportData;
    }

    /**
     * Check PSM tables for missing values <br>
     * Update protein list <br>
     * Get all genes <br>
     * Create index for each PSM file <br>
     * Maintain extra info of each PSM for final report
     *
     * @throws IOException if an I/O error occurs
     */
    public void checkPsmAndBuildIndex() throws IOException {
        Collections.sort(parameters.fNameLi);

        for (File psmFile : parameters.fileList) {
            checkPsmFile(psmFile);
        }
    }

    /**
     * Collect protein information from protein.tsv files <br>
     * - Protein ID <br>
     * - Organism <br>
     * - Indistinguishable Proteins <br>
     */
    public void collectProteinInfo() throws IOException {
        for (File proteinFile : parameters.proteinFileList) {
            try (BufferedReader reader = new BufferedReader(new FileReader(proteinFile))) {
                String titleLine = reader.readLine();
                ProteinIndex proteinIndex = parameters.proteinIndexMap.get(proteinFile.getAbsolutePath());
                mapColumnIndex(titleLine, proteinIndex);

                String line;
                while ((line = reader.readLine()) != null) {
                    String[] fields = line.split("\t");
                    reportData.updateExtraProteinInfo(fields, proteinIndex);
                }
            }
        }
    }

    /**
     * Preprocess PSM files <br>
     * - Read PSM file and preprocess each line <br>
     * - Calculate TMT threshold <br>
     * - Collapse PSM lines based on the criteria <br>
     * - Select best PSM if required <br>
     *
     * @return list of PSM objects
     */
    public List<Psm> parseAndFilterPsmFiles() {
        List<Psm> psmList = new ArrayList<>();
        for (File psmFile : parameters.fileList) {
            try {
                Psm psm = new Psm(parameters, psmFile);
                // parse PSM file and calculate tmt threshold
                List<Double> tmtIntensities = new ArrayList<>();
                psm.readPsmFile(psmFile, tmtIntensities);
                double tmtThreshold = calculateTmtThreshold(psm.getPsmRecords().size(), tmtIntensities);

                // collapse PSM lines based on the criteria
                Map<String, List<PsmRecord>> psmMap = new HashMap<>();
                psmMap.put("NotUsed", new ArrayList<>());
                collapsePsm(psm, psmMap, tmtThreshold);

                // select best PSM if required
                selectBestPsm(psmMap);

                // keep only used PSMs
                psm.filterUnUsedPsm(psmMap);
                psmList.add(psm);
            } catch (IOException e) {
                System.err.println("Error processing PSM file: " + psmFile.getAbsolutePath());
                e.printStackTrace();
            }
        }
        return psmList;
    }

    // region checkPsmAndBuildIndex helper methods=============================================
    private void checkPsmFile(File psmFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(psmFile))) {
            // process header and column index
            String title = reader.readLine();

            // checking existence of the resolution and SNR columns
            if (title.contains("Resolution ") && title.contains("SNR ")) {
                parameters.addIsobaricFilter = true;
            } else if (title.contains("Resolution ") || title.contains("SNR ")) {
                System.err.println("Error: Both Resolution and SNR columns should be present in the PSM file.");
                System.exit(1);
            }

            getColumnIndex(title, psmFile);
            Index index = parameters.indMap.get(psmFile.getAbsolutePath());

            // process data
            int totalCount = 0;
            int ms1MissingCount = 0;
            int geneMissingCount = 0;

            String line;
            while ((line = reader.readLine()) != null) {
                totalCount++;
                String[] fields = line.split("\t");

                if (isFieldMissing(fields, index, "ms1")) {
                    ms1MissingCount++;
                }
                if (isFieldMissing(fields, index, "gene")) {
                    geneMissingCount++;
                } else {
                    updateAllGenesList(fields, index);
                }
                reportData.updateExtraPsmInfo(fields, index);
            }

            checkAndWarn(totalCount, ms1MissingCount, psmFile, "ms1");
            checkAndWarn(totalCount, geneMissingCount, psmFile, "gene");
        }
    }

    /**
     * check if the field is missing by value
     *
     * @param fields array of fields
     * @param index  index for fields
     * @param type   type of field to check
     * @return true if the field is missing
     */
    private boolean isFieldMissing(String[] fields, Index index, String type) {
        switch (type) {
            case "ms1":
                return Double.parseDouble(fields[index.ms1IntIndex]) == 0;
            case "gene":
                return fields[index.genecIndex].trim().isEmpty();
            default:
                System.err.println("Invalid type: " + type);
                return false;
        }
    }

    private void updateAllGenesList(String[] fields, Index index) {
        String gene = fields[index.genecIndex].trim();
        parameters.allGeneSet.add(gene);
    }

    private void checkAndWarn(int totalCount, int missingCount, File psmFile, String type) {
        if (totalCount == missingCount) {
            switch (type) {
                case "ms1":
                    parameters.ms1Int = false;
                    System.out.println("Warning: All MS1 Intensities in " + psmFile.getAbsolutePath() +
                            " are 0. TMT-Integrator will use summed MS2 reporter ion intensity instead of " +
                            "MS1 ion intensity as reference intensity for abundance calculation.");
                    break;
                case "gene":
                    parameters.geneflag = true;
                    System.out.println("Warning: No gene report is generated because all gene symbols in " +
                            psmFile.getAbsolutePath() + " are missing. The protein ID will be used as gene in other reports.");
                    break;
            }
        }
    }

    /**
     * Process the header line and map column index to field name
     *
     * @param title   header line
     * @param psmFile PSM file
     */
    private void getColumnIndex(String title, File psmFile) {
        String[] columns = title.split("\t");
        Index index = parameters.indMap.getOrDefault(psmFile.getAbsolutePath(), new Index());

        // map column index to field name
        mapColumnIndex(columns, index);
        // find num of channels and check gene flag
        findChannels(columns, index);
        checkReferenceColumns(columns, index);
        adjustRefIndex(columns, index);

        validateColumns(index);
    }

    private void checkReferenceColumns(String[] columns, Index index) {
        int refCount = 0;
        StringBuilder refErrors = new StringBuilder();

        for (int i = index.abnIndex; i < index.abnIndex + parameters.channelNum; i++) {
            String column = columns[i].trim();
            if (column.contains(parameters.refTag)) {
                refCount++;
                refErrors.append(i + 1).append(", ");
                index.refIndex = i;
            }
        }

        if (refCount == 0 && parameters.add_Ref < 0) {
            String message = "TMT-Integrator can't find the reference channel matching \"" + parameters.refTag
                    + "\" from columns " + String.join(",", columns) + ".\n"
                    + "Please check if the reference tag is correctly defined in the parameter file.\n";
            throw new IllegalArgumentException(message);
        }

        if (refCount > 1 && parameters.add_Ref < 0) {
            if (refErrors.length() > 0) {
                refErrors.setLength(refErrors.length() - 2); // remove the last comma and space
            }
            String message = "There are more than one \"" + parameters.refTag + "\" in the columns "
                    + String.join(",", columns) + ".\n"
                    + "Repeated reference tag at column: " + refErrors + ".\n"
                    + "Please make sure the reference tag is unique among all the columns.\n";
            throw new IllegalArgumentException(message);
        }
    }

    private void mapColumnIndex(String[] columns, Index index) {
        for (int i = 0; i < columns.length; i++) {
            String column = columns[i].trim();
            switch (column) {
                case "Spectrum":
                    index.spectrumIndex = i;
                    break;
                case "Peptide":
                    index.pepcIndex = i;
                    index.peptideIndex = i;
                    break;
                case "Modified Peptide":
                    index.modifiedPeptideIndex = i;
                    break;
                case "Extended Peptide":
                    index.extpepIndex = i;
                    break;
                case "Charge":
                    index.chargeIndex = i;
                    break;
                case "Retention":
                    index.rtIndex = i;
                    break;
                case "Observed M/Z":
                    index.observedMzIndex = i;
                    break;
                case "Calculated Peptide Mass":
                case "Peptide Mass":
                    index.pepMassIndex = i;
                    break;
                case "PeptideProphet Probability":
                case "Probability":
                    index.pepProbcIndex = i;
                    break;
                case "Number of Enzymatic Termini":
                    index.numEnzyTermi = i;
                    break;
                case "Protein Start":
                    index.protsIndex = i;
                    break;
                case "Protein End":
                    index.proteIndex = i;
                    break;
                case "Intensity":
                    index.ms1IntIndex = i;
                    break;
                case "Assigned Modifications":
                    index.assignedModcIndex = i;
                    break;
                case "Observed Modifications":
                    index.observedModIndex = i;
                    break;
                case "Purity":
                    index.purityIndex = i;
                    break;
                case "Is Unique":
                    index.isUniquecIndex = i;
                    break;
                case "Protein":
                    index.proteincIndex = i;
                    break;
                case "Protein ID":
                    index.proteinIDcIndex = i;
                    break;
                case "Entry Name":
                    index.entryNameIndex = i;
                    break;
                case "Gene":
                    index.genecIndex = i;
                    break;
                case "Protein Description":
                    index.proteinDescIndex = i;
                    break;
                case "Mapped Genes":
                    index.mapGeneIndex = i;
                    break;
                case "Mapped Proteins":
                    index.mappedProteinsIndex = i;
                    break;
                case "Phospho Site Localization":
                    index.ptmLocalcIndex = i;
                    break;
                case "Total Glycan Composition":
                    index.glycoCompositionIndex = i;
                    break;
                case "Glycan q-value":
                    index.glycoQvalIndex = i;
                    break;
                default:
                    // handle cases require complex evaluation
                    handleComplexCases(column, i, index);
            }
        }
    }

    private void mapColumnIndex(String title, ProteinIndex proteinIndex) {
        // Only mapping necessary columns
        String[] columns = title.split("\t");
        for (int i = 0; i < columns.length; i++) {
            String column = columns[i].trim();
            switch (column) {
                case "Protein ID":
                    proteinIndex.proteinIDIdx = i;
                    break;
                case "Organism":
                    proteinIndex.organismIdx = i;
                    break;
                case "Indistinguishable Proteins":
                    proteinIndex.indistinguishableProteinsIdx = i;
                    break;
            }
        }
    }

    private void handleComplexCases(String column, int columnIdx, Index index) {
        if (!parameters.columntag.isEmpty() && column.contains(parameters.columntag) && !column.contains("Best Localization")) {
            index.ptmLocalcIndex = columnIdx;
        }
    }

    private void findChannels(String[] columns, Index index) {
        // find the offset of the channel columns:
        //   first column name start with "Intensity ", "Resolution ", "SNR "
        for (int i = 0; i < columns.length; i++) {
            if (index.abnIndex < 0 &&columns[i].startsWith("Intensity ")) {
                index.abnIndex = i;
            }
            if (parameters.addIsobaricFilter) {
                if (index.resOffset < 0 && columns[i].startsWith("Resolution ")) {
                    index.resOffset = i;
                }
                if (index.snrOffset < 0 && columns[i].startsWith("SNR ")) {
                    index.snrOffset = i;
                }
            }
        }

        if (parameters.geneflag) {
            index.genecIndex = index.proteincIndex;
        }
    }

    private void adjustRefIndex(String[] columns, Index index) {
        int t = 0;
        int cnum = 0;
        for (int i = index.abnIndex; i < index.abnIndex + parameters.channelNum; i++) {
            if (notNaColumn(columns[i])) {
                cnum++;
            } else if (i < index.refIndex) {
                t++;
            }
        }
        index.refIndex -= t;
        index.usedChannelNum = cnum;
    }

    private boolean notNaColumn(String value) {
        return !value.trim().equalsIgnoreCase("Intensity NA");
    }

    private void validateColumns(Index index) {
        checkIndex(index.pepcIndex, "Peptide");
        checkIndex(index.pepProbcIndex, "Probability");
        checkIndex(index.assignedModcIndex, "Assigned Modifications");
        if (parameters.minSiteProb > 0) {
            checkIndex(index.ptmLocalcIndex, parameters.columntag);
        }
        checkIndex(index.proteinIDcIndex, "Protein");
        checkIndex(index.genecIndex, "Gene");
        checkIndex(index.isUniquecIndex, "Is Unique");
        checkIndex(index.rtIndex, "Retention");
        checkIndex(index.ms1IntIndex, "Intensity");
        checkIndex(index.purityIndex, "Purity");
        checkIndex(index.peptideIndex, "Peptide");
        checkIndex(index.chargeIndex, "Charge");
        checkIndex(index.observedMzIndex, "Observed M/Z");
        checkIndex(index.pepMassIndex, "Calculated Peptide Mass");
        checkIndex(index.modifiedPeptideIndex, "Modified Peptide");
    }

    private void checkIndex(int index, String columnName) {
        if (index < 0) {
            String message = "";
            if (parameters.minSiteProb > 0) {
                message += "Missing the site localization column.\n";
            }
            message += "Missing '" + columnName + "' column. Please check if the column is in the psm tables.";
            throw new IllegalArgumentException(message);
        }
    }
    // endregion==========================================================

    // region parseAndFilterPsmFiles helper methods==============================
    private double calculateTmtThreshold(int numRecords, List<Double> tmtIntensities) {
        Collections.sort(tmtIntensities);
        int tmtThresholdIndex = (int) (numRecords * parameters.minPercent);
        return tmtIntensities.isEmpty() ? 0 : tmtIntensities.get(tmtThresholdIndex);
    }

    private void collapsePsm(Psm psm, Map<String, List<PsmRecord>> psmMap, double tmtThreshold) {
        Set<String> newModTagSet = new HashSet<>();
        for (PsmRecord psmRecord : psm.getPsmRecords()) {
            processPsmRecord(psmRecord, psmMap, newModTagSet, tmtThreshold);
        }
        parameters.modTagSet.addAll(newModTagSet);
    }

    private void processPsmRecord(PsmRecord psmRecord, Map<String, List<PsmRecord>> psmMap, Set<String> modTags, double tmtThreshold) {
        psmRecord.updateFlags(modTags);
        if (psmRecord.isPassCriteria(tmtThreshold)) {
            psmRecord.setUsed(true);
            String key = psmRecord.generatePsmKey();
            // update psmMap
            List<PsmRecord> psmList = psmMap.computeIfAbsent(key, k -> new ArrayList<>());
            psmList.add(psmRecord);
        } else {
            psmMap.get("NotUsed").add(psmRecord);
        }
    }

    private void selectBestPsm(Map<String, List<PsmRecord>> psmMap) {
        if (!parameters.bestPsm) {
            return;
        }
        for (String key : psmMap.keySet()) {
            if (!key.equals("NotUsed")) {
                List<PsmRecord> psmList = psmMap.get(key);
                int bestPsmIndex = findBestPsmIdx(psmList);
                updatePsmList(psmList, bestPsmIndex);
            }
        }
    }

    private int findBestPsmIdx(List<PsmRecord> psmList) {
        double maxIntensity = 0;
        int bestPsmIndex = -1;
        for (int i = 0; i < psmList.size(); i++) {
            double tmtIntensity = psmList.get(i).getSumTmtIntensity();
            if (tmtIntensity > maxIntensity) {
                maxIntensity = tmtIntensity;
                bestPsmIndex = i;
            }
        }
        return bestPsmIndex;
    }

    private void updatePsmList(List<PsmRecord> psmList, int bestPsmIndex) {
        if (bestPsmIndex < 0) {
            return;
        }

        for (int i = 0; i < psmList.size(); i++) {
            psmList.get(i).setUsed(i == bestPsmIndex);
        }
    }
    // endregion==========================================================

}
