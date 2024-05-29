package tmtintegrator.integrator;

import tmtintegrator.pojo.PsmEntry;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.utils.Utils;

import java.io.*;
import java.util.*;

public class PsmPreProcessor {

    private final Parameters parameters;
    private final List<String> proteinList; // TODO: no usage

    // TODO: Just an example here to avoid magic numbers
    //   Better to rewrite parameters and other data structures to use enums or constants
    private enum ReferenceType {
        SUMMATION,
        AVERAGE,
        MEDIAN;

        public static ReferenceType fromInteger(int value) {
            switch (value) {
                case 0:
                    return SUMMATION;
                case 1:
                    return AVERAGE;
                case 2:
                    return MEDIAN;
                default:
                    throw new IllegalArgumentException("Invalid reference type value: " + value);
            }
        }
    }

    public PsmPreProcessor(Parameters parameters) {
        this.parameters = parameters;
        this.proteinList = new ArrayList<>();
    }

    /**
     * Check PSM tables for missing values <br>
     * Update protein list <br>
     * Get all genes <br>
     * Create index for each PSM file
     *
     * @throws IOException if an I/O error occurs
     */
    public void checkPsmAndBuildIndex() throws IOException {
        Collections.sort(parameters.fNameLi); // TODO: looks unnecessary

        for (File psmFile : parameters.FileLi) {
            checkPsmFile(psmFile);
        }
    }

    /**
     * Update PSM files based on the criteria, select best PSM if required, and print PSM files
     *
     */
    public void updatePsmFiles() {
        for (File psmFile : parameters.FileLi) {
            try {
                Index index = parameters.indMap.get(psmFile.getAbsolutePath());
                // read PSM file and preprocess each line
                List<Double> tmtIntensities = new ArrayList<>();
                List<String> processedLines = readPsmFile(psmFile, index, tmtIntensities);
                double tmtThreshold = calculateTmtThreshold(processedLines, tmtIntensities);
                updateIndex(index);

                // collapse PSM lines based on the criteria
                Map<String, List<String>> psmMap = new TreeMap<>(); // TODO: HashMap?
                psmMap.put("NotUsed", new ArrayList<>());
                collapsePsmLines(processedLines, psmMap, tmtThreshold, index);

                // select best PSM if required
                selectBestPsm(psmMap, index);

                // print intermediate psm file
                printPsmFile(psmFile, index, psmMap);

            } catch (IOException e) {
                System.err.println("Error processing PSM file: " + psmFile.getAbsolutePath());
                e.printStackTrace();
            }
        }
    }

    // region checkPsmAndBuildIndex helper methods=============================================
    private void checkPsmFile(File psmFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(psmFile))) {
            // process header and column index
            String title = reader.readLine();
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
                updateProteins(fields, index);
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
                // TODO: log error
                return false; // This should never happen
        }
    }

    private void updateAllGenesList(String[] fields, Index index) {
        String gene = fields[index.genecIndex].trim();
        if (!parameters.AllGeneLi.contains(gene)) {
            parameters.AllGeneLi.add(gene);
        }
    }

    private void updateProteins(String[] fields, Index index) {
        if (!proteinList.contains(fields[index.proteincIndex])) {
            proteinList.add(fields[index.proteincIndex]);
        }
        if (!parameters.ppMap.containsKey(fields[index.proteinIDcIndex])) {
            parameters.ppMap.put(fields[index.proteinIDcIndex], fields[index.proteincIndex]);
        }
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

        checkReferenceColumns(columns, index);
        // map column index to field name
        mapColumnIndex(columns, index);
        // find num of channels and check gene flag
        findChannels(columns, index);
        adjustRefIndex(columns, index);

        validateColumns(index);
    }

    private void checkReferenceColumns(String[] columns, Index index) {
        int refCount = 0;
        StringBuilder refErrors = new StringBuilder();

        for (int i = 0; i < columns.length; i++) {
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
                case "PeptideProphet Probability":
                case "Probability":
                    index.pepProbcIndex = i;
                    break;
                case "Peptide":
                    index.pepcIndex = i;
                    index.peptideIndex = i;
                    break;
                case "Assigned Modifications":
                    index.assignedModcIndex = i;
                    break;
                case "Phospho Site Localization":
                    index.ptmLocalcIndex = i;
                    break;
                case "Protein ID":
                    index.proteinIDcIndex = i;
                    break;
                case "Protein":
                    index.proteincIndex = i;
                    break;
                case "Gene":
                    index.genecIndex = i;
                    break;
                case "Is Unique":
                    index.isUniquecIndex = i;
                    break;
                case "Retention":
                    index.rtIndex = i;
                    break;
                case "Intensity":
                    index.ms1IntIndex = i;
                    break;
                case "Purity":
                    index.purityIndex = i;
                    break;
                case "Charge":
                    index.chargeIndex = i;
                    break;
                case "Observed M/Z":
                    index.observedMzIndex = i;
                    break;
                case "Calculated Peptide Mass":
                case "Peptide Mass":
                    index.pepMassIndex = i;
                    break;
                case "Mapped Genes":
                    index.mapGeneIndex = i;
                    break;
                case "Modified Peptide":
                    index.modifiedPeptideIndex = i;
                    break;
                case "Number of Enzymatic Termini":
                    index.numEnzyTermi = i;
                    break;
                case "Total Glycan Composition":
                    index.glycoCompositionIndex = i;
                    break;
                case "Glycan q-value":
                    index.glycoQvalIndex = i;
                    break;
                case "Observed Modifications":
                    index.observedModIndex = i;
                    break;
                case "Protein Start":
                    index.protsIndex = i;
                    break;
                case "Protein End":
                    index.proteIndex = i;
                    break;
                case "Extended Peptide":
                    index.extpepIndex = i;
                    break;
                default:
                    // handle cases require complex evaluation
                    handleComplexCases(column, i, index);
            }
        }
    }

    private void handleComplexCases(String column, int columnIdx, Index index) {
        if (!parameters.columntag.isEmpty() && column.contains(parameters.columntag) && !column.contains("Best Localization")) {
            index.ptmLocalcIndex = columnIdx;
        }
    }

    private void findChannels(String[] columns, Index index) {
        index.abnIndex = columns.length - parameters.channelNum;
        index.flength = columns.length;

        if (parameters.geneflag) {
            index.genecIndex = index.proteincIndex;
        }
    }

    private void adjustRefIndex(String[] columns, Index index) {
        int t = 0;
        int cnum = 0;
        for (int i = index.abnIndex; i < columns.length; i++) {
            if (notNaColumn(columns[i])) {
                cnum++;
            } else if (i < index.refIndex) {
                t++;
            }
        }
        index.refIndex -= t;
        index.totLen = (parameters.add_Ref < 0) ? cnum + 1 : cnum + 2;
        index.plexNum = (parameters.add_Ref < 0) ? cnum : cnum + 1;
        index.refIndex = (parameters.add_Ref < 0) ? index.refIndex : index.abnIndex + cnum;
    }

    private boolean notNaColumn(String value) {
        return !value.trim().equalsIgnoreCase("na");
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

    // region updatePsmFiles helper methods==============================
    /**
     * Read PSM file and preprocess each line
     *
     * @param psmFile PSM file
     * @param index   index for PSM file
     * @return list of processed PSM lines
     * @throws IOException if an I/O error occurs
     */
    private List<String> readPsmFile(File psmFile, Index index, List<Double> tmtIntensities) throws IOException {
        List<String> allPsmLines = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(psmFile))) {
            reader.readLine(); // skip header
            String line;
            while ((line = reader.readLine()) != null) {
                allPsmLines.add(processLine(line, index, tmtIntensities));
            }
        }
        return allPsmLines;
    }

    private String processLine(String line, Index index, List<Double> tmtIntensities) {
        String[] fields = line.split("\t");
        StringBuilder newPsm = new StringBuilder();
        double sumTmtIntensity = 0;

        for (int i = 0; i < index.abnIndex; i++) {
            newPsm.append(fields[i]).append("\t");
        }
        newPsm.append("false\t"); // TODO: Assuming "false" for isUsed column
        // append the rest of the fields
        for (int i = index.abnIndex; i < fields.length; i++) {
            newPsm.append(fields[i]).append("\t");
        }

        for (int i = index.abnIndex; i < index.flength; i++) {
            try {
                sumTmtIntensity += Double.parseDouble(fields[i]);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Invalid TMT intensity value: " + fields[i]);
            }
        }
        tmtIntensities.add(sumTmtIntensity);
        newPsm.append(sumTmtIntensity); // TODO: Assuming sum is appended at the end

        return newPsm.toString();
    }

    private double calculateTmtThreshold(List<String> processedLines, List<Double> tmtIntensities) {
        Collections.sort(tmtIntensities);
        int tmtThresholdIndex = (int) (processedLines.size() * parameters.minPercent);
        return tmtIntensities.isEmpty() ? 0 : tmtIntensities.get(tmtThresholdIndex);
    }

    private void updateIndex(Index index) {
        index.isUsedIndex = index.abnIndex;
        index.abnIndex++;
        index.refIndex++;
        index.flength++;
    }

    /**
     * Collapse PSM lines based on the criteria
     *
     * @param processedLines list of processed PSM lines
     * @param psmMap         map of PSM lines (filtered by criteria)
     * @param tmtThreshold   TMT threshold
     * @param index          index for PSM file
     */
    private void collapsePsmLines(List<String> processedLines, Map<String, List<String>> psmMap, double tmtThreshold, Index index) {
        List<String> newModTagList = new ArrayList<>();
        for (String psm : processedLines) {
            processPsmEntry(psm, psmMap, newModTagList, tmtThreshold, index);
        }
        parameters.modTagLi.addAll(newModTagList);
    }

    private void processPsmEntry(String psm, Map<String, List<String>> psmMap, List<String> newModTagList, double tmtThreshold, Index index) {
        PsmEntry psmEntry = new PsmEntry(psm, parameters, index);
        psmEntry.parsePsmEntry();
        psmEntry.checkConfigurations(newModTagList);
        if (psmEntry.isPassCriteria(tmtThreshold)) {
            String newPsm = psmEntry.getProcessedPsm();
            String key = psmEntry.generatePsmKey();
            // update psmMap
            List<String> psmList = psmMap.getOrDefault(key, new ArrayList<>());
            psmList.add(newPsm);
            psmMap.put(key, psmList);
        } else {
            psmMap.get("NotUsed").add(psm + "\t" + psmEntry.getGeneCategory());
        }
    }

    private void selectBestPsm(Map<String, List<String>> psmMap, Index index) {
        if (!parameters.bestPsm) {
            return;
        }
        for (String key : psmMap.keySet()) {
            if (!key.equals("NotUsed")) {
                List<String> psmList = psmMap.get(key);
                int bestPsmIndex = findBestPsm(psmList);
                updatePsmList(psmList, bestPsmIndex, index);
            }
        }

    }

    private int findBestPsm(List<String> psmList) {
        double maxIntensity = 0;
        int bestPsmIndex = -1;
        for (int i = 0; i < psmList.size(); i++) {
            String[] fields = psmList.get(i).split("\t");
            // Index for TMT intensity: fields.length - 2
            // TODO: magic number
            double tmtIntensity = Double.parseDouble(fields[fields.length - 2]);
            if (tmtIntensity > maxIntensity) {
                maxIntensity = tmtIntensity;
                bestPsmIndex = i;
            }
        }
        return bestPsmIndex;
    }

    private void updatePsmList(List<String> psmList, int bestPsmIndex, Index index) {
        if (bestPsmIndex < 0) {
            return;
        }

        List<String> newPsmList = new ArrayList<>();
        for (int i = 0; i < psmList.size(); i++) {
            String[] fields = psmList.get(i).split("\t");
            fields[index.isUsedIndex] = i == bestPsmIndex ? "true" : "false";
            String newPsm = String.join("\t", fields);
            newPsmList.add(newPsm);
        }
        psmList.clear();
        psmList.addAll(newPsmList);
    }

    /**
     * Print PSMs to intermediate .ti file, add one column "IS Used (TMT-I)"
     * @param psmFile PSM file
     * @param index index for PSM file (column index)
     * @param psmMap map of PSM lines (filtered by criteria)
     * @throws IOException if an I/O error occurs
     */
    private void printPsmFile(File psmFile, Index index, Map<String, List<String>> psmMap) throws IOException {
        String newPath = psmFile.getAbsolutePath().replace(".tsv", ".ti");
        // get title line
        String title = getTitleLine(psmFile);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(newPath))) {
            String newHeaderLine = writeHeader(writer, title, index, psmFile);

            // Split the header for use in filtering columns during printing
            String[] headers = newHeaderLine.split("\t");

            // Print the PSM lines filtered by headers
            for (List<String> psmList : psmMap.values()) {
                for (String psm : psmList) {
                    writePsm(writer, psm, headers, index);
                }
            }
        }

    }

    private String getTitleLine(File psmFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(psmFile))) {
            return reader.readLine();
        }
    }

    private String writeHeader(BufferedWriter writer, String title, Index index, File psmFile) throws IOException {
        String[] titles = title.split("\t");
        StringBuilder headerBuilder = new StringBuilder();

        for (int i = 0; i < index.abnIndex - 1; i++) {
            writer.write(titles[i] + "\t");
            headerBuilder.append(titles[i]).append("\t");
        }

        writer.write("Is Used (TMT-I)\t");
        headerBuilder.append("Is Used (TMT-I)\t");

        for (int i = index.abnIndex - 1; i < (index.abnIndex + parameters.channelNum - 1); i++) {
            headerBuilder.append(titles[i]).append("\t");
            if (notNaColumn(titles[i])) {
                writer.write(titles[i] + "\t");
            }
        }

        if (parameters.add_Ref >= 0) {
            String virtualReference = "Virtual_Reference_" + psmFile.getParentFile().getName();
            writer.write(virtualReference + "\t");
            headerBuilder.append(virtualReference).append("\t");
            parameters.refTag = "Virtual_Reference";
        }
        writer.newLine();
        return headerBuilder.toString();
    }

    private void writePsm(BufferedWriter writer, String psm, String[] headers, Index index) throws IOException {
        String[] fields = psm.split("\t");
        int printNum = index.abnIndex + parameters.channelNum;
        for (int i = 0; i < printNum; i++) {
            if (notNaColumn(headers[i])) {
                writer.write(fields[i] + "\t");
            }
        }

        if (parameters.add_Ref >= 0) {
            double refAbundance = calculateRefAbundance(fields, headers, index, printNum);
            writer.write(refAbundance + "\t");
        }
        writer.newLine();
    }

    private double calculateRefAbundance(String[] fields, String[] headers, Index index, int printNum) {
        ReferenceType method = ReferenceType.fromInteger(parameters.add_Ref); // TODO: magic number
        switch (method) {
            case SUMMATION:
                return calculateSummation(fields, headers, index, printNum);
            case AVERAGE:
                return calculateAverage(fields, headers, index, printNum);
            case MEDIAN:
                return calculateMedian(fields, headers, index, printNum);
        }
        return 0;
    }

    private double calculateSummation(String[] fields, String[] headers, Index index, int printNum) {
        double refAbundance = 0;
        for (int i = index.abnIndex; i < printNum; i++) {
            if (notNaColumn(headers[i])) {
                double value = Utils.tryParseDouble(fields[i]);
                refAbundance += value > 0 ? value : 0;
            }
        }
        return refAbundance;
    }

    private double calculateAverage(String[] fields, String[] headers, Index index, int printNum) {
        double refAbundance = 0;
        int count = 0;
        for (int i = index.abnIndex; i < printNum; i++) {
            if (notNaColumn(headers[i])) {
                double value = Utils.tryParseDouble(fields[i]);
                if (value > 0) {
                    refAbundance += value;
                    count++;
                }
            }
        }
        return count > 0 ? refAbundance / count : 0;
    }

    private double calculateMedian(String[] fields, String[] headers, Index index, int printNum) {
        List<Double> refAbundances = new ArrayList<>();
        for (int i = index.abnIndex; i < printNum; i++) {
            if (notNaColumn(headers[i])) {
                double value = Utils.tryParseDouble(fields[i]);
                if (value > 0) {
                    refAbundances.add(value);
                }
            }
        }
        return refAbundances.isEmpty() ? 0 : Utils.takeMedian(refAbundances);
    }
    // endregion==========================================================

}
