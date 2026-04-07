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

package tmtintegrator.integrator;

import static tmtintegrator.utils.Utils.myPrint;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.ProteinIndex;
import tmtintegrator.pojo.Subplex;
import tmtintegrator.pojo.psm.Psm;
import tmtintegrator.pojo.psm.PsmRecord;
import tmtintegrator.utils.ReportData;
import tmtintegrator.utils.Utils;

/**
 * Preprocess PSM files, check for missing values, create index for each PSM file
 */
public class PsmPreProcessor {

    private final Parameters parameters;
    private final ReportData reportData;

    public PsmPreProcessor(Parameters parameters, ReportData reportData) {
        this.parameters = parameters;
        this.reportData = reportData;
    }

    public void checkPsmAndBuildIndex() throws IOException {
        Collections.sort(parameters.fNameLi);

        for (File psmFile : parameters.fileList) {
            checkPsmFile(psmFile);
        }

        if (parameters.isHyperplexing) {
            readCompanionPsmFiles();
        }

        myPrint("Unique genes found across all files: " + parameters.allGeneSet.size(), "INFO");
    }

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

    public List<Psm> preprocessPsm() {
        List<Psm> psmList = new ArrayList<>();
        int totalPsmRead = 0;
        int totalPsmRetained = 0;
        for (File psmFile : parameters.fileList) {
            try {
                Psm psm = new Psm(parameters, psmFile);
                List<Double> tmtIntensities = new ArrayList<>();
                psm.readPsmFile(psmFile, tmtIntensities);
                psm.adjustPlexChannelOffsets();

                int psmReadCount = psm.getPsmRecords().size();
                double tmtThreshold = calculateTmtThreshold(psmReadCount, tmtIntensities);

                Map<String, List<PsmRecord>> psmMap = new HashMap<>();
                psmMap.put("NotUsed", new ArrayList<>());
                collapsePsm(psm, psmMap, tmtThreshold);

                selectBestPsm(psmMap);
                psm.filterUnUsedPsm(psmMap);

                int psmRetainedCount = psm.getPsmRecords().size();

                // normalize each plex's channels
                normalizeData(psm);

                if (parameters.ms1Int) {
                    psm.useMS1Intensity();
                }

                psmList.add(psm);

                totalPsmRead += psmReadCount;
                totalPsmRetained += psmRetainedCount;
                myPrint(psmFile.getName() + ": " + psmReadCount + " PSMs read, "
                        + psmRetainedCount + " retained (" + (psmReadCount - psmRetainedCount) + " filtered)"
                        + ", min_percent=" + parameters.minPercent
                        + " (bottom " + String.format("%.0f", parameters.minPercent * 100) + "% TMT intensity cutoff="
                        + String.format("%.1f", tmtThreshold) + ")", "INFO");
            } catch (IOException e) {
                myPrint("Error processing PSM file: " + psmFile.getAbsolutePath(), "ERROR");
                e.printStackTrace();
            }
        }
        myPrint("Total: " + totalPsmRead + " PSMs read, " + totalPsmRetained + " retained across " + psmList.size() + " files", "INFO");
        return psmList;
    }

    // region checkPsmAndBuildIndex helper methods

    /**
     * For hyperplexing, read companion psm.tsv files to populate extra report columns
     * (Entry Name, Protein Description, Mapped Genes, Mapped Proteins) that are absent from psm_pair.tsv.
     */
    private void readCompanionPsmFiles() throws IOException {
        for (File psmPairFile : parameters.fileList) {
            String psmPath = psmPairFile.getAbsolutePath().replace("psm_pair.tsv", "psm.tsv");
            File psmFile = new File(psmPath);
            if (!psmFile.exists()) {
                myPrint("Companion psm.tsv not found: " + psmPath, "WARN");
                continue;
            }
            try (BufferedReader reader = new BufferedReader(new FileReader(psmFile))) {
                String title = reader.readLine();
                Index tempIndex = new Index();
                mapColumnIndex(title.split("\t"), tempIndex);

                String line;
                while ((line = reader.readLine()) != null) {
                    String[] fields = line.split("\t", -1);
                    reportData.updateExtraColumnsFromCompanionPsm(fields, tempIndex);
                }
            }
        }
    }
    private void checkPsmFile(File psmFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(psmFile))) {
            String title = reader.readLine();

            if (parameters.isHyperplexing) {
                // Validate that expected label columns exist
                for (String label : parameters.hyperplexLabelNames) {
                    if (!title.contains(label + " Intensity ")) {
                        throw new IllegalArgumentException("Hyperplexing label '" + label + "' configured but PSM file '"
                                + psmFile.getName() + "' lacks '" + label + " Intensity' columns.");
                    }
                }
            } else {
                if (title.contains("Resolution ") && title.contains("SNR ")) {
                    parameters.addIsobaricFilter = true;
                } else if (title.contains("Resolution ") || title.contains("SNR ")) {
                    myPrint("Both Resolution and SNR columns should be present in the PSM file.", "ERROR");
                    System.exit(1);
                }
            }

            getColumnIndex(title, psmFile);
            Index index = parameters.indMap.get(psmFile.getAbsolutePath());

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

            myPrint(psmFile.getName() + ": " + totalCount + " PSMs"
                    + (ms1MissingCount > 0 ? ", " + ms1MissingCount + " MS1-missing" : "")
                    + (geneMissingCount > 0 ? ", " + geneMissingCount + " gene-missing" : ""), "INFO");
        }
    }

    private boolean isFieldMissing(String[] fields, Index index, String type) {
        switch (type) {
            case "ms1":
                double sum = 0;
                for (int idx : index.labelMs1IntIndices) {
                    if (idx >= 0 && idx < fields.length) {
                        String val = fields[idx].trim();
                        if (!val.isEmpty()) {
                            sum += Double.parseDouble(val);
                        }
                    }
                }
                return sum == 0;
            case "gene":
                return fields[index.genecIndex].trim().isEmpty();
            default:
                myPrint("Invalid type: " + type, "ERROR");
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
                    myPrint("All MS1 Intensities in " + psmFile.getAbsolutePath() +
                            " are 0. TMT-Integrator will use summed MS2 reporter ion intensity instead of " +
                            "MS1 ion intensity as reference intensity for abundance calculation.", "WARN");
                    break;
                case "gene":
                    parameters.geneflag = true;
                    myPrint("No gene report is generated because all gene symbols in " +
                            psmFile.getAbsolutePath() + " are missing. The protein ID will be used as gene in other reports.", "WARN");
                    break;
            }
        }
    }

    private void getColumnIndex(String title, File psmFile) {
        String[] columns = title.split("\t");
        Index index = parameters.indMap.getOrDefault(psmFile.getAbsolutePath(), new Index());

        mapColumnIndex(columns, index);

        // For hyperplexing, map standard indices to their hyperplexing equivalents
        if (parameters.isHyperplexing) {
            index.pepcIndex = index.labelFreeModifiedPeptideIndex;
            for (int idx : index.labelSpectrumIndices) {
                if (idx >= 0) {
                    index.spectrumIndex = idx;
                    break;
                }
            }
        }

        findChannels(columns, index);

        // check reference columns for each plex
        for (int i = 0; i < parameters.subplexes.size(); i++) {
            Index.SubplexIndex si = index.subplexIndices.get(i);
            Subplex subplex = parameters.subplexes.get(i);
            checkReferenceColumns(columns, si, si.abnIndex, si.abnIndex + subplex.channelCount, subplex.refTag);
        }

        // adjust ref indices for NA channels in each plex
        for (int i = 0; i < parameters.subplexes.size(); i++) {
            adjustPlexRefIndex(columns, index, i);
        }

        validateColumns(index);
    }

    private void checkReferenceColumns(String[] columns, Index.SubplexIndex si, int startIdx, int endIdx, String refTag) {
        int refCount = 0;
        StringBuilder refErrors = new StringBuilder();

        for (int i = startIdx; i < endIdx; i++) {
            String column = columns[i].trim();
            if (column.contains(refTag)) {
                refCount++;
                refErrors.append(i + 1).append(", ");
                si.refIndex = i;
            }
        }

        if (refCount == 0 && parameters.add_Ref < 0) {
            throw new IllegalArgumentException("TMT-Integrator can't find the reference channel matching \"" + refTag
                    + "\" from columns " + String.join(",", columns) + ".\n"
                    + "Please check if the reference tag is correctly defined in the parameter file.\n");
        }

        if (refCount > 1 && parameters.add_Ref < 0) {
            if (refErrors.length() > 0) {
                refErrors.setLength(refErrors.length() - 2);
            }
            throw new IllegalArgumentException("There are more than one \"" + refTag + "\" in the columns "
                    + String.join(",", columns) + ".\n"
                    + "Repeated reference tag at column: " + refErrors + ".\n"
                    + "Please make sure the reference tag is unique among all the columns.\n");
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
                    index.labelMs1IntIndices[0] = i;
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
                case "Label Free Modified Peptide":
                    index.labelFreeModifiedPeptideIndex = i;
                    break;
                case "Light Spectrum":
                    index.labelSpectrumIndices[0] = i;
                    break;
                case "Medium Spectrum":
                    index.labelSpectrumIndices[1] = i;
                    break;
                case "Heavy Spectrum":
                    index.labelSpectrumIndices[2] = i;
                    break;
                case "Light Assigned Modifications":
                    index.labelAssignedModIndices[0] = i;
                    break;
                case "Medium Assigned Modifications":
                    index.labelAssignedModIndices[1] = i;
                    break;
                case "Heavy Assigned Modifications":
                    index.labelAssignedModIndices[2] = i;
                    break;
                case "Light Intensity":
                    index.labelMs1IntIndices[0] = i;
                    break;
                case "Medium Intensity":
                    index.labelMs1IntIndices[1] = i;
                    break;
                case "Heavy Intensity":
                    index.labelMs1IntIndices[2] = i;
                    break;
                default:
                    handleComplexCases(column, i, index);
            }
        }
    }

    private void mapColumnIndex(String title, ProteinIndex proteinIndex) {
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
        if (parameters.isHyperplexing) {
            findHyperplexingChannels(columns, index);
        } else {
            findStandardChannels(columns, index);
        }

        if (parameters.geneflag) {
            index.genecIndex = index.proteincIndex;
        }
    }

    private void findStandardChannels(String[] columns, Index index) {
        int resOffset = -1;
        int snrOffset = -1;
        for (int i = 0; i < columns.length; i++) {
            if (columns[i].startsWith("Intensity ") && index.subplexIndices.isEmpty()) {
                int offset = i;
                for (Subplex subplex : parameters.subplexes) {
                    Index.SubplexIndex si = new Index.SubplexIndex();
                    si.abnIndex = offset;
                    index.subplexIndices.add(si);
                    offset += subplex.channelCount;
                }
            }
            if (parameters.addIsobaricFilter) {
                if (resOffset < 0 && columns[i].startsWith("Resolution ")) {
                    resOffset = i;
                }
                if (snrOffset < 0 && columns[i].startsWith("SNR ")) {
                    snrOffset = i;
                }
            }
        }

        for (Index.SubplexIndex si : index.subplexIndices) {
            si.resOffset = resOffset;
            si.snrOffset = snrOffset;
        }
    }

    private void findHyperplexingChannels(String[] columns, Index index) {
        int subplexOffset = 0;
        for (int labelIdx = 0; labelIdx < parameters.hyperplexLabelNames.size(); labelIdx++) {
            String label = parameters.hyperplexLabelNames.get(labelIdx);
            int numSubplexes = parameters.hyperplexLabelSubplexCounts.get(labelIdx);

            // Find base column offsets for this label
            String intensityPrefix = label + " Intensity ";
            String resPrefix = label + " Resolution ";
            String snrPrefix = label + " SNR ";

            int abnBase = -1, resBase = -1, snrBase = -1;
            for (int i = 0; i < columns.length; i++) {
                String col = columns[i].trim();
                if (abnBase < 0 && col.startsWith(intensityPrefix) && col.length() > intensityPrefix.length()) {
                    abnBase = i;
                }
                if (resBase < 0 && col.startsWith(resPrefix)) {
                    resBase = i;
                }
                if (snrBase < 0 && col.startsWith(snrPrefix)) {
                    snrBase = i;
                }
            }

            if (abnBase < 0) {
                throw new IllegalArgumentException("No '" + intensityPrefix + "' columns found for configured label '" + label + "'.");
            }

            // Create SubplexIndex entries with cumulative offsets within this label
            int channelOffset = 0;
            for (int s = 0; s < numSubplexes; s++) {
                Index.SubplexIndex si = new Index.SubplexIndex();
                si.abnIndex = abnBase + channelOffset;
                si.resOffset = resBase >= 0 ? resBase + channelOffset : -1;
                si.snrOffset = snrBase >= 0 ? snrBase + channelOffset : -1;
                index.subplexIndices.add(si);
                channelOffset += parameters.subplexes.get(subplexOffset + s).channelCount;
            }
            subplexOffset += numSubplexes;
        }

        // detect isobaric filter from per-label Resolution/SNR
        if (index.subplexIndices.stream().anyMatch(si -> si.resOffset >= 0 && si.snrOffset >= 0)) {
            parameters.addIsobaricFilter = true;
        }
    }

    private void adjustPlexRefIndex(String[] columns, Index index, int plexIdx) {
        Index.SubplexIndex si = index.subplexIndices.get(plexIdx);
        int channelCount = parameters.subplexes.get(plexIdx).channelCount;
        int t = 0;
        int cnum = 0;
        int endIdx = si.abnIndex + channelCount;
        for (int i = si.abnIndex; i < endIdx; i++) {
            if (!Psm.isNaChannel(columns[i])) {
                cnum++;
            } else if (i < si.refIndex) {
                t++;
            }
        }
        si.rawRefIndex = si.refIndex; // preserve raw index for reading from unadjusted PSM fields
        si.refIndex -= t;
        si.usedChannelNum = cnum;
    }

    private void validateColumns(Index index) {
        if (parameters.isHyperplexing) {
            validateHyperplexingColumns(index);
        } else {
            validateStandardColumns(index);
        }
    }

    private void validateStandardColumns(Index index) {
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
        checkIndex(index.labelMs1IntIndices[0], "Intensity");
        checkIndex(index.purityIndex, "Purity");
        checkIndex(index.chargeIndex, "Charge");
        checkIndex(index.observedMzIndex, "Observed M/Z");
        checkIndex(index.pepMassIndex, "Calculated Peptide Mass");
        checkIndex(index.modifiedPeptideIndex, "Modified Peptide");
    }

    private void validateHyperplexingColumns(Index index) {
        checkIndex(index.labelFreeModifiedPeptideIndex, "Label Free Modified Peptide");
        boolean hasAnyLabelSpectrum = false;
        for (int idx : index.labelSpectrumIndices) {
            if (idx >= 0) {
                hasAnyLabelSpectrum = true;
                break;
            }
        }
        if (!hasAnyLabelSpectrum) {
            throw new IllegalArgumentException("Missing Light/Medium/Heavy Spectrum columns for hyperplexing format.");
        }
        checkIndex(index.pepProbcIndex, "Probability");
        if (parameters.minSiteProb > 0) {
            checkIndex(index.ptmLocalcIndex, parameters.columntag);
        }
        checkIndex(index.proteinIDcIndex, "Protein");
        checkIndex(index.genecIndex, "Gene");
        checkIndex(index.isUniquecIndex, "Is Unique");
        checkIndex(index.rtIndex, "Retention");
        checkIndex(index.purityIndex, "Purity");
        checkIndex(index.observedMzIndex, "Observed M/Z");
        checkIndex(index.pepMassIndex, "Calculated Peptide Mass");
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
    // endregion

    // region parseAndFilterPsmFiles helper methods
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

    private void normalizeData(Psm psm) {
        // normalize each plex's channels independently
        for (int i = 0; i < parameters.subplexes.size(); i++) {
            psm.setActivePlex(i);
            if (parameters.abn_type == 0) {
                Utils.logNormalizeData(psm);
            }
            if (parameters.psmNorm) {
                Utils.rtNormalizeData(psm);
            }
        }
    }
    // endregion
}
