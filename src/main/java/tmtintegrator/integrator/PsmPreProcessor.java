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
        for (File psmFile : parameters.fileList) {
            try {
                Psm psm = new Psm(parameters, psmFile);
                List<Double> tmtIntensities = new ArrayList<>();
                psm.readPsmFile(psmFile, tmtIntensities);
                psm.adjustPlexChannelOffsets();

                double tmtThreshold = calculateTmtThreshold(psm.getPsmRecords().size(), tmtIntensities);

                Map<String, List<PsmRecord>> psmMap = new HashMap<>();
                psmMap.put("NotUsed", new ArrayList<>());
                collapsePsm(psm, psmMap, tmtThreshold);

                selectBestPsm(psmMap);
                psm.filterUnUsedPsm(psmMap);

                // normalize each plex's channels
                normalizeData(psm);

                if (parameters.ms1Int) {
                    psm.useMS1Intensity();
                }

                psmList.add(psm);
            } catch (IOException e) {
                System.err.println("Error processing PSM file: " + psmFile.getAbsolutePath());
                e.printStackTrace();
            }
        }
        return psmList;
    }

    // region checkPsmAndBuildIndex helper methods
    private void checkPsmFile(File psmFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(psmFile))) {
            String title = reader.readLine();

            if (title.contains("Resolution ") && title.contains("SNR ")) {
                parameters.addIsobaricFilter = true;
            } else if (title.contains("Resolution ") || title.contains("SNR ")) {
                System.err.println("Error: Both Resolution and SNR columns should be present in the PSM file.");
                System.exit(1);
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
        }
    }

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

    private void getColumnIndex(String title, File psmFile) {
        String[] columns = title.split("\t");
        Index index = parameters.indMap.getOrDefault(psmFile.getAbsolutePath(), new Index());

        mapColumnIndex(columns, index);
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
        for (int i = 0; i < columns.length; i++) {
            if (columns[i].startsWith("Intensity ") && index.subplexIndices.isEmpty()) {
                // compute abnIndex for each plex cumulatively
                int offset = i;
                for (Subplex subplex : parameters.subplexes) {
                    Index.SubplexIndex si = new Index.SubplexIndex();
                    si.abnIndex = offset;
                    index.subplexIndices.add(si);
                    offset += subplex.channelCount;
                }
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

    private void adjustPlexRefIndex(String[] columns, Index index, int plexIdx) {
        Index.SubplexIndex si = index.subplexIndices.get(plexIdx);
        int channelCount = parameters.subplexes.get(plexIdx).channelCount;
        int t = 0;
        int cnum = 0;
        int endIdx = si.abnIndex + channelCount;
        for (int i = si.abnIndex; i < endIdx; i++) {
            if (notNaColumn(columns[i])) {
                cnum++;
            } else if (i < si.refIndex) {
                t++;
            }
        }
        si.refIndex -= t;
        si.usedChannelNum = cnum;
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
