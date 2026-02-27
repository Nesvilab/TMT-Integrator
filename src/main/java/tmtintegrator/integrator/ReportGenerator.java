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

import tmtintegrator.constants.Constants;
import tmtintegrator.constants.GroupBy;
import tmtintegrator.constants.NormType;
import tmtintegrator.constants.ReportType;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.ReportInfo;
import tmtintegrator.pojo.Subplex;
import tmtintegrator.utils.ReportData;
import tmtintegrator.utils.Utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;

/**
 * Generate reports for the integrated data
 */
public class ReportGenerator {
    private final Parameters parameters;
    private final ReportData reportData;
    private final GroupBy groupBy;
    private final NormType normType;
    private final List<Map<String, Map<String, double[]>>> plexAbundanceMaps;

    public ReportGenerator(Parameters parameters, ReportData reportData, GroupBy groupBy, NormType normType,
                           List<Map<String, Map<String, double[]>>> plexAbundanceMaps) {
        this.parameters = parameters;
        this.reportData = reportData;
        this.groupBy = groupBy;
        this.normType = normType;
        this.plexAbundanceMaps = plexAbundanceMaps;
    }

    public void generateReports() throws IOException {
        if (parameters.abn_type == 0) {
            if (normType == NormType.ALL_NORM || normType == NormType.SL_IRS) {
                report(ReportType.RATIO_ABUNDANCE);
            } else {
                report(ReportType.RATIO);
                report(ReportType.RATIO_ABUNDANCE);
            }
        } else {
            report(ReportType.RAW_ABUNDANCE);
        }
    }

    // region helper methods
    private void report(ReportType type) throws IOException {
        String groupTag = getGroupTag();
        String protNormTag = getProtNormTag();

        Files.createDirectories(Paths.get(parameters.reportPath));

        String reportPath = generateReportPath(type, groupTag, protNormTag);

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(reportPath))) {
            writeTitle(writer);
            writeRatiosOrAbundances(writer, type);
        }
    }

    private String getGroupTag() {
        switch (groupBy) {
            case GENE: return "gene";
            case PROTEIN_ID: return "protein";
            case PEPTIDE: return "peptide";
            case MULTI_PHOSPHO_SITE: return "multi-site";
            case SINGLE_PHOSPHO_SITE: return "single-site";
            case MULTI_MASS_GLYCO: return "multi-mass";
            case MODIFIED_PEPTIDE: return "modified-peptide";
            default: throw new RuntimeException("Invalid groupBy: " + groupBy);
        }
    }

    private String getProtNormTag() {
        switch (normType) {
            case NONE: return "None";
            case MC: return "MD";
            case GN: return "GN";
            case SL_IRS: return "SL+IRS";
            default: throw new RuntimeException("Invalid protNorm: " + normType);
        }
    }

    private String generateReportPath(ReportType type, String groupTag, String protNormTag) {
        switch (type) {
            case RATIO:
                return parameters.reportPath + File.separator + "ratio_" + groupTag + "_" + protNormTag + ".tsv";
            case RATIO_ABUNDANCE:
                return parameters.reportPath + File.separator + "abundance_" + groupTag + "_" + protNormTag + ".tsv";
            case RAW_ABUNDANCE:
                return parameters.reportPath + File.separator + "raw2abundance_" + groupTag + "_" + protNormTag + ".tsv";
            default:
                throw new IllegalArgumentException("Invalid report type: " + type);
        }
    }

    private void writeTitle(BufferedWriter writer) throws IOException {
        switch (groupBy) {
            case PROTEIN_ID:
                writer.write(ReportInfo.HEADER_PROTEIN);
                break;
            case PEPTIDE:
            case MULTI_PHOSPHO_SITE:
            case SINGLE_PHOSPHO_SITE:
            case MULTI_MASS_GLYCO:
                writer.write(ReportInfo.HEADER);
                break;
            case GENE:
                writer.write("Index\tNumberPSM\tProteinID\tMaxPepProb");
                break;
            case MODIFIED_PEPTIDE:
                writer.write(ReportInfo.HEADER_MODIFIED_PEPTIDE);
                break;
            default:
                throw new IllegalArgumentException("Invalid groupBy: " + groupBy);
        }
        writer.write("\tReferenceIntensity");

        for (String fileName : parameters.fNameLi) {
            Index index = parameters.indMap.get(fileName);
            String[] titles = parameters.titleMap.get(fileName).split("\t");
            for (int p = 0; p < parameters.subplexes.size(); p++) {
                Index.SubplexIndex si = index.subplexIndices.get(p);
                Subplex subplex = parameters.subplexes.get(p);
                writeChannelsTitle(writer, si.abnIndex, si.abnIndex + si.usedChannelNum, subplex.refTag, titles);
            }
        }

        if (parameters.print_RefInt) {
            for (String fileName : parameters.fNameLi) {
                Index index = parameters.indMap.get(fileName);
                File file = new File(fileName);
                String[] titles = parameters.titleMap.get(fileName).split("\t");
                String parentDir = file.getParent().substring(file.getParent().lastIndexOf(File.separator) + 1);
                for (int p = 0; p < parameters.subplexes.size(); p++) {
                    Index.SubplexIndex si = index.subplexIndices.get(p);
                    writer.write("\tRefInt_" + (parameters.add_Ref < 0
                            ? titles[si.refIndex].replace(" Abundance", "").replace("Intensity ", "") : parentDir));
                }
            }
        }
        writer.newLine();
    }

    private void writeChannelsTitle(BufferedWriter writer, int startIdx, int endIdx, String refTag, String[] titles) throws IOException {
        for (int i = startIdx; i < endIdx; i++) {
            if (!titles[i].contains(refTag)) {
                writer.write("\t" + titles[i].replace(" Abundance", "").replace("Intensity ", ""));
            }
        }
    }

    private void writeRatiosOrAbundances(BufferedWriter writer, ReportType type) throws IOException {
        double globalMinRefInt = Double.MAX_VALUE;
        for (Map<String, Map<String, double[]>> plexMap : plexAbundanceMaps) {
            globalMinRefInt = Math.min(globalMinRefInt, Utils.calculateGlobalMinRefInt(plexMap));
        }

        // Use first plex map as the key source (all plexes share the same group keys)
        Map<String, Map<String, double[]>> firstPlexMap = plexAbundanceMaps.get(0);

        for (String groupKey : firstPlexMap.keySet()) {
            boolean isPrint = true;

            double avgAbundance = Utils.calculateAvgAbundance(
                    plexAbundanceMaps.stream()
                            .map(m -> m.getOrDefault(groupKey, Map.of()))
                            .collect(java.util.stream.Collectors.toList()),
                    globalMinRefInt, parameters);

            String[] keyParts = groupKey.split("\t");

            double maxPepProb;
            if (groupBy == GroupBy.MODIFIED_PEPTIDE) {
                maxPepProb = Double.parseDouble(keyParts[keyParts.length - 2]);
            } else {
                maxPepProb = Double.parseDouble(keyParts[keyParts.length - 1]);
            }
            if (maxPepProb < parameters.max_pep_prob_thres || avgAbundance <= 0) {
                isPrint = false;
            }

            String displayKey = groupKey;
            if (groupBy == GroupBy.SINGLE_PHOSPHO_SITE) {
                displayKey = updateGroupKeyForSingleSite(groupKey);
            } else if (groupBy == GroupBy.MULTI_PHOSPHO_SITE) {
                keyParts[4] = keyParts[4].replace(".", "");
                displayKey = String.join("\t", keyParts);
            }

            if (isPrint) {
                if (groupBy == GroupBy.GENE) {
                    writer.write(displayKey.replace("%", "_"));
                } else {
                    reportData.writeReportInfo(writer, displayKey, groupBy);
                }
                writeData(writer, type, avgAbundance, groupKey);
                writer.newLine();
            }
        }
    }

    private String updateGroupKeyForSingleSite(String groupKey) {
        String[] keyParts = groupKey.split("\t");
        Matcher matcher = Constants.KEY_PATTERN.matcher(keyParts[0]);
        int site;
        if (matcher.matches()) {
            site = Integer.parseInt(matcher.group(3));
        } else {
            throw new IllegalArgumentException("Error: cannot parse site from " + keyParts[0]);
        }

        int offset = site - Integer.parseInt(keyParts[5]);
        String pepSeq = keyParts[4].replace(".", "");

        int dotIndex = keyParts[4].indexOf(".");
        if (dotIndex == -1) {
            throw new IllegalArgumentException("Error: cannot find dot in " + keyParts[4]);
        }

        int startIndex = dotIndex + offset - 7;
        int endIndex = dotIndex + offset + 8;

        if (startIndex < 0) {
            pepSeq = "_".repeat(-startIndex) + pepSeq;
            endIndex -= startIndex;
            startIndex = 0;
        }
        if (endIndex > pepSeq.length()) {
            pepSeq += "_".repeat(endIndex - pepSeq.length());
        }

        String sequenceWindow = pepSeq.substring(startIndex, startIndex + 7)
                + Character.toLowerCase(pepSeq.charAt(startIndex + 7))
                + pepSeq.substring(startIndex + 8, endIndex);
        keyParts[4] = sequenceWindow;
        return String.join("\t", keyParts);
    }

    private void writeData(BufferedWriter writer, ReportType type, double avgAbundance, String groupKey) throws IOException {
        writeRefInt(writer, type, avgAbundance);

        StringBuilder abnBuilder = new StringBuilder();
        for (String fileName : parameters.fNameLi) {
            Index index = parameters.indMap.get(fileName);

            for (int p = 0; p < plexAbundanceMaps.size(); p++) {
                Map<String, Map<String, double[]>> plexMap = plexAbundanceMaps.get(p);
                Index.SubplexIndex si = index.subplexIndices.get(p);
                int refIdx = parameters.add_Ref < 0 ? si.refIndex - si.abnIndex : -1;

                double[] nanArray = new double[si.usedChannelNum + 1];
                Arrays.fill(nanArray, Double.NaN);

                Map<String, double[]> fileMap = plexMap.getOrDefault(groupKey, Map.of());
                double[] medians = fileMap.getOrDefault(fileName, nanArray);

                recordRefAbundances(abnBuilder, medians, type);
                writeValues(writer, medians, refIdx, type, avgAbundance);
            }
        }

        if (parameters.print_RefInt) {
            writer.write(abnBuilder.toString());
        }
    }

    private void writeRefInt(BufferedWriter writer, ReportType type, double avgAbundance) throws IOException {
        if (normType != NormType.SL_IRS) {
            if (type == ReportType.RAW_ABUNDANCE) {
                writer.write("\t" + avgAbundance);
            } else {
                writer.write("\t" + (parameters.log2transformed ? Utils.log2(avgAbundance) : avgAbundance));
            }
        } else {
            writer.write("\t" + avgAbundance);
        }
    }

    private void recordRefAbundances(StringBuilder abnBuilder, double[] medians, ReportType type) {
        if (medians[medians.length - 1] > 0) {
            double value = medians[medians.length - 1];
            if (normType != NormType.SL_IRS) {
                if (type == ReportType.RAW_ABUNDANCE) {
                    abnBuilder.append("\t").append(value);
                } else {
                    abnBuilder.append("\t").append(parameters.log2transformed ? Utils.log2(value) : value);
                }
            } else {
                abnBuilder.append("\t").append(value);
            }
        } else {
            abnBuilder.append("\tNA");
        }
    }

    private void writeValues(BufferedWriter writer, double[] medians, int refIndex,
                             ReportType type, double avgAbundance) throws IOException {
        if (normType != NormType.SL_IRS) {
            switch (type) {
                case RATIO:
                    writeRatioValues(writer, medians, refIndex);
                    break;
                case RATIO_ABUNDANCE:
                    writeRatioAbnValues(writer, medians, refIndex, avgAbundance);
                    break;
                case RAW_ABUNDANCE:
                    writeAbnValues(writer, medians, refIndex);
                    break;
                default:
                    throw new IllegalArgumentException("Invalid report type: " + type);
            }
        } else {
            writeAbnValues(writer, medians, refIndex);
        }
    }

    private void writeRatioValues(BufferedWriter writer, double[] medians, int refIndex) throws IOException {
        for (int i = 0; i < medians.length - 1; i++) {
            if (i != refIndex) {
                if (Double.isNaN(medians[i])) {
                    writer.write("\tNA");
                } else {
                    writer.write("\t" + (parameters.log2transformed ? medians[i] : Utils.pow2(medians[i])));
                }
            }
        }
    }

    private void writeRatioAbnValues(BufferedWriter writer, double[] medians, int refIndex,
                                     double avgAbundance) throws IOException {
        for (int i = 0; i < medians.length - 1; i++) {
            if (i != refIndex) {
                if (Double.isNaN(medians[i])) {
                    writer.write("\tNA");
                } else {
                    double abn = parameters.log2transformed ? Utils.log2(Utils.pow2(medians[i]) * avgAbundance)
                            : Utils.pow2(medians[i]) * avgAbundance;
                    writer.write("\t" + abn);
                }
            }
        }
    }

    private void writeAbnValues(BufferedWriter writer, double[] medians, int refIndex) throws IOException {
        for (int i = 0; i < medians.length - 1; i++) {
            if (i != refIndex) {
                if (Double.isNaN(medians[i])) {
                    writer.write("\tNA");
                } else {
                    writer.write("\t" + medians[i]);
                }
            }
        }
    }

    // endregion
}
