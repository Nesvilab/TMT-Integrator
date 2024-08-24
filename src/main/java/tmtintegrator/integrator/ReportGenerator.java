package tmtintegrator.integrator;

import tmtintegrator.constants.Constants;
import tmtintegrator.constants.GroupBy;
import tmtintegrator.constants.NormType;
import tmtintegrator.constants.ReportType;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.ReportInfo;
import tmtintegrator.utils.ReportData;
import tmtintegrator.utils.Utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
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
    private final Map<String, Map<String, double[]>> groupAboundanceMap; // <groupKey, <filename, abundance[]>>

    public ReportGenerator(Parameters parameters, ReportData reportData, GroupBy groupBy, NormType normType,
                           Map<String, Map<String, double[]>> groupAboundanceMap) {
        this.parameters = parameters;
        this.reportData = reportData;
        this.groupBy = groupBy;
        this.normType = normType;
        this.groupAboundanceMap = groupAboundanceMap;
    }


    public void generateReports() throws IOException {
        if (parameters.abn_type == 0) {
            if (normType == NormType.ALL_NORM || normType == NormType.SL_IRS) {
                report(ReportType.RATIO_ABUNDANCE); // report abundances
            } else {
                report(ReportType.RATIO); // report ratios
                report(ReportType.RATIO_ABUNDANCE); // report abundances
            }
        } else {
            report(ReportType.RAW_ABUNDANCE); // report abundances
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
            case GENE:
                return "gene";
            case PROTEIN_ID:
                return "protein";
            case PEPTIDE:
                return "peptide";
            case MULTI_PHOSPHO_SITE:
                return "multi-site";
            case SINGLE_PHOSPHO_SITE:
                return "single-site";
            case MULTI_MASS_GLYCO:
                return "multi-mass";
            default:
                throw new RuntimeException("Invalid groupBy: " + groupBy);
        }
    }

    private String getProtNormTag() {
        switch (normType) {
            case NONE:
                return "None";
            case MC:
                return "MD"; // FIXME: why not "MC"?
            case GN:
                return "GN";
            case SL_IRS:
                return "SL+IRS";
            default:
                throw new RuntimeException("Invalid protNorm: " + normType);
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
            default:
                throw new IllegalArgumentException("Invalid groupBy: " + groupBy);
        }
        writer.write("\tReferenceIntensity");

        for (String fileName : parameters.fNameLi) {
            Index index = parameters.indMap.get(fileName);
            String[] titles = parameters.titleMap.get(fileName).split("\t");
            for (int i = index.abnIndex; i < titles.length; i++) {
                if (!titles[i].contains(parameters.refTag)) {
                    writer.write("\t" + titles[i].replace(" Abundance", ""));
                }
            }
        }

        if (parameters.print_RefInt) {
            for (String fileName : parameters.fNameLi) {
                Index index = parameters.indMap.get(fileName);
                File file = new File(fileName);
                String[] titles = parameters.titleMap.get(fileName).split("\t");
                String parentDir = file.getParent().substring(file.getParent().lastIndexOf(File.separator) + 1);
                writer.write("\tRefInt_" + (parameters.add_Ref < 0
                        ? titles[index.refIndex].replace(" Abundance", "") : parentDir));
            }
        }
        writer.newLine();
    }

    private void writeRatiosOrAbundances(BufferedWriter writer, ReportType type) throws IOException {
        double globalMinRefInt = Utils.calculateGlobalMinRefInt(groupAboundanceMap, parameters);
        for (Map.Entry<String, Map<String, double[]>> groupEntry : groupAboundanceMap.entrySet()) {
            String groupKey = groupEntry.getKey();
            Map<String, double[]> fileAbundanceMap = groupEntry.getValue();
            boolean isPrint = true;

            // M2: Calculate average abundance (using global minimum reference intensity)
            double avgAbundance = Utils.calculateAvgAbundance(fileAbundanceMap, globalMinRefInt, parameters);

            String[] keyParts = groupKey.split("\t");
            double maxPepProb = Double.parseDouble(keyParts[keyParts.length - 1]);
            if (avgAbundance <= 0 || maxPepProb < parameters.max_pep_prob_thres) {
                isPrint = false;
            }

            if (groupBy == GroupBy.SINGLE_PHOSPHO_SITE) {
                groupKey = updateGroupKeyForSingleSite(groupKey);
            }

            if (isPrint) {
                if (groupBy == GroupBy.GENE) {
                    writer.write(groupKey.replace("%", "_"));
                } else {
                    // peptide and site level, protein level
                    reportData.writeReportInfo(writer, groupKey, groupBy);
                }
                writeData(writer, type, avgAbundance, fileAbundanceMap);
                writer.newLine();
            }
        }
    }

    private String updateGroupKeyForSingleSite(String groupKey) {
        String[] keyParts = groupKey.split("\t");
        Matcher matcher = Constants.KEY_PATTERN.matcher(keyParts[0]);
        int site;
        if (matcher.matches()) {
            site = Integer.parseInt(matcher.group(3)); // Match the site number, the last match group in the PATTERN (\\d+)
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

    private void writeData(BufferedWriter writer, ReportType type, double avgAbundance,
                           Map<String, double[]> fileAbundanceMap) throws IOException {
        if (normType != NormType.SL_IRS) {
            if (type == ReportType.RAW_ABUNDANCE) {
                writer.write("\t" + avgAbundance);
            } else {
                writer.write("\t" + (parameters.log2transformed ? Utils.log2(avgAbundance) : avgAbundance));
            }
        } else {
            writer.write("\t" + avgAbundance);
        }

        StringBuilder abnBuilder = new StringBuilder();
        for (String fileName : parameters.fNameLi) {
            Index index = parameters.indMap.get(fileName);
            int refIndex = index.refIndex - index.abnIndex;
            double[] nanArray = new double[index.totLen];
            Arrays.fill(nanArray, Double.NaN);
            double[] medians = fileAbundanceMap.getOrDefault(fileName, nanArray);

            // record reference abundances
            recordRefAbundances(abnBuilder, medians, index, type);

            writeValues(writer, medians, refIndex, type, avgAbundance);
        }

        if (parameters.print_RefInt) {
            writer.write(abnBuilder.toString());
        }
    }

    private void recordRefAbundances(StringBuilder abnBuilder, double[] medians, Index index, ReportType type) {
        if (medians[index.plexNum] > 0) {
            double value = medians[index.plexNum];
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
        if (normType != NormType.SL_IRS) { // NONE, MC, GN
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
        } else { // SL, IRS, SL+IRS
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
