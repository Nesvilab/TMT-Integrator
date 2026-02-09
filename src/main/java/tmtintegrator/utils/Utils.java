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

package tmtintegrator.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import tmtintegrator.constants.Constants;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.Ratio;
import tmtintegrator.pojo.psm.Psm;
import tmtintegrator.pojo.psm.PsmRecord;

/**
 * Utility class for common arithmetic operations, and psm processing.
 */
public final class Utils {

    private static final Pattern varModPattern = Pattern.compile("([0-9]+)([A-Z])\\(([0-9.-]+)\\)");
    private static final Pattern nTermModPattern = Pattern.compile("N-term\\(([0-9.-]+)\\)");
    private static final Pattern cTermModPattern = Pattern.compile("C-term\\(([0-9.-]+)\\)");

    private Utils() {
        // private constructor to prevent instantiation
        throw new AssertionError("The Utils class cannot be instantiated");
    }

    public static boolean matchLabels(String assignedModifications, float[] labels, float tol) {
        if (labels.length == 1 && labels[0] == 0) { // if the label mass is unknown, do not filter the PSMs based on the label
            return true;
        }

        float mod;
        Matcher m1 = nTermModPattern.matcher(assignedModifications);
        Matcher m2 = varModPattern.matcher(assignedModifications);
        Matcher m3 = cTermModPattern.matcher(assignedModifications);

        while (m1.find()) {
            mod = Float.parseFloat(m1.group(1));
            if (matchMass(mod, labels, tol)) {
                return true;
            }
        }

        while (m2.find()) {
            mod = Float.parseFloat(m2.group(3));
            if (matchMass(mod, labels, tol)) {
                return true;
            }
        }

        while (m3.find()) {
            mod = Float.parseFloat(m3.group(1));
            if (matchMass(mod, labels, tol)) {
                return true;
            }
        }

        return false;
    }

    private static boolean matchMass(float a, float[] b, float t) {
        for (float f : b) {
            if (Math.abs(a - f) < t) {
                return true;
            }
        }
        return false;
    }

    public static double tryParseDouble(String input) {
        try {
            return Double.parseDouble(input);
        } catch (NumberFormatException e) {
            return Double.NaN;
        }
    }

    public static int tryParseInt(String input) {
        try {
            return Integer.parseInt(input);
        } catch (NumberFormatException e) {
            return Integer.MIN_VALUE;
        }
    }

    public static double takeMedian(List<Double> numbers) {
        if (numbers == null || numbers.isEmpty()) {
            return Double.NaN;
        }

        Collections.sort(numbers);
        int size = numbers.size();
        int middleIndex = size / 2;

        if (size % 2 == 0) {
            return (numbers.get(middleIndex - 1) + numbers.get(middleIndex)) / 2.0;
        } else {
            return numbers.get(middleIndex);
        }
    }

    public static double takeWeightedMedian(List<Ratio> ratioList) {
        if (ratioList == null || ratioList.isEmpty()) {
            return Double.NaN;
        } else if (ratioList.size() == 1) {
            return ratioList.get(0).ratio;
        } else if (ratioList.size() == 2) {
            return (ratioList.get(0).ratio + ratioList.get(1).ratio) / 2.0;
        }

        // 1. Calculate weights and normalize
        takeWeightsAndNormalize(ratioList);

        // 3. Sort ratios by ratio values
        ratioList.sort(Comparator.comparingDouble(r -> r.ratio));

        // 4. Find the weighted median
        return findWeightedMedian(ratioList);
    }

    public static double log2(double x) {
        return Math.log(x) / Math.log(2);
    }

    public static double pow2(double x) {
        return Math.pow(2, x);
    }

    /**
     * Take log and normalize PSM.
     *
     * @param psm     PSM file
     * @param isTmt35 whether the experiment is TMT 35-plex
     */
    public static void logNormalizeData(Psm psm, boolean isTmt35) {
        List<PsmRecord> psmRecords = psm.getPsmRecords();
        // log normalize PSM
        for (PsmRecord psmRecord : psmRecords) {
            logNormalizePsm(psmRecord, isTmt35);
        }
    }

    /**
     * Normalize PSM by retention time.
     *
     * @param psm PSM files
     */
    public static void rtNormalizeData(Psm psm, boolean isTmt35) {
        List<PsmRecord> psmRecords = psm.getPsmRecords();
        Index index = psm.getIndex();
        // rt normalize PSM
        rtNormalizePsm(psmRecords, index, isTmt35);
    }

    public static String[] getLRStrings(String extPep) {
        int firstPeriodIdx = extPep.indexOf(".");
        if (firstPeriodIdx == -1) {
            return new String[]{"", ""};
        }
        int startIdx = Math.max(0, firstPeriodIdx - 7);
        int lastPeriodIdx = extPep.indexOf(".", firstPeriodIdx + 1);
        if (lastPeriodIdx == -1) {
            return new String[]{extPep.substring(startIdx, firstPeriodIdx), ""};
        }
        int endIdx = Math.min(extPep.length(), lastPeriodIdx + 8);
        return new String[]{extPep.substring(startIdx, firstPeriodIdx), extPep.substring(lastPeriodIdx + 1, endIdx)};
    }

    /**
     * Compute the interquartile range (IQR) of a list of ratios.
     *
     * @param ratios list of ratios
     * @return lower and upper bounds of IQR {lower bound, upper bound}
     */
    public static double[] computeIQR(List<Double> ratios) {
        Collections.sort(ratios);
        int n = ratios.size() / 2;
        List<Double> q1List = new ArrayList<>(ratios.subList(0, n));
        List<Double> q3List = new ArrayList<>(ratios.subList(n, ratios.size()));

        double q1 = takeMedian(q1List);
        double q3 = takeMedian(q3List);
        double iqr = q3 - q1;
        return new double[]{q1 - 1.5 * iqr, q3 + 1.5 * iqr};
    }

    public static double calculateGlobalMedian(Map<String, double[]> medianMap) {
        List<Double> channelValues = new ArrayList<>();
        for (double[] medianValues : medianMap.values()) {
            for (double median : medianValues) {
                if (!Double.isNaN(median)) {
                    channelValues.add(median);
                }
            }
        }
        return takeMedian(channelValues);
    }

    public static double calculateGlobalMinRefInt(Map<String, Map<String, double[]>> groupAbundanceMap) {
        double globalMinRefInt = Double.MAX_VALUE;
        for (Map<String, double[]> fileAbundanceMap : groupAbundanceMap.values()) {
            for (Map.Entry<String, double[]> entry : fileAbundanceMap.entrySet()) {
                double[] medianValues = entry.getValue();
                double refInt = medianValues[medianValues.length - 1]; // total reference intensity at the end
                if (refInt > 0 && refInt < globalMinRefInt) {
                    globalMinRefInt = refInt;
                }
            }
        }
        return globalMinRefInt;
    }

    /**
     * M2: Calculate average abundance (using global min reference intensity).
     *
     * @param fileAbundanceMap map of file to abundance values
     * @param globalMinRefInt  global minimum reference intensity
     * @return average abundance
     */
    public static double calculateAvgAbundance(Map<String, double[]> fileAbundanceMap, double globalMinRefInt, Parameters parameters) {
        double sumAbundance = 0;
        for (String filename : parameters.fNameLi) {
            Index index = parameters.indMap.get(filename);
            double[] medians = fileAbundanceMap.getOrDefault(filename, new double[index.usedChannelNum + 1]);
            double totalRefInt = medians[medians.length - 1]; // total reference intensity at the end
            sumAbundance += (totalRefInt > 0) ? totalRefInt : globalMinRefInt;
        }
        return sumAbundance / parameters.fNameLi.size();
    }

    public static double calculateAvgAbundance(Map<String, double[]> fileAbundanceMap,
                                               Map<String, double[]> fileDAbundanceMap,
                                               double globalMinRefInt, Parameters parameters) {
        double sumAbundance = 0;
        for (String filename : parameters.fNameLi) {
            Index index = parameters.indMap.get(filename);
            double[] medians = fileAbundanceMap.getOrDefault(filename, new double[index.usedChannelNum + 1]);
            double totalRefInt = medians[medians.length - 1]; // total reference intensity at the end
            if (parameters.isTmt35) {
                double[] dMedians = fileDAbundanceMap.getOrDefault(filename, new double[index.usedDChannelNum + 1]);
                totalRefInt += dMedians[dMedians.length - 1];
            }
            sumAbundance += (totalRefInt > 0) ? totalRefInt : globalMinRefInt;
        }
        return sumAbundance / parameters.fNameLi.size();
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
     * @author dpolasky
     */
    public static String getAssignedModIndex(String inputMod, String glycanComposition, boolean useGlycanComposition) {
        String mod;
        if (useGlycanComposition && !glycanComposition.isEmpty()) {
            mod = inputMod.substring(inputMod.indexOf("(") - 1, inputMod.indexOf("("));

            // remove extra whitespace
            // original format "modTags % 123.45" -> "modTags%123.45"
            glycanComposition = glycanComposition.replaceAll("\\s+", "");

            mod = String.format("%s(%s)", mod, glycanComposition);
        } else if (inputMod.toLowerCase().startsWith("n-term") || inputMod.toLowerCase().startsWith("c-term")) {
            return inputMod;
        } else {
            mod = inputMod.substring(inputMod.indexOf("(") - 1);
        }
        return mod;
    }

    public static boolean isModTagMatch(String targetMod, Set<String> modTagSet) {
        Matcher nonGlycoMatcher = Constants.MOD_TAG_PATTERN.matcher(targetMod);
        if (nonGlycoMatcher.matches()) {
            String targetSeq = nonGlycoMatcher.group(Constants.AA_GROUP);
            double targetMass = Double.parseDouble(nonGlycoMatcher.group(Constants.MASS_GROUP));
            return checkMatch(targetMod, modTagSet, Constants.MOD_TAG_PATTERN, targetSeq, targetMass);
        } else {
            Matcher glycoMatcher = Constants.GLYCAN_MOD_TAG_PATTERN.matcher(targetMod);
            if (glycoMatcher.matches()) {
                String targetSeq = glycoMatcher.group(Constants.AA_GROUP);
                double targetMass = Double.parseDouble(glycoMatcher.group(Constants.MASS_GROUP));
                return checkMatch(targetMod, modTagSet, Constants.GLYCAN_MOD_TAG_PATTERN, targetSeq, targetMass);
            }
        }
        return false;
    }

    // region helper methods
    private static boolean checkMatch(String targetMod, Set<String> modTagSet, Pattern pattern, String targetSeq, double targetMass) {
        for (String modTag : modTagSet) {
            if (modTag.equals(targetMod)) {
                return true;
            }
            Matcher matcher = pattern.matcher(modTag);
            if (matcher.matches()) {
                String seq = matcher.group(Constants.AA_GROUP);
                double mass = Double.parseDouble(matcher.group(Constants.MASS_GROUP));
                if (seq.equals(targetSeq) && Math.abs(mass - targetMass) <= Constants.MASS_TOLERANCE) {
                    return true;
                }
            }
        }
        return false;
    }

    private static void takeWeightsAndNormalize(List<Ratio> ratioList) {
        double sum = 0;
        double pow = Constants.POWER_NUM;
        // Precursor intensity with power
        for (Ratio ratio : ratioList) {
            ratio.weight = Math.pow(ratio.preInt, pow);
            sum += ratio.weight;
        }
        // Normalize weights
        for (Ratio ratio : ratioList) {
            ratio.weight /= sum;
        }
    }

    private static double findWeightedMedian(List<Ratio> ratioList) {
        double cumulativeWeight = 0.0;
        for (Ratio ratio : ratioList) {
            cumulativeWeight += ratio.weight;
            if (cumulativeWeight >= 0.5) {
                return ratio.ratio;
            }
        }
        ratioList.sort(Comparator.comparingDouble(r -> r.weight));
        return ratioList.get(ratioList.size() - 1).ratio;
    }

    private static void logNormalizePsm(PsmRecord psmRecord, boolean isTmt35) {
        double refIntensity = psmRecord.getRefIntensity();
        double logRefInt = refIntensity > 0 ? Utils.log2(refIntensity) : 0;
        List<Double> channels = psmRecord.getChannels();
        for (int i = 0; i < channels.size(); i++) {
            double intensity = channels.get(i);
            intensity = intensity > 0 ? Utils.log2(intensity) - logRefInt : Double.NaN;
            channels.set(i, intensity);
        }

        if (isTmt35) {
            double refDIntensity = psmRecord.getDRefIntensity();
            double logRefDInt = refDIntensity > 0 ? Utils.log2(refDIntensity) : 0;
            List<Double> dChannels = psmRecord.getDChannels();
            for (int i = 0; i < dChannels.size(); i++) {
                double intensity = dChannels.get(i);
                intensity = intensity > 0 ? Utils.log2(intensity) - logRefDInt : Double.NaN;
                dChannels.set(i, intensity);
            }
        }
    }

    private static void rtNormalizePsm(List<PsmRecord> psmRecords, Index index, boolean isTmt35) {
        double minRt = Double.MAX_VALUE;
        double maxRt = Double.MIN_VALUE;

        // find min and max retention time
        for (PsmRecord psmRecord : psmRecords) {
            double rt = psmRecord.getRetention();
            if (rt < minRt) {
                minRt = rt;
            }
            if (rt > maxRt) {
                maxRt = rt;
            }
        }

        NavigableMap<Double, List<PsmRecord>> binMap = createBins(minRt, maxRt, Constants.BIN_NUM);

        // assign PSMs to bins
        for (PsmRecord psmRecord : psmRecords) {
            double rt = psmRecord.getRetention();
            double binStart = binMap.floorKey(rt); // find the bin start time
            binMap.get(binStart).add(psmRecord);
        }

        // normalize PSMs in each bin
        for (Map.Entry<Double, List<PsmRecord>> entry : binMap.entrySet()) {
            List<PsmRecord> binPsmList = entry.getValue();
            normalizePsmList(binPsmList, index, isTmt35);
        }
    }

    /**
     * Create bins for retention time.
     *
     * @param minRt  minimum retention time
     * @param maxRt  maximum retention time
     * @param binNum number of bins
     * @return map of bin start time to list of PSMs
     */
    private static NavigableMap<Double, List<PsmRecord>> createBins(double minRt, double maxRt, double binNum) {
        NavigableMap<Double, List<PsmRecord>> binMap = new TreeMap<>();
        double binWidth = (maxRt - minRt) / binNum;
        double binStart = minRt;
        while (binStart <= maxRt) {
            binMap.put(binStart, new ArrayList<>());
            binStart += binWidth;
        }

        return binMap;
    }

    private static void normalizePsmList(List<PsmRecord> psmRecords, Index index, boolean isTmt35) {
        // calculate median ratio for each channel
        double[] medianValues = calculateMedianByChannel(psmRecords, index, false);
        // subtract median ratio from each channel
        adjustRatiosByMedian(psmRecords, medianValues, false);

        if (isTmt35) {
            double[] medianDValues = calculateMedianByChannel(psmRecords, index, true);
            adjustRatiosByMedian(psmRecords, medianDValues, true);
        }
    }

    private static double[] calculateMedianByChannel(List<PsmRecord> psmRecords, Index index, boolean isTmt35) {
        int size = isTmt35 ? index.usedDChannelNum : index.usedChannelNum;
        double[] medianValues = new double[size];
        for (int j = 0; j < size; j++) {
            List<Double> channelValues = new ArrayList<>();
            for (PsmRecord psmRecord : psmRecords) {
                double intensity = isTmt35 ? psmRecord.getDChannels().get(j) : psmRecord.getChannels().get(j);
                if (!Double.isNaN(intensity)) {
                    channelValues.add(intensity);
                }
            }
            medianValues[j] = Utils.takeMedian(channelValues);
        }
        return medianValues;
    }

    private static void adjustRatiosByMedian(List<PsmRecord> psmRecords, double[] medianValues, boolean isTmt35) {
        for (PsmRecord psmRecord : psmRecords) {
            List<Double> channels = isTmt35 ? psmRecord.getDChannels() : psmRecord.getChannels();
            for (int j = 0; j < channels.size(); j++) {
                double intensity = channels.get(j);
                if (!Double.isNaN(intensity)) {
                    intensity -= medianValues[j];
                    channels.set(j, intensity);
                }
            }
        }
    }
    // endregion
}
