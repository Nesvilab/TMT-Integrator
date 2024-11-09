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
import tmtintegrator.pojo.psm.PsmRecord;

/**
 * Utility class for common arithmetic operations, and psm processing.
 *
 * @author rogerli on 05/2024
 */
public final class Utils {

    private Utils() {
        // private constructor to prevent instantiation
        throw new AssertionError("The Utils class cannot be instantiated");
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
     * Create bins for retention time.
     *
     * @param minRt  minimum retention time
     * @param maxRt  maximum retention time
     * @param binNum number of bins
     * @return map of bin start time to list of PSMs
     */
    public static NavigableMap<Double, List<PsmRecord>> createBins(double minRt, double maxRt, double binNum) {
        NavigableMap<Double, List<PsmRecord>> binMap = new TreeMap<>();
        double binWidth = (maxRt - minRt) / binNum;
        double binStart = minRt;
        while (binStart <= maxRt) {
            binMap.put(binStart, new ArrayList<>());
            binStart += binWidth;
        }

        return binMap;
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
            double[] medians = fileAbundanceMap.getOrDefault(filename, new double[index.totLen]);
            sumAbundance += (medians[index.plexNum] > 0) ? medians[index.plexNum] : globalMinRefInt;
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

    // region helper methods
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
    // endregion
}
