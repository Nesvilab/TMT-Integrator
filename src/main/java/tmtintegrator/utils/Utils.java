package tmtintegrator.utils;

import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.Ratio;

import java.util.*;

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
            return Double.NaN; // TODO: need test
        }
    }

    public static int tryParseInt(String input) {
        try {
            return Integer.parseInt(input);
        } catch (NumberFormatException e) {
            return -1; // FIXME: use Integer.MIN_VALUE
        }
    }

    public static double takeMedian(List<Double> numbers) {
        if (numbers == null || numbers.isEmpty()) {
            return -9999; // FIXME: Double.NaN is better
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
            return -9999; // FIXME: Double.NaN is better
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
    public static NavigableMap<Double, List<String>> createBins(double minRt, double maxRt, double binNum) {
        NavigableMap<Double, List<String>> binMap = new TreeMap<>();
        double binWidth = (maxRt - minRt) / binNum;
        // create bins as empty lists
//        for (int i = 0; i < binNum; i++) {
//            double binStart = minRt + i * binWidth;
//            binMap.put(binStart, new ArrayList<>());
//        }

        // FIXME: This creates one more bin than expected (including a bin starting from maxRt)
        //   just to keep the decimal precision, should be refactored to code above
        double binStart = minRt;
        while (binStart <= maxRt) {
            binMap.put(binStart, new ArrayList<>());
            binStart += binWidth;
        }

        return binMap;
    }

    // FIXME: This method is kept for double decimal precision, should be refactored to use createBins
    public static NavigableMap<Double, List<String>> createBinsToBeRefactor(double minRt, double maxRt, int binNum) {
        NavigableMap<Double, List<String>> binMap = new TreeMap<>();
        //region Generate bins
        double gap = (maxRt - minRt) / binNum;
        double sRt = minRt;
        double eRt = maxRt + gap;
        while (sRt < eRt) {
            binMap.put(sRt, new ArrayList<>());
            sRt += gap;
        }
        binMap.put((maxRt + gap), new ArrayList<>());

        return binMap;
    }

    /**
     * Make sure there are <= 7 amino acids on both sides of the extended peptide sequence.
     *
     * @param extPep extended peptide sequence
     * @return refined extended peptide sequence
     * @author Fengchao
     */
    public static String refineExtendedSequence(String extPep) {
        int firstPeriodIdx = extPep.indexOf(".");
        if (firstPeriodIdx == -1) {
            return extPep;
        }
        int startIdx = Math.max(0, firstPeriodIdx - 7);
        int lastPeriodIdx = extPep.indexOf(".", firstPeriodIdx + 1);
        if (lastPeriodIdx == -1) {
            return extPep.substring(startIdx);
        }
        int endIdx = Math.min(extPep.length(), lastPeriodIdx + 8);
        return extPep.substring(startIdx, endIdx);
    }

    public static double[][] convertTo2DArray(List<String> psmList, Index index) {
        double[][] ratio2DValues = new double[psmList.size()][index.plexNum];
        for (int i = 0; i < psmList.size(); i++) {
            String[] fields = psmList.get(i).split("\t");
            for (int j = index.abnIndex; j < fields.length; j++) {
                try {
                    ratio2DValues[i][j - index.abnIndex] = Double.parseDouble(fields[j]);
                } catch (NumberFormatException e) {
                    ratio2DValues[i][j - index.abnIndex] = Double.NaN;
                }
            }
        }
        return ratio2DValues;
    }

    /**
     * Compute the interquartile range (IQR) of a list of ratios.
     *
     * @param ratios list of ratios
     * @return lower and upper bounds of IQR {lower bound, upper bound}
     */
    public static double[] computeIQR(List<Double> ratios) {
        Collections.sort(ratios);
        // TODO: In old code, q1List range [0, size / 4] and q3List range [size - size / 4 - 1, size - 1]
        //   which seems incorrect according to IQR definition, should be [0, size / 2] and [(size + 1) / 2, size - 1]
        int n = ratios.size() / 4;
        List<Double> q1List = new ArrayList<>(ratios.subList(0, n + 1));
        List<Double> q3List = new ArrayList<>(ratios.subList(ratios.size() - n - 1, ratios.size()));

        double q1 = takeMedian(q1List);
        double q3 = takeMedian(q3List);
        double iqr = q3 - q1;
        return new double[]{q1 - 1.5 * iqr, q3 + 1.5 * iqr};
    }

    /**
     * Update PSM ratios with new values.
     *
     * @param psmList       list of PSMs
     * @param ratio2DValues 2D array of ratios
     * @param index         index of fields
     * @return updated list of PSMs
     */
    public static List<String> updatePsmRatios(List<String> psmList, double[][] ratio2DValues, Index index) {
        List<String> newPsmList = new ArrayList<>();
        for (int i = 0; i < psmList.size(); i++) {
            StringBuilder psmBuilder = new StringBuilder();
            String[] fields = psmList.get(i).split("\t");
            for (int j = 0; j < fields.length; j++) {
                psmBuilder.append(j < index.abnIndex ? (fields[j] + "\t") : (ratio2DValues[i][j - index.abnIndex] + "\t"));
            }
            newPsmList.add(psmBuilder.toString().trim());
        }
        return newPsmList;
    }

    public static double calculateGlobalMedian(Map<String, double[]> medianMap) {
        List<Double> channelValues = new ArrayList<>();
        for (double[] medianValues : medianMap.values()) {
            for (double median : medianValues) {
                if (median != -9999) { // FIXME: !Double.isNaN(median) is better
                    channelValues.add(median);
                }
            }
        }
        return takeMedian(channelValues);
    }

    public static double calculateGlobalMinRefInt(Map<String, Map<String, double[]>> groupAbundanceMap, Parameters parameters) {
        double globalMinRefInt = Double.MAX_VALUE;
        for (Map<String, double[]> fileAbundanceMap : groupAbundanceMap.values()) {
            for (Map.Entry<String, double[]> entry : fileAbundanceMap.entrySet()) {
                String filename = entry.getKey();
                double[] medianValues = entry.getValue();
                Index index = parameters.indMap.get(filename);
                double refInt = medianValues[index.plexNum];

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
        double avgAbundance = 0;
        for (String filename : parameters.fNameLi) {
            Index index = parameters.indMap.get(filename);
            double[] medians = fileAbundanceMap.getOrDefault(filename, new double[index.totLen]);
            avgAbundance += (medians[index.plexNum] > 0) ? medians[index.plexNum] : globalMinRefInt;
        }
        return avgAbundance / parameters.fNameLi.size();
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
            // if using composition for index, read it from glycan composition column. Still get AA site from assigned mods
            mod = inputMod.substring(inputMod.indexOf("(") - 1, inputMod.indexOf("("));
            mod = String.format("%s(%s)", mod, glycanComposition);
        } else {
            // read mass from assigned mod, as would do for any other mod
            mod = inputMod.substring(inputMod.indexOf("(") - 1);
        }
        return mod;
    }

    // region helper methods
    private static void takeWeightsAndNormalize(List<Ratio> ratioList) {
        double sum = 0;
        double pow = 1; // FIXME: should be configurable
        // Precursor intensity with power
        for (Ratio ratio : ratioList) {
            ratio.weight = Math.pow(ratio.preInt, pow); // FIXME: pow == 1?
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
        // Fall back to the ratio with the highest weight (should not happen)
        ratioList.sort(Comparator.comparingDouble(r -> r.weight));
        return ratioList.get(ratioList.size() - 1).ratio;
    }
    // endregion
}
