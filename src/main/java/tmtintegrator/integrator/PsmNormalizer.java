package tmtintegrator.integrator;

import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.utils.Utils;

import java.util.*;

/**
 * Utility class to normalize PSMs.
 *
 * @author rogerli on 05/2024
 */

public class PsmNormalizer {

    private final Parameters parameters;
    public static final int BIN_NUM = 10;

    public PsmNormalizer(Parameters parameters) {
        this.parameters = parameters;
    }

    /**
     * Take log and normalize PSM.
     *
     * @param fileMap map of file path to list of PSM entries
     */
    public void logNormalizeData(Map<String, List<String>> fileMap) {
        for (Map.Entry<String, List<String>> entry : fileMap.entrySet()) {
            String filePath = entry.getKey();
            List<String> psmList = entry.getValue();
            Index index = parameters.indMap.get(filePath);
            List<String> normPsmList = logNormalizePsm(psmList, index);
            psmList.clear();
            psmList.addAll(normPsmList);
        }
    }

    /**
     * Normalize PSM by retention time.
     *
     * @param fileMap map of file path to list of PSM entries
     */
    public void rtNormalizeData(Map<String, List<String>> fileMap) {
        for (Map.Entry<String, List<String>> entry : fileMap.entrySet()) {
            String filePath = entry.getKey();
            List<String> psmList = entry.getValue();
            Index index = parameters.indMap.get(filePath);
            List<String> normPsmList = rtNormalizePsm(psmList, index);
            psmList.clear();
            psmList.addAll(normPsmList);
        }
    }

    // region helper methods
    private List<String> logNormalizePsm(List<String> psmList, Index index) {
        List<String> normPsmList = new ArrayList<>();
        for (String psm : psmList) {
            // calculate log2 ratio values
            String[] fields = psm.split("\t");
            double refValue = Double.parseDouble(fields[index.refIndex]) > 0 ?
                    Utils.log2(Double.parseDouble(fields[index.refIndex])) : 0;
            double[] ratioValues = new double[index.plexNum];
            for (int i = index.abnIndex; i < fields.length; i++) {
                ratioValues[i - index.abnIndex] = Double.parseDouble(fields[i]) > 0 ?
                        Utils.log2(Double.parseDouble(fields[i])) - refValue : -9999; // FIXME: Double.NaN is better
            }
            // generate normalized PSM
            StringBuilder normPsm = new StringBuilder();
            for (int i = 0; i < fields.length; i++) {
                normPsm.append(i < index.abnIndex ? (fields[i] + "\t") : (ratioValues[i - index.abnIndex] + "\t"));
            }
            normPsmList.add(normPsm.toString().trim());
        }
        return normPsmList;
    }

    private List<String> rtNormalizePsm(List<String> psmList, Index index) {
        double minRt = Double.MAX_VALUE;
        double maxRt = Double.MIN_VALUE;

        // find min and max retention time
        for (String psm : psmList) {
            double rt = extractRetentionTime(psm, index.rtIndex);
            if (rt < minRt) {
                minRt = rt;
            }
            if (rt > maxRt) {
                maxRt = rt;
            }
        }

        // define bins of retention time
        NavigableMap<Double, List<String>> binMap = Utils.createBins(minRt, maxRt, BIN_NUM);
//        NavigableMap<Double, List<String>> binMap = Utils.createBinsToBeRefactor(minRt, maxRt, BIN_NUM); // FIXME: use createBins

        // assign PSMs to bins
        for (String psm : psmList) {
            double rt = extractRetentionTime(psm, index.rtIndex);
            double binStart = binMap.floorKey(rt); // find the bin start time
            binMap.get(binStart).add(psm);
        }

        // normalize PSMs in each channel(bin)
        List<String> normPsmList = new ArrayList<>();
        for (Map.Entry<Double, List<String>> entry : binMap.entrySet()) {
            List<String> binPsmList = entry.getValue();
            normPsmList.addAll(normalizePsmList(binPsmList, index));
        }

        return normPsmList;
    }

    private double extractRetentionTime(String psm, int rtIndex) {
        String[] fields = psm.split("\t");
        return Double.parseDouble(fields[rtIndex]);
    }

    private List<String> normalizePsmList(List<String> psmList, Index index) {
        // convert to 2D array
        double[][] ratio2DValues = Utils.convertTo2DArray(psmList, index);

        // calculate median ratio for each channel
        double[] medianValues = calculateMedianByChannel(ratio2DValues, index);

        // subtract median ratio from each channel
        adjustRationsByMedian(ratio2DValues, medianValues);

        // update new ratios
        return Utils.updatePsmRatios(psmList, ratio2DValues, index);
    }

    private double[] calculateMedianByChannel(double[][] ratio2DValues, Index index) {
        double[] medianValues = new double[index.plexNum];
        for (int j = 0; j < index.plexNum; j++) {
            List<Double> channelValues = new ArrayList<>();
            for (double[] row : ratio2DValues) {
                if (row[j] != -9999) { // FIXME: !Double.isNaN(row[j]) is better
                    channelValues.add(row[j]);
                }
            }
            medianValues[j] = Utils.takeMedian(channelValues);
        }
        return medianValues;
    }

    private void adjustRationsByMedian(double[][] ratio2DValues, double[] medianValues) {
        for (double[] row : ratio2DValues) {
            for (int j = 0; j < row.length; j++) {
                if (row[j] != -9999) { // FIXME: !Double.isNaN(row[j]) is better
                    row[j] -= medianValues[j];
                }
            }
        }
    }
    // endregion
}
