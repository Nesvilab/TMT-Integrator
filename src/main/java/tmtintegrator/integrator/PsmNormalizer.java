package tmtintegrator.integrator;

import tmtintegrator.constants.Constants;
import tmtintegrator.constants.NormType;
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
    private final NormType normType; // protein normalization type
    private Map<String, Map<String, double[]>> groupAbundanceMap; // <groupKey(proteinId), <fileName, abundance>>

    public PsmNormalizer(Parameters parameters, NormType normType) {
        this.parameters = parameters;
        this.normType = normType;
    }

    public void setGroupAbundanceMap(Map<String, Map<String, double[]>> groupAbundanceMap) {
        this.groupAbundanceMap = groupAbundanceMap;
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

    /**
     * Normalize protein abundance.
     */
    public void proteinNormalize() {
        if (normType == NormType.MC || normType == NormType.GN) {
            // preform median centering first
            Map<String, double[]> protMedianMap = getProtMedianMap(false);
            double globalMedian = Utils.calculateGlobalMedian(protMedianMap);
            subtractProt(protMedianMap, -1, -1, true); // subtract protein ratios

            if (normType == NormType.GN) {
                // for GN, perform variance scaling
                Map<String, double[]> absProtMedianMap = getProtMedianMap(true);
                double globalAbsMedian = Utils.calculateGlobalMedian(absProtMedianMap);
                subtractProt(absProtMedianMap, globalMedian, globalAbsMedian, false); // subtract protein ratios
            }
        } else if (normType == NormType.SL_IRS) {
            // Convert ratio to abundance if not raw-based
            if (parameters.abn_type == 0) {
                ratioToAbundance();
            }
            slNormalize(); // sample loading
            irsNormalize(); // internal reference scaling
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
        NavigableMap<Double, List<String>> binMap = Utils.createBins(minRt, maxRt, Constants.BIN_NUM);

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

    private Map<String, double[]> getProtMedianMap(boolean useAbsValue) {
        Map<String, double[]> protMedianMap = new TreeMap<>(); // <filename, medianValues> TODO: HashMap?
        // get the median
        for (String filename : parameters.titleMap.keySet()) {
            Index index = parameters.indMap.get(filename);
            double[] medianValues = new double[index.plexNum];
            List<double[]> mediansList = new ArrayList<>();
            // take all abundance values
            for (Map<String, double[]> fileAbundanceMap : groupAbundanceMap.values()) {
                if (fileAbundanceMap.containsKey(filename)) {
                    mediansList.add(fileAbundanceMap.get(filename));
                }
            }
            // calculate median
            for (int j = 0; j < medianValues.length; j++) {
                List<Double> channelValues = new ArrayList<>();
                for (double[] medians : mediansList) {
                    if (medians[j] != -9999) { // FIXME: !Double.isNaN(medians[j]) is better
                        channelValues.add(useAbsValue ? Math.abs(medians[j]) : medians[j]);
                    }
                }
                medianValues[j] = Utils.takeMedian(channelValues);
            }
            protMedianMap.put(filename, medianValues);
        }
        return protMedianMap;
    }

    private void subtractProt(Map<String, double[]> protMedianMap, double globalMedian, double globalAbsMedian, boolean forMC) {
        for (Map.Entry<String, Map<String, double[]>> entry : groupAbundanceMap.entrySet()) {
            Map<String, double[]> fileAbundanceMap = entry.getValue();
            for (String filename : fileAbundanceMap.keySet()) {
                double[] protMedianValues = protMedianMap.get(filename);
                double[] medianValues = fileAbundanceMap.get(filename);
                for (int j = 0; j < protMedianValues.length; j++) {
                    if (medianValues[j] != -9999) { // FIXME: !Double.isNaN(medianValues[j]) is better
                        if (forMC) {
                            medianValues[j] -= protMedianValues[j];
                        } else { // for GN variance scaling
                            // medianValues[j] = (medianValues[j] / protMedianValues[j]) * globalAbsMedian + globalMedian;
                            medianValues[j] = (medianValues[j] / protMedianValues[j]) * globalAbsMedian;
                        }
                    }
                }
            }
        }
    }

    private void ratioToAbundance() {
        double globalMinRefInt = Utils.calculateGlobalMinRefInt(groupAbundanceMap, parameters);
        for (Map<String, double[]> fileAbundanceMap : groupAbundanceMap.values()) {
            double avgAbundance = Utils.calculateAvgAbundance(fileAbundanceMap, globalMinRefInt, parameters);
            convertRatioToAbundance(fileAbundanceMap, avgAbundance);
        }
    }

    private void convertRatioToAbundance(Map<String, double[]> fileAbundanceMap, double avgAbundance) {
        for (String filename : parameters.fNameLi) {
            double[] medianValues = fileAbundanceMap.get(filename);
            if (medianValues != null) {
                Index index = parameters.indMap.get(filename);
                int refIndex = index.refIndex - index.abnIndex;
                for (int j = 0; j < index.plexNum; j++) {
                    if (medianValues[j] != -9999 && j != refIndex) { // FIXME: !Double.isNaN(medianValues[j]) is better
                        medianValues[j] = Utils.log2(Utils.pow2(medianValues[j]) * avgAbundance);
                    }
                }
            }
        }
    }

    /**
     * Sample loading normalization.
     */
    private void slNormalize() {
        // 1. get the summation of all channels
        Map<String, double[]> sumMap = getSumMap();
        // 2. compute sum average
        double sumAvg = calculateSumAvg(sumMap);
        // 3. adjust intensity using factors
        adjustIntensity(sumMap, sumAvg);
    }

    private Map<String, double[]> getSumMap() {
        Map<String, double[]> sumMap = new TreeMap<>(); // <filename, sumValues> TODO: HashMap?
        // Reference channel should be excluded from the summation
        for (String filename : parameters.titleMap.keySet()) {
            Index index = parameters.indMap.get(filename);
            List<double[]> mediansList = new ArrayList<>();
            for (Map<String, double[]> fileAbundanceMap : groupAbundanceMap.values()) {
                if (fileAbundanceMap.containsKey(filename)) {
                    mediansList.add(fileAbundanceMap.get(filename));
                }
            }
            double[] sumValues = new double[index.totLen];
            for (int j = 0; j < sumValues.length; j++) {
                double sum = 0;
                for (double[] medians : mediansList) {
                    if (medians[j] != -9999) { // FIXME: !Double.isNaN(medians[j]) is better
                        sum += medians[j];
                    }
                }
                sumValues[j] = sum;
            }
            sumMap.put(filename, sumValues);
        }
        return sumMap;
    }

    private double calculateSumAvg(Map<String, double[]> sumMap) {
        double sumAvg = 0;
        int count = 0;
        for (double[] sumValues : sumMap.values()) {
            for (double sum : sumValues) {
                if (sum > 0) {
                    sumAvg += sum;
                    count++;
                }
            }
        }
        return sumAvg / count;
    }

    private void adjustIntensity(Map<String, double[]> sumMap, double sumAvg) {
        for (String filename : parameters.titleMap.keySet()) {
            Index index = parameters.indMap.get(filename);
            double[] sumValues = sumMap.get(filename);
            for (Map<String, double[]> fileAbundanceMap : groupAbundanceMap.values()) {
                if (fileAbundanceMap.containsKey(filename)) {
                    double[] medians = fileAbundanceMap.get(filename);
                    for (int j = 0; j < index.totLen; j++) {
                        if (medians[j] != -9999) { // FIXME: !Double.isNaN(medians[j]) is better
                            medians[j] *= (sumAvg / sumValues[j]);
                        }
                    }
                }
            }
        }
    }

    /**
     * Internal reference scaling normalization.
     */
    private void irsNormalize() {
        double globalMinRefInt = Utils.calculateGlobalMinRefInt(groupAbundanceMap, parameters);
        for (Map<String, double[]> fileAbundanceMap : groupAbundanceMap.values()) {
            double avgRefInt = calculateAvgRefInt(fileAbundanceMap);
            adjustIntensityForIRS(fileAbundanceMap, avgRefInt, globalMinRefInt);
        }
    }

    private double calculateAvgRefInt(Map<String, double[]> fileAbundanceMap) {
        double sumRefInt = 0;
        int count = 0;
        for (String filename : parameters.titleMap.keySet()) {
            if (fileAbundanceMap.containsKey(filename)) {
                double[] medians = fileAbundanceMap.get(filename);
                double refInt = medians[medians.length - 1];
                if (refInt > 0) {
                    sumRefInt += Math.log(refInt);
                    count++;
                }
            }
        }
        return Math.exp(sumRefInt / count);
    }

    private void adjustIntensityForIRS(Map<String, double[]> fileAbundanceMap, double avgRefInt, double globalMinRefInt) {
        for (String filename : parameters.titleMap.keySet()) {
            if (fileAbundanceMap.containsKey(filename)) {
                double[] medians = fileAbundanceMap.get(filename);
                double refInt = medians[medians.length - 1];
                double factor = (refInt > 0) ? avgRefInt / refInt : globalMinRefInt;
                for (int j = 0; j < parameters.channelNum; j++) {
                    if (medians[j] != -9999) { // FIXME: !Double.isNaN(medians[j]) is better
                        medians[j] *= factor;
                    }
                }
            }
        }
    }
    // endregion
}
