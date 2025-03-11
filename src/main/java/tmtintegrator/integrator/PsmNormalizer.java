package tmtintegrator.integrator;

import tmtintegrator.constants.Constants;
import tmtintegrator.constants.NormType;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.psm.Psm;
import tmtintegrator.pojo.psm.PsmRecord;
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
    private Map<String, Map<String, double[]>> groupAbundanceMap; // <groupKey, <fileName, abundance>>

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
     * @param psmList list of PSM files
     */
    public void logNormalizeData(List<Psm> psmList) {
        for (Psm psm : psmList) {
            List<PsmRecord> psmRecords = psm.getPsmRecords();
            // log normalize PSM
            for (PsmRecord psmRecord : psmRecords) {
                logNormalizePsm(psmRecord);
            }
        }
    }

    /**
     * Normalize PSM by retention time.
     *
     * @param psmList list of PSM files
     */
    public void rtNormalizeData(List<Psm> psmList) {
        for (Psm psm : psmList) {
            List<PsmRecord> psmRecords = psm.getPsmRecords();
            Index index = psm.getIndex();
            // rt normalize PSM
            rtNormalizePsm(psmRecords, index);
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
    private void logNormalizePsm(PsmRecord psmRecord) {
        double refIntensity = parameters.isTmt35 ? psmRecord.getRefIntensity() : psmRecord.getMS2Intensity(); // FIXME 01: to be removed
        double logRefInt = refIntensity > 0 ? Utils.log2(refIntensity) : 0;
        List<Double> channels = psmRecord.getChannels();
        for (int i = 0; i < channels.size(); i++) {
            double intensity = channels.get(i);
            intensity = intensity > 0 ? Utils.log2(intensity) - logRefInt : Double.NaN;
            channels.set(i, intensity);
        }
    }

    private void rtNormalizePsm(List<PsmRecord> psmRecords, Index index) {
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

        NavigableMap<Double, List<PsmRecord>> binMap = Utils.createBins(minRt, maxRt, Constants.BIN_NUM);

        // assign PSMs to bins
        for (PsmRecord psmRecord : psmRecords) {
            double rt = psmRecord.getRetention();
            double binStart = binMap.floorKey(rt); // find the bin start time
            binMap.get(binStart).add(psmRecord);
        }

        // normalize PSMs in each channel(bin)
        for (Map.Entry<Double, List<PsmRecord>> entry : binMap.entrySet()) {
            List<PsmRecord> binPsmList = entry.getValue();
            normalizePsmList(binPsmList, index);
        }
    }

    private void normalizePsmList(List<PsmRecord> psmRecords, Index index) {
        // calculate median ratio for each channel
        double[] medianValues = calculateMedianByChannel(psmRecords, index);

        // subtract median ratio from each channel
        adjustRatiosByMedian(psmRecords, medianValues);
    }

    private double[] calculateMedianByChannel(List<PsmRecord> psmRecords, Index index) {
        double[] medianValues = new double[index.usedChannelNum];
        for (int j = 0; j < index.usedChannelNum; j++) {
            List<Double> channelValues = new ArrayList<>();
            for (PsmRecord psmRecord : psmRecords) {
                double intensity = psmRecord.getChannels().get(j);
                if (!Double.isNaN(intensity)) {
                    channelValues.add(intensity);
                }
            }
            medianValues[j] = Utils.takeMedian(channelValues);
        }
        return medianValues;
    }

    private void adjustRatiosByMedian(List<PsmRecord> psmRecords, double[] medianValues) {
        for (PsmRecord psmRecord : psmRecords) {
            List<Double> channels = psmRecord.getChannels();
            for (int j = 0; j < channels.size(); j++) {
                double intensity = channels.get(j);
                if (!Double.isNaN(intensity)) {
                    intensity -= medianValues[j];
                    channels.set(j, intensity);
                }
            }
        }
    }

    private Map<String, double[]> getProtMedianMap(boolean useAbsValue) {
        Map<String, double[]> protMedianMap = new HashMap<>();
        // get the median
        for (String filename : parameters.fNameLi) {
            Index index = parameters.indMap.get(filename);
            double[] medianValues = new double[index.usedChannelNum];
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
                    if (!Double.isNaN(medians[j])) {
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
                    if (!Double.isNaN(medianValues[j])) {
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
        double globalMinRefInt = Utils.calculateGlobalMinRefInt(groupAbundanceMap);
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
                int refIndex = parameters.add_Ref < 0 ? index.refIndex - index.abnIndex : -1;
                for (int j = 0; j < index.usedChannelNum; j++) {
                    if (!Double.isNaN(medianValues[j]) && j != refIndex) {
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
        Map<String, double[]> sumMap = new HashMap<>();
        // Reference channel should be excluded from the summation
        for (String filename : parameters.titleMap.keySet()) {
            Index index = parameters.indMap.get(filename);
            List<double[]> mediansList = new ArrayList<>();
            for (Map<String, double[]> fileAbundanceMap : groupAbundanceMap.values()) {
                if (fileAbundanceMap.containsKey(filename)) {
                    mediansList.add(fileAbundanceMap.get(filename));
                }
            }
            double[] sumValues = new double[index.usedChannelNum + 1];
            for (int j = 0; j < sumValues.length; j++) {
                double sum = 0;
                for (double[] medians : mediansList) {
                    if (!Double.isNaN(medians[j])) {
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
                    for (int j = 0; j < index.usedChannelNum + 1; j++) {
                        if (!Double.isNaN(medians[j])) {
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
        double globalMinRefInt = Utils.calculateGlobalMinRefInt(groupAbundanceMap);
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
                Index index = parameters.indMap.get(filename);
                double refInt = medians[medians.length - 1];
                double factor = (refInt > 0) ? avgRefInt / refInt : globalMinRefInt;
                for (int j = 0; j < index.usedChannelNum; j++) {
                    if (!Double.isNaN(medians[j])) {
                        medians[j] *= factor;
                    }
                }
            }
        }
    }
    // endregion
}
