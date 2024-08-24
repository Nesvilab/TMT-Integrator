package tmtintegrator.integrator;


import tmtintegrator.constants.GroupBy;
import tmtintegrator.constants.NormType;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.utils.ReportData;

import java.io.IOException;
import java.util.List;
import java.util.Map;

/**
 * Integrator to run the PSM processing pipeline
 */

public class Integrator {
    private final Parameters parameters;
    private final ReportData reportData;
    private long phaseStartTime;

    public Integrator(Parameters parameters, ReportData reportData) {
        this.parameters = parameters;
        this.reportData = reportData;
    }

    public void run(int groupByVal, int protNormVal) throws IOException {
        GroupBy groupBy = GroupBy.fromValue(groupByVal);
        NormType protNorm = NormType.fromValue(protNormVal);

        // Check if single-site, need to do multi-site first
        boolean secondProcess = false;
        if (groupBy == GroupBy.SINGLE_PHOSPHO_SITE) {
            groupBy = GroupBy.MULTI_PHOSPHO_SITE;
            secondProcess = true;
        }

        phaseStartTime = System.currentTimeMillis();

        // Load all psm.tsvs in FileMap
        PsmFileLoader loader = new PsmFileLoader(parameters);
        Map<String, List<String>> fileMap = loader.loadPsmFiles(groupBy); // key: file path, value: list of psm (0 is title)
        printTime("LoadPsms");

        // Normalize data
        PsmNormalizer normalizer = new PsmNormalizer(parameters, protNorm);
        if (parameters.abn_type == 0) {
            normalizer.logNormalizeData(fileMap);
        }
        printTime("Take log and normalize data");

        // PSM normalization
        if (parameters.psmNorm) {
            normalizer.rtNormalizeData(fileMap);
        }
        printTime("PSM normalization");

        // Group PSM
        PsmProcessor processor = new PsmProcessor(parameters, groupBy);
        processor.groupPsm(fileMap);
        printTime("Group PSM");

        // Remove outliers
        if (parameters.outlierRemoval) {
            processor.removeOutlier();
        }
        printTime("Remove outliers");

        // Collapse
        processor.collapse();
        Map<String, Map<String, double[]>> groupAbundanceMap = processor.getGroupAbundanceMap();
        printTime("Collapse");

        // Protein normalization
        if (protNorm == NormType.MC || protNorm == NormType.GN || protNorm == NormType.SL_IRS) {
            normalizer.setGroupAbundanceMap(groupAbundanceMap);
            normalizer.proteinNormalize();
            printTime("Protein normalization");
        }

        // Second process
        if (secondProcess) {
            processor.generateSingleSite();
            groupBy = GroupBy.SINGLE_PHOSPHO_SITE;
            printTime("Second process for single-site");
        }

        // Generate reports
        ReportGenerator reporter = new ReportGenerator(parameters, reportData, groupBy, protNorm, groupAbundanceMap);
        reporter.generateReports();

        printTime("Generate reports");
    }

    /**
     * Print time taken for a phase and update phaseStartTime
     * @param phase phase name
     */
    private void printTime(String phase) {
        long endTime = System.currentTimeMillis();
        System.out.println(phase + ": " + (endTime - phaseStartTime) + " ms");
        phaseStartTime = endTime;
    }
}
