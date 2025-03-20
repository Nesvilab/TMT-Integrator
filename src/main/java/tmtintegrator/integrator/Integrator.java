package tmtintegrator.integrator;


import tmtintegrator.constants.GroupBy;
import tmtintegrator.constants.NormType;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.psm.Psm;
import tmtintegrator.utils.ReportData;
import tmtintegrator.utils.Utils;

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

    public void run(int groupByVal, int protNormVal, List<Psm> psmList) throws IOException {
        phaseStartTime = System.currentTimeMillis();

        GroupBy groupBy = GroupBy.fromValue(groupByVal);
        NormType protNorm = NormType.fromValue(protNormVal);

        // Check if single-site, need to do multi-site first
        boolean secondProcess = false;
        if (groupBy == GroupBy.SINGLE_PHOSPHO_SITE) {
            groupBy = GroupBy.MULTI_PHOSPHO_SITE;
            secondProcess = true;
        }

        for (Psm psm : psmList) {
            psm.backup();
            psm.analyzeByGroup(groupBy);
        }

        // Normalize and quantify
        if (parameters.isTmt35) {
            System.out.println("Processing non-deuterium channels for TMT 35-plex:");
        }

        if (parameters.abn_type == 0) {
            Utils.logNormalizeData(psmList, parameters.isTmt35);
            printTime("Take log and normalize data");
        }

        Map<String, Map<String, double[]>> groupAbundanceMap = normAndQuantification(groupBy, protNorm, secondProcess, psmList);

        Map<String, Map<String, double[]>> dGroupAbundanceMap = null;
        if (parameters.isTmt35) {
            System.out.println("Processing deuterium channels for TMT 35-plex:");
            for (Psm psm : psmList) {
                psm.setDChannels();
            }
            dGroupAbundanceMap = normAndQuantification(groupBy, protNorm, secondProcess, psmList);
        }

        // reset groupBy for single-site
        if (secondProcess) {
            groupBy = GroupBy.SINGLE_PHOSPHO_SITE;
        }

        for (Psm psm : psmList) {
            psm.reset(); // reset all data for next groupBy
        }

        printTime("Analyzing by group " + groupBy.name());

        // Generate reports
        ReportGenerator reporter = new ReportGenerator(parameters, reportData, groupBy, protNorm, groupAbundanceMap, dGroupAbundanceMap);
        reporter.generateReports();

        printTime("Generate reports");
    }

    private Map<String, Map<String, double[]>> normAndQuantification(GroupBy groupBy, NormType protNorm, boolean secondProcess, List<Psm> psmList) {
        // Normalize data
        PsmNormalizer normalizer = new PsmNormalizer(parameters, protNorm);

        if (parameters.psmNorm) {
            normalizer.rtNormalizeData(psmList);
            printTime("PSM normalization");
        }

        // Group PSM
        PsmProcessor processor = new PsmProcessor(parameters, groupBy);
        processor.groupPsm(psmList);
        printTime("Group PSM");

        // Remove outliers
        if (parameters.outlierRemoval) {
            processor.removeOutlier();
            printTime("Remove outliers");
        }

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
            printTime("Second process for single-site");
        }

        return groupAbundanceMap;
    }

    /**
     * Print time taken for a phase and update phaseStartTime
     *
     * @param phase phase name
     */
    private void printTime(String phase) {
        long endTime = System.currentTimeMillis();
        System.out.println(phase + ": " + (endTime - phaseStartTime) + " ms");
        phaseStartTime = endTime;
    }
}
