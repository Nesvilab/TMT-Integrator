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

import tmtintegrator.constants.GroupBy;
import tmtintegrator.constants.NormType;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.psm.Psm;
import tmtintegrator.utils.ReportData;

import java.io.IOException;
import java.util.ArrayList;
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

        // Process each plex
        List<Map<String, Map<String, double[]>>> plexAbundanceMaps = new ArrayList<>();
        for (int i = 0; i < parameters.subplexes.size(); i++) {
            System.out.println("Processing plex " + (i + 1) + " (reference channel: " + parameters.subplexes.get(i).refTag + "):");
            for (Psm psm : psmList) {
                psm.setActivePlex(i);
            }
            plexAbundanceMaps.add(normAndQuantification(groupBy, protNorm, secondProcess, psmList));
        }

        // reset groupBy for single-site
        if (secondProcess) {
            groupBy = GroupBy.SINGLE_PHOSPHO_SITE;
        }

        for (Psm psm : psmList) {
            psm.reset();
        }

        printTime("Analyzing by group " + groupBy.name());

        // Generate reports
        ReportGenerator reporter = new ReportGenerator(parameters, reportData, groupBy, protNorm, plexAbundanceMaps);
        reporter.generateReports();

        printTime("Generate reports");
    }

    private Map<String, Map<String, double[]>> normAndQuantification(GroupBy groupBy, NormType protNorm, boolean secondProcess, List<Psm> psmList) {
        PsmProcessor processor = new PsmProcessor(parameters, groupBy);
        processor.groupPsm(psmList);
        printTime("Group PSM");

        if (parameters.outlierRemoval) {
            processor.removeOutlier();
            printTime("Remove outliers");
        }

        processor.collapse();
        Map<String, Map<String, double[]>> groupAbundanceMap = processor.getGroupAbundanceMap();
        printTime("Collapse");

        if (protNorm == NormType.MC || protNorm == NormType.GN || protNorm == NormType.SL_IRS) {
            PsmNormalizer normalizer = new PsmNormalizer(parameters, protNorm);
            normalizer.setGroupAbundanceMap(groupAbundanceMap);
            normalizer.proteinNormalize();
            printTime("Protein normalization");
        }

        if (secondProcess) {
            processor.generateSingleSite();
            printTime("Second process for single-site");
        }

        return groupAbundanceMap;
    }

    private void printTime(String phase) {
        long endTime = System.currentTimeMillis();
        System.out.println(phase + ": " + (endTime - phaseStartTime) + " ms");
        phaseStartTime = endTime;
    }
}
