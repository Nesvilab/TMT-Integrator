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

package tmtintegrator;

import static tmtintegrator.utils.Utils.myPrint;

import java.io.IOException;
import java.util.List;

import tmtintegrator.config.ArgumentParser;
import tmtintegrator.config.ConfigLoader;
import tmtintegrator.integrator.Integrator;
import tmtintegrator.integrator.PsmPreProcessor;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.psm.Psm;
import tmtintegrator.utils.ReportData;

public class TMTIntegrator {

    private final static String APP_NAME = "TMT Integrator";
    private final static String APP_VERSION = "6.1.3";
    private final Parameters param;
    private final ReportData reportData;

    public TMTIntegrator(Parameters param) {
        this.param = param;
        this.reportData = new ReportData();
    }

    public static void main(String[] args) throws IOException {

        myPrint(APP_NAME + " " + APP_VERSION, "INFO");

        // Process command line arguments
        ArgumentParser argumentParser = new ArgumentParser(args);

        // Load parameters from the YAML file
        ConfigLoader configLoader = new ConfigLoader();
        configLoader.loadParameters(argumentParser.getYamlFile());
        if (argumentParser.isValidateOnly()) {
            myPrint("Validating parameters only", "INFO");
            return;
        }

        // Load input files
        configLoader.loadFileList(argumentParser.getInputFiles());

        try {
            TMTIntegrator integrator = new TMTIntegrator(configLoader.getParameters());
            integrator.run();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    // region helper methods
    private void run() throws IOException {
        param.printAllParameters();

        // region check PSM tables, get genes, and build index
        PsmPreProcessor processor = new PsmPreProcessor(param, reportData);
        processor.checkPsmAndBuildIndex();
        processor.collectProteinInfo();
        // endregion

        // region preprocess PSM files
        List<Psm> psmList = processor.preprocessPsm();
        myPrint("Preprocessing finished.", "INFO");
        // endregion

        // region integrate PSM files
        if (inValidAbundanceType()) {
            myPrint("For raw-based abundance reports, TMT-Integrator only supports " +
                    "(1) no normlization, and (2) sample loading and internal reference scaling (SL+IRS).", "WARN");
            return;
        }
        // options for groupBy and protNorm
        int startGroupBy = (!param.geneflag) ? 0 : 1;
        int endGroupBy = (param.glycoflag) ? 5 : 4;
        if (param.minSiteProb < 0) { // not ptm data
            endGroupBy = 2;
        }
        int normalizationOptions = 2;

        if (param.protNorm >= 0) {
            processGroupBy(startGroupBy, endGroupBy, psmList);
        } else {
            processProtNorm(startGroupBy, endGroupBy, normalizationOptions, psmList);
        }
    }

    private boolean inValidAbundanceType() {
        // TODO: Maybe a buggy condition here, type 1 need to be reimplemented
        return param.abn_type == 1 && (param.protNorm >= 1 && param.protNorm < 3);
    }

    private void processGroupBy(int startGroupBy, int endGroupBy, List<Psm> psmList) throws IOException {
        Integrator integrator = new Integrator(param, reportData);
        if (param.groupBy >= 0) {
            integrator.run(param.groupBy, param.protNorm, psmList);
        } else {
            for (int i = startGroupBy; i <= endGroupBy; i++) {
                myPrint("Start to process GroupBy: " + i, "INFO");
                integrator.run(i, param.protNorm, psmList);
            }

            if (endGroupBy == 2 && param.modAA.isEmpty() && "none".equals(param.modTagSet.stream().findFirst().orElse(null))) {
                // currently the modified peptide reports are only generated when on mod tags are specified because of the integral PSM filtering under multiple criteria
                int j = 6;
                myPrint("Start to process GroupBy: " + j, "INFO");
                integrator.run(j, param.protNorm, psmList);
            }
        }
    }

    // FIXME 01: buggy logic, need to be reimplemented for abn_type == 1
    private void processProtNorm(int startGroupBy, int endGroupBy, int normalizationOptions, List<Psm> psmList) throws IOException {
        Integrator integrator = new Integrator(param, reportData);
        if (param.groupBy >= 0) {
            for (int i = 0; i <= normalizationOptions; i++) {
                if (param.abn_type == 1 && i < 1) {
                    myPrint("Start to process protNorm=" + i, "INFO");
                    integrator.run(param.groupBy, i, psmList);
                } else if (param.abn_type == 0) {
                    integrator.run(param.groupBy, i, psmList);
                }
            }
        } else {
            for (int i = startGroupBy; i <= endGroupBy; i++) {
                for (int j = 0; j <= normalizationOptions; j++) {
                    if ((param.abn_type == 1) && (i < 1 || i >= 3)) {
                        myPrint("Start to process GroupBy: " + i + "_protNorm=" + j, "INFO");
                        integrator.run(i, j, psmList);
                    } else if (param.abn_type == 0) {
                        integrator.run(i, j, psmList);
                    }
                }
            }
        }
    }
}
