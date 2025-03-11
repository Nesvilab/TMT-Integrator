package tmtintegrator;

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
    private final static String APP_VERSION = "6.1.0";
    private final Parameters param;
    private final ReportData reportData;

    public TMTIntegrator(Parameters param) {
        this.param = param;
        this.reportData = new ReportData();
    }

    public static void main(String[] args) throws IOException {

        System.out.println(APP_NAME + " " + APP_VERSION);

        long startTime = System.currentTimeMillis();

        // Process command line arguments
        ArgumentParser argumentParser = new ArgumentParser(args);

        // Load parameters from the YAML file
        long startLoadTime = System.currentTimeMillis();
        ConfigLoader configLoader = new ConfigLoader();
        configLoader.loadParameters(argumentParser.getYamlFile());
        if (argumentParser.isValidateOnly()) {
            System.out.println("Validating parameters only");
            return;
        }

        // Load input files
        configLoader.loadFileList(argumentParser.getInputFiles());
        long endLoadTime = System.currentTimeMillis();
        System.out.println("Parameter Loading: " + (endLoadTime - startLoadTime) + " ms");

        try {
            TMTIntegrator integrator = new TMTIntegrator(configLoader.getParameters());
            integrator.run();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        long endTime = System.currentTimeMillis();
        System.out.println("Execution time: " + (endTime - startTime) + " ms");
    }

    // region helper methods
    private void run() throws IOException {
        // region check PSM tables, get genes, and build index
        long start = System.currentTimeMillis();
        PsmPreProcessor processor = new PsmPreProcessor(param, reportData);
        processor.checkPsmAndBuildIndex();
        processor.collectProteinInfo();
        System.out.println("Check PSM tables, get genes, and build index: " + (System.currentTimeMillis() - start) + " ms");
        // endregion

        // region preprocess PSM files
        start = System.currentTimeMillis();
        List<Psm> psmList = processor.parseAndFilterPsmFiles();
        System.out.println("Parsing and Filtering: " + (System.currentTimeMillis() - start) + " ms");
        System.out.println("Preprocessing finished.\n");
        // endregion

        // region integrate PSM files
        if (inValidAbundanceType()) {
            System.out.println("For raw-based abundance reports, TMT-Integrator only supports " +
                    "(1) no normlization, and (2) sample loading and internal reference scaling (SL+IRS).");
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

    // TODO: logic review required
    private void processGroupBy(int startGroupBy, int endGroupBy, List<Psm> psmList) throws IOException {
        Integrator integrator = new Integrator(param, reportData);
        if (param.groupBy >= 0) {
            integrator.run(param.groupBy, param.protNorm, psmList);
        } else {
            for (int i = startGroupBy; i <= endGroupBy; i++) {
                System.out.println("Start to process GroupBy: " + i);
                integrator.run(i, param.protNorm, psmList);
                System.out.println("-----------------------------------------------------------------------");
            }
        }
    }

    // FIXME 01: buggy logic, need to be reimplemented for abn_type == 1
    private void processProtNorm(int startGroupBy, int endGroupBy, int normalizationOptions, List<Psm> psmList) throws IOException {
        Integrator integrator = new Integrator(param, reportData);
        if (param.groupBy >= 0) {
            for (int i = 0; i <= normalizationOptions; i++) {
                if (param.abn_type == 1 && i < 1) {
                    System.out.println("Start to process protNorm=" + i);
                    integrator.run(param.groupBy, i, psmList);
                    System.out.println("-----------------------------------------------------------------------");
                } else if (param.abn_type == 0) {
                    integrator.run(param.groupBy, i, psmList);
                    System.out.println("-----------------------------------------------------------------------");
                }
            }
        } else {
            for (int i = startGroupBy; i <= endGroupBy; i++) {
                for (int j = 0; j <= normalizationOptions; j++) {
                    if ((param.abn_type == 1) && (i < 1 || i >= 3)) {
                        System.out.println("Start to process GroupBy: " + i + "_protNorm=" + j);
                        integrator.run(i, j, psmList);
                        System.out.println("-----------------------------------------------------------------------");
                    } else if (param.abn_type == 0) {
                        integrator.run(i, j, psmList);
                        System.out.println("-----------------------------------------------------------------------");
                    }
                }
            }
        }
    }
}