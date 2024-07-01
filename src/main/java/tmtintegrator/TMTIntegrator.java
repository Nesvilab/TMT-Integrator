package tmtintegrator;

import tmtintegrator.config.ArgumentParser;
import tmtintegrator.config.ConfigLoader;
import tmtintegrator.integrator.Integrator;
import tmtintegrator.integrator.PsmPreProcessor;
import tmtintegrator.pojo.Parameters;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TMTIntegrator {

    private final static String APP_NAME = "TMT Integrator";
    private final static String APP_VERSION = "6.0.0";
    private final Parameters param;
    private List<String> proteinLi; // TODO: no usage

    public TMTIntegrator(Parameters param) {
        this.param = param;
        this.proteinLi = new ArrayList<>();
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
            // TODO: handle exception, log error
            e.printStackTrace();
        }

        long endTime = System.currentTimeMillis();
        System.out.println("Execution time: " + (endTime - startTime) + " ms");
    }

    // region helper methods
    private void run() throws IOException {
        try {
            // region check PSM tables, get genes, and build index
            long start = System.currentTimeMillis();
            PsmPreProcessor processor = new PsmPreProcessor(param);
            processor.checkPsmAndBuildIndex();
            System.out.println("Check PSM tables, get genes, and build index: " + (System.currentTimeMillis() - start) + " ms");
            // endregion

            // region preprocess PSM files
            start = System.currentTimeMillis();
            processor.updatePsmFiles();
            System.out.println("Update PSM files: " + (System.currentTimeMillis() - start) + " ms");
            System.out.println("Preprocessing finished.\n");
            // endregion

            // region integrate PSM files
            if (inValidAbundanceType()) {
                System.out.println("For raw-based abundance reports, TMT-Integrator only supports " +
                        "(1) no normlization, and (2) sample loading and internal reference scaling (SL+IRS).");
                return;
            }
            // options for groupBy and protNorm
            // FIXME: magic numbers
            int startGroupBy = (!param.geneflag) ? 0 : 1;
            int endGroupBy = (param.glycoflag) ? 5 : 4;
            if (param.minSiteProb < 0) { // not ptm data
                endGroupBy = 2;
            }
            int normalizationOptions = 2;

            if (param.protNorm >= 0) {
                processGroupBy(startGroupBy, endGroupBy);
            } else {
                processProtNorm(startGroupBy, endGroupBy, normalizationOptions);
            }
            // endregion
        } finally {
            cleanUp();
        }
    }

    private boolean inValidAbundanceType() {
        // TODO: Maybe a buggy condition here, type 1 need to be reimplemented
        // TODO: magic numbers
        return param.abn_type == 1 && (param.protNorm >= 1 && param.protNorm < 3);
    }

    // TODO: logic review required
    private void processGroupBy(int startGroupBy, int endGroupBy) throws IOException {
        Integrator integrator = new Integrator(param);
        if (param.groupBy >= 0) {
            integrator.run(param.groupBy, param.protNorm);
        } else {
            for (int i = startGroupBy; i <= endGroupBy; i++) {
                System.out.println("Start to process GroupBy: " + i);
                integrator.run(i, param.protNorm);
                System.out.println("-----------------------------------------------------------------------");
            }
        }
    }

    // FIXME 01: buggy logic, need to be reimplemented for abn_type == 1
    private void processProtNorm(int startGroupBy, int endGroupBy, int normalizationOptions) throws IOException {
        Integrator integrator = new Integrator(param);
        if (param.groupBy >= 0) {
            for (int i = 0; i <= normalizationOptions; i++) {
                if (param.abn_type == 1 && i < 1) {
                    System.out.println("Start to process protNorm=" + i);
                    integrator.run(param.groupBy, i);
                    System.out.println("-----------------------------------------------------------------------");
                } else if (param.abn_type == 0) {
                    integrator.run(param.groupBy, i);
                    System.out.println("-----------------------------------------------------------------------");
                }
            }
        } else {
            for (int i = startGroupBy; i <= endGroupBy; i++) {
                for (int j = 0; j <= normalizationOptions; j++) {
                    if ((param.abn_type == 1) && (i < 1 || i >= 3)) {
                        System.out.println("Start to process GroupBy: " + i + "_protNorm=" + j);
                        integrator.run(i, j);
                        System.out.println("-----------------------------------------------------------------------");
                    } else if (param.abn_type == 0) {
                        integrator.run(i, j);
                        System.out.println("-----------------------------------------------------------------------");
                    }
                }
            }
        }
    }

    private void cleanUp() {
        for (File file : param.fileList) {
            String tempPath = file.getAbsolutePath().substring(0, file.getAbsolutePath().lastIndexOf(".")) + ".ti";
            File tempFile = new File(tempPath);
            if (tempFile.delete()) {
                System.out.println("Deleted: " + tempPath);
            } else {
                System.err.println("Failed to delete: " + tempPath);
            }
        }
    }
    // endregion
}