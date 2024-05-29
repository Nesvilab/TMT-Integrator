package tmtintegrator;

import tmtintegrator.integrator.PsmPreProcessor;
import tmtintegrator.pojo.Parameters;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TMTIntegrator {
    private final Parameters param;
    private List<String> proteinLi; // TODO: no usage

    public TMTIntegrator(Parameters param) {
        this.param = param;
        this.proteinLi = new ArrayList<>();
    }

    public void run() throws IOException {
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
            // TODO: magic numbers
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

    /**
     * Helper method for getting the part of a modification to use as the index. Supports using the glycan composition
     * instead of mass if specified.
     * Note: if observed mods is empty, default to using the mass value instead
     *
     * @param inputMod             single Assigned modification string
     * @param glycanComposition    contents of the corresponding glycan composition column (only needed for glyco mode)
     * @param useGlycanComposition whether to use the glycan composition or mass as the index
     * @return mod string to use as index
     */
    public static String getAssignedModIndex(String inputMod, String glycanComposition, boolean useGlycanComposition) {
        String mod;
        if (useGlycanComposition && glycanComposition.length() > 0) {
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
    private boolean inValidAbundanceType() {
        // TODO: Maybe a buggy condition here, type 1 need to be reimplemented
        // TODO: magic numbers
        return param.abn_type == 1 && (param.protNorm >= 1 && param.protNorm < 3);
    }

    // TODO: logic review required
    private void processGroupBy(int startGroupBy, int endGroupBy) throws IOException {
        if (param.groupBy >= 0) {
            integrate.run(param, param.groupBy, param.protNorm);
        } else {
            for (int i = startGroupBy; i <= endGroupBy; i++) {
                System.out.println("Start to process GroupBy=" + i);
                integrate.run(param, i, param.protNorm);
                System.out.println("-----------------------------------------------------------------------");
            }
        }
    }

    // TODO: buggy logic, need to be reimplemented
    private void processProtNorm(int startGroupBy, int endGroupBy, int normalizationOptions) throws IOException {
        if (param.groupBy >= 0) {
            for (int i = 0; i <= normalizationOptions; i++) {
                if (param.abn_type == 1 && i < 1) {
                    System.out.println("Start to process protNorm=" + i);
                    integrate.run(param, param.groupBy, i);
                    System.out.println("-----------------------------------------------------------------------");
                } else if (param.abn_type == 0) {
                    integrate.run(param, param.groupBy, i);
                    System.out.println("-----------------------------------------------------------------------");
                }
            }
        } else {
            for (int i = startGroupBy; i <= endGroupBy; i++) {
                for (int j = 0; j <= normalizationOptions; j++) {
                    if ((param.abn_type == 1) && (i < 1 || i >= 3)) {
                        System.out.println("Start to process GroupBy=" + i + "_protNorm=" + j);
                        integrate.run(param, i, j);
                        System.out.println("-----------------------------------------------------------------------");
                    } else if (param.abn_type == 0) {
                        integrate.run(param, i, j);
                        System.out.println("-----------------------------------------------------------------------");
                    }
                }
            }
        }
    }

    private void cleanUp() {
        for (File file : param.FileLi) {
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