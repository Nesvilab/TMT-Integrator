package tmtintegrator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;

public class TMTIntegrator
{
    private static final NumberFormat formatter = new DecimalFormat("#0.00000");

    private ds_Parameters param;
    private List<String> proteinLi; // TODO: no usage

    public TMTIntegrator(ds_Parameters param) {
        this.param = param;
        this.proteinLi = new ArrayList<>();
    }

    public void run() throws IOException
    {
            // region check PSM tables, get genes, and build index
            long start = System.currentTimeMillis();
            PsmProcessor psmProcessor = new PsmProcessor(param);
            psmProcessor.checkPsmAndBuildIndex();
            System.out.println("Check PSM tables, get genes, and build index: " + (System.currentTimeMillis() - start) + " ms");
            // endregion

            // region preprocess PSM files
            start = System.currentTimeMillis();
            psmProcessor.updatePsmFiles();
            System.out.println("Update PSM files: " + (System.currentTimeMillis() - start) + " ms");
            System.out.println("Preprocessing finished.\n");
            // endregion

            int start_gop = (!param.geneflag) ? 0: 1;
            int end_gop = (param.glycoflag)? 5: 4; //group options
            if(param.minSiteProb<0){//not ptm data
                end_gop = 2;
            }
            int nop = 2; //Normalization options

            if ((param.groupBy>=0) && (param.protNorm>=0)){
                if((param.abn_type==1) && (param.protNorm>1 || param.protNorm<3)){
                    System.out.println("For raw-based abundance reports, TMT-Integrator only supports " +
                            "(1) no normlization, and (2) sample loading and internal reference scaling (SL+IRS).");
                }
                else{
                    integrate.run(param, param.groupBy, param.protNorm);
                }
            }
            else if((param.groupBy<0) && (param.protNorm>=0)){
                if((param.abn_type==1) && (param.protNorm>=1 && param.protNorm<3)){
                    System.out.println("For raw-based abundance reports, TMT-Integrator only supports " +
                            "(1) no normlization, and (2) sample loading and internal reference scaling (SL+IRS).");
                }
                else{
                    for(int i=start_gop;i<=end_gop;i++) {//groupby
                        System.out.println("Start to process GroupBy="+i);
                        integrate.run(param, i, param.protNorm);
                        System.out.println("-----------------------------------------------------------------------");
                    }
                }
            }
            else if((param.groupBy>=0) && (param.protNorm<0)){
                for(int i=0;i<=nop;i++) {//protNorm
                    if((param.abn_type==1) && (i<1 || i>=3)){
                        System.out.println("Start to process protNorm="+i);
                        integrate.run(param, param.groupBy, i);
                        System.out.println("-----------------------------------------------------------------------");
                    }
                    else if (param.abn_type==0){
                        integrate.run(param, param.groupBy, i);
                        System.out.println("-----------------------------------------------------------------------");
                    }
                }
            }
            else if((param.groupBy<0) && (param.protNorm<0)){
                for(int i=start_gop; i<=end_gop; i++){//groupby
                    for(int j=0; j<=nop; j++) {//protNorm
                        if((param.abn_type==1) && (i<1 || i>=3)){
                            System.out.println("Start to process GroupBy="+i+"_protNorm="+j);
                            integrate.run(param, i, j);
                            System.out.println("-----------------------------------------------------------------------");
                        }
                        else if (param.abn_type==0){
                            integrate.run(param, i, j);
                            System.out.println("-----------------------------------------------------------------------");
                        }
                    }
                }
            }

            //region delete ti files

            for(File f : param.FileLi)
            {
                String NewPath = f.getAbsolutePath().substring(0, f.getAbsolutePath().lastIndexOf("."))+".ti";
                File file = new File(NewPath);
                file.delete();
            }

            //endregion
    }

    /**
     * Helper method for getting the part of a modification to use as the index. Supports using the glycan composition
     * instead of mass if specified.
     * Note: if observed mods is empty, default to using the mass value instead
     * @param inputMod single Assigned modification string
     * @param glycanComposition contents of the corresponding glycan composition column (only needed for glyco mode)
     * @param useGlycanComposition whether to use the glycan composition or mass as the index
     * @return mod string to use as index
     */
    public static String getAssignedModIndex(String inputMod, String glycanComposition, boolean useGlycanComposition) {
        String mod;
        if (useGlycanComposition && glycanComposition.length() > 0) {
            // if using composition for index, read it from glycan composition column. Still get AA site from assigned mods
            mod = inputMod.substring(inputMod.indexOf("(")-1, inputMod.indexOf("("));
            mod = String.format("%s(%s)", mod, glycanComposition);
        } else {
            // read mass from assigned mod, as would do for any other mod
            mod = inputMod.substring(inputMod.indexOf("(")-1);
        }
        return mod;
    }
}