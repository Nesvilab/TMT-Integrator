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
    public static final String name = "TMT-Integrator";
    public static final String version = "5.0.1";
    private static final NumberFormat formatter = new DecimalFormat("#0.00000");

    private static ds_Parameters param = new ds_Parameters();
    private static List<String> proteinLi = new ArrayList<String>();

    public static void main(String[] args) throws IOException
    {
        long totalStart = System.currentTimeMillis();

        System.out.printf("%s v%s%n", name, version);
        File YamlFile = args[0].contains(".yml") ? new File(args[0]) : null ;
        if(YamlFile==null) {
            System.out.println("Please input a yml file.");
            return;
        }
        else if(args[1].equalsIgnoreCase("--ValParam"))  //validate parameters
        {
            LoadParam(YamlFile);
        }
        else
        {
            long start = System.currentTimeMillis();

            LoadParam(YamlFile); //Load & check parameters

            for (int i = 1 ; i < args.length ; i++){
                param.FileLi.add(new File(args[i]));
            }

            for(int i=0;i<param.FileLi.size();i++)
            {
                param.fNameLi.add(param.FileLi.get(i).getAbsolutePath());

                ds_Index indObj = new ds_Index();
                param.indMap.put(param.FileLi.get(i).getAbsolutePath(), indObj);
            }
            Collections.sort(param.fNameLi);
            CheckPSMs(param.FileLi);//Check PSM tables
            System.out.println("Load parameters and check PSMs--- " + formatter.format((System.currentTimeMillis() - start) / (1000d * 60)) + " min.");

            start = System.currentTimeMillis();
            GetAllGenes(param.FileLi); //Get All Genes
            for(File PsmF : param.FileLi){UpdateColumns(PsmF, param.bestPsm); }
            System.out.println("UpdateColumns--- " + formatter.format((System.currentTimeMillis() - start) / (1000d * 60)) + " min.");

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

            System.out.println("Total run time--- " + formatter.format((System.currentTimeMillis() - totalStart) / (1000d * 60)) + " min.");

            System.out.println("Finish!!!");
        }

    }

    private static void CheckPSMs(List<File> FileLi) throws IOException
    {
        for(File psmF : FileLi)
        {
            BufferedReader br = new BufferedReader(new FileReader(psmF.getAbsolutePath()));
            String title = br.readLine();
            GetColumnIndex(title, psmF.getAbsolutePath());
            ds_Index indObj = param.indMap.get(psmF.getAbsolutePath());
            String line = "";
            int total=0;
            int ms1count=0;
            int genecount=0;
            while ((line = br.readLine()) != null)
            {
                String[] strAry = line.split("\t");
                ms1count += (Double.valueOf(strAry[indObj.ms1IntIndex])==0)? 1:0;
                genecount += strAry[indObj.genecIndex].trim().length()==0?1:0;
                total+=1;

                if(!proteinLi.contains(strAry[indObj.proteincIndex])){
                    proteinLi.add(strAry[indObj.proteincIndex]);
                }
                if(!param.ppMap.containsKey(strAry[indObj.proteinIDcIndex])){
                    param.ppMap.put(strAry[indObj.proteinIDcIndex], strAry[indObj.proteincIndex]);
                }
            }
            br.close();

            if(total==ms1count){
                param.ms1Int = false;
                System.out.println("Warning: All MS1 Intensities in "+ psmF.getAbsolutePath() +" are 0. TMT-Integrator will use summed MS2 reporter ion intensity instead of MS1 ion intensity as reference intensity for abundance calculation.");
                break;
            }
            if(total==genecount){
                param.geneflag = true;
                System.out.println("Warning: No gene report is generated because all gene symbols in "+ psmF.getAbsolutePath() +" are missing. The protein ID will be used as gene in other reports.");
                break;
            }
        }
    }

    private static void LoadParam(File YamlFile) throws IOException
    {
        BufferedReader br = new BufferedReader(new FileReader(YamlFile.getAbsolutePath()));
        String line = "";
        while ((line = br.readLine()) != null)
        {
            try
            {
                if (line.startsWith("#")) {
                    continue;   // skip commented lines
                }
                if(line.contains(":"))
                {
                    String header =  line.substring(0, line.indexOf(":")).trim();
                    String value = line.substring(line.indexOf(":")+1, line.length()).trim();
                    value = value.contains("#")? value.substring(0, value.indexOf("#")).trim() : value.trim();
                    if(header.equals("protein_database"))
                    {
                        //param.fastaF= new File(value);
                    }
                    else if(header.equals("output"))
                    {
                        param.reportPath = value;
                    }
                    else if(header.equals("combined_protein"))
                    {
                        param.combinedF = new File(value);
                    }
                    else if(header.equals("channel_num"))
                    {
                        param.channelNum = Integer.parseInt(value);
                    }
                    else if(header.equals("ref_tag"))
                    {
                        param.refTag = value;
                    }
                    else if(header.equals("groupby")){
                        param.groupBy = Integer.parseInt(value);
                    }
                    else if(header.equals("psm_norm"))
                    {
                        param.psmNorm = Boolean.parseBoolean(value);
                    }
                    else if(header.equals("outlier_removal"))
                    {
                        param.outlierRemoval = Boolean.parseBoolean(value);
                    }
                    else if(header.equals("prot_norm"))
                    {
                        param.protNorm = Integer.parseInt(value);
                    }
                    else if(header.equals("min_pep_prob")){
                        param.minPepProb = Float.parseFloat(value);
                    }
                    else if(header.equals("min_purity"))
                    {
                        param.minPurity = Float.parseFloat(value);
                    }
                    else if(header.equals("min_percent"))
                    {
                        param.minPercent = Float.parseFloat(value);
                    }
                    else if(header.equals("unique_pep"))
                    {
                        param.uniquePep = Boolean.parseBoolean(value);
                    }
                    else if(header.equals("unique_gene"))
                    {
                        param.uniqueGene = Integer.parseInt(value);
                    }
                    else if(header.equals("best_psm"))
                    {
                        param.bestPsm = Boolean.parseBoolean(value);
                    }
                    else if(header.equals("prot_exclude"))
                    {
                        param.protExcludeAry = (value.trim().length()==0)? "none".split(","):value.split(",");
                    }
                    else if(header.equals("allow_overlabel"))
                    {
                        param.allow_overlabel = Boolean.parseBoolean(value);
                    }
                    else if(header.equals("allow_unlabeled"))
                    {
                        param.allow_unlabeled = Boolean.parseBoolean(value);
                    }
                    else if(header.equals("mod_tag"))
                    {
                        value = (value.trim().length()==0)?"none":value;
                        if(value.equalsIgnoreCase("none"))
                        {
                            param.modTagLi.add(value);
                        }
                        else
                        {
                            String modAA="";
                            String[] strAry = value.split(",");
                            String targetmass = "";
                            for(int i=0;i<strAry.length;i++)
                            {
                                //double mass = (strAry[i].contains("(") && strAry[i].contains(")"))? Math.floor(Double.valueOf(strAry[i].substring(strAry[i].indexOf("(")+1,strAry[i].indexOf(")")))*10000)/10000:-1;
                                String mass = (strAry[i].contains("(") && strAry[i].contains(")"))?strAry[i].substring(strAry[i].indexOf("(")+1,strAry[i].indexOf(")")):"";
                                String tmp = strAry[i].contains("(")?strAry[i].substring(0,strAry[i].indexOf("("))+"("+mass+")": strAry[i];
                                param.modTagLi.add(tmp.trim());

                                if((i==0) && (mass!="")){
                                    targetmass=mass;
                                }

                                if(tmp.equalsIgnoreCase("n-glyco")||tmp.equalsIgnoreCase("o-glyco")){
                                    param.glycoflag = true;
                                }

                                if(strAry[i].contains("("))
                                {
                                    String AA=strAry[i].substring(strAry[i].indexOf('(')-1, strAry[i].indexOf('('));
                                    if(!modAA.contains(AA))
                                    {
                                        modAA +=  AA+"|";
                                    }
                                }
                            }
                            param.modAA=modAA!=""?modAA.substring(0,modAA.length()-1):"";
                            param.columntag = targetmass!=""?(param.modAA.replace("|","")+":"+targetmass.substring(0, targetmass.indexOf(".")+3)):"";
                        }
                    }
                    else if(header.equals("min_site_prob"))
                    {
                        param.minSiteProb = Float.parseFloat(value);
                    }
                    else if(header.equals("ms1_int"))
                    {
                        param.ms1Int = Boolean.parseBoolean(value);
                    }
                    else if(header.equals("top3_pep"))
                    {
                        param.top3Pep = Boolean.parseBoolean(value);
                    }
                    else if(header.equals("print_RefInt"))
                    {
                        param.print_RefInt = Boolean.parseBoolean(value);
                    }
                    else if(header.equals("add_Ref"))
                    {
                        param.add_Ref = Integer.parseInt(value);
                    }
                    else if(header.equals("max_pep_prob_thres"))
                    {
                        param.max_pep_prob_thres = Double.parseDouble(value);
                    }
                    else if(header.equals("min_ntt"))
                    {
                        param.min_ntt = Integer.parseInt(value);
                    }
                    else if(header.equals("abn_type"))
                    {
                        param.abn_type = Integer.parseInt(value);
                    }
                    else if(header.equals("aggregation_method"))
                    {
                        param.aggregation_method = Integer.parseInt(value);
                    }
                    else if(header.equals("glyco_qval"))
                    {
                        param.glycoQval = Float.parseFloat(value);
                    }
                    else if(header.equals("use_glycan_composition"))
                    {
                        param.useGlycoComposition = Boolean.parseBoolean(value);
                    }
                    else if(header.equals("prefix"))
                    {
                        param.prefix = value;
                    }
                    else if(header.equals("log2transformed"))
                    {
                        param.log2transformed = Boolean.parseBoolean(value);
                    }
                }
            }
            catch(Exception e)
            {
                System.out.println("Error at: " + line);
                e.printStackTrace();
                System.exit(1);
            }
        }
        br.close();
    }

    private static void GetAllGenes(List<File> FileLi)  throws IOException
    {
        for(int i=0; i<FileLi.size(); i++)
        {
            try
            {
                BufferedReader br = new BufferedReader(new FileReader(FileLi.get(i).getAbsolutePath()));
                String title = br.readLine();
                GetColumnIndex(title, FileLi.get(i).getAbsolutePath());
                ds_Index indObj = param.indMap.get(FileLi.get(i).getAbsolutePath());
                String line = "";
                while ((line = br.readLine()) != null)
                {
                    String[] strAry = line.split("\t");
                    String gene = strAry[indObj.genecIndex].trim();
                    if(!param.AllGeneLi.contains(gene))
                    {
                        param.AllGeneLi.add(gene);
                    }
                }
                br.close();
            }
            catch (Exception e)
            {
                System.out.println("Error at: "+FileLi.get(i).getAbsolutePath()+" "+e);
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    private static void GetColumnIndex(String title, String fName)
    {
        String[] tAry = title.split("\t");
        boolean refexit = false;

        String ref_error = "";
        int RefNum = 0;
        for(int i = 0; i < tAry.length; i++){
            String str = tAry[i].trim();
            if (str.contains(param.refTag)) {
                ref_error += (i+1) +", ";
                RefNum += 1;
            }
        }
        if(ref_error.length()>0){
            ref_error=ref_error.substring(0,ref_error.lastIndexOf(","));
        }

        ds_Index indObj = param.indMap.get(fName);
        for(int i = 0; i < tAry.length; i++){
            String str = tAry[i].trim();
            if(param.add_Ref>=0)
            {
                indObj.refIndex = tAry.length;
                refexit = true;
            }
            else if(str.contains(param.refTag))
            {
                indObj.refIndex = i;
                refexit = true;
            }

            if(str.equals("PeptideProphet Probability")){
                indObj.pepProbcIndex = i;
            }
            if(str.equals("Peptide")){
                indObj.pepcIndex = i;
            }
            if(str.equals("Assigned Modifications")){
                indObj.assignedModcIndex = i;
            }
            if(str.equals("Phospho Site Localization") || (!param.columntag.isEmpty() && str.contains(param.columntag) && !str.contains("Best Localization"))){
                indObj.ptmLocalcIndex = i;
            }
            if(str.equals("Protein ID")){
                indObj.proteinIDcIndex = i;
            }
            if(str.equals("Protein")){
                indObj.proteincIndex = i;
            }
            if(str.equals("Gene")){
                indObj.genecIndex = i;
            }
            if(str.equals("Is Unique")){
                indObj.isUniquecIndex = i;
            }
            if(str.equals("Retention")){
                indObj.rtIndex = i;
            }
            if(str.equals("Intensity")){
                indObj.ms1IntIndex = i;
            }
            if(str.equals("Purity")){
                indObj.purityIndex = i;
            }
            if(str.equals("Peptide")){
                indObj.peptideIndex = i;
            }
            if(str.equals("Charge")){
                indObj.chargeIndex = i;
            }
            if(str.equals("Observed M/Z")){
                indObj.observedMzIndex = i;
            }
            if(str.equals("Calculated Peptide Mass") || str.equals("Peptide Mass")){
                indObj.pepMassIndex = i;
            }
            if(str.equals("Mapped Genes")){
                indObj.mapGeneIndex = i;
            }
            if(str.equals("Modified Peptide")){
                indObj.modifiedPeptideIndex = i;
            }
            if(str.equals("Number of Enzymatic Termini")){
                indObj.numEnzyTermi = i;
            }
            if(str.equals("Glycan q-value")){
                indObj.glycoQvalIndex = i;
            }
            if(str.equals("Observed Modifications")){
                indObj.observedModIndex = i;
            }
            if(str.equals("Protein Start")){
                indObj.protsIndex = i;
            }
            if(str.equals("Protein End")){
                indObj.proteIndex = i;
            }
            if(str.equals("Extended Peptide")){
                indObj.extpepIndex = i;
            }
        }

        //find the number of channel
        indObj.abnIndex = tAry.length-param.channelNum;
        indObj.flength = tAry.length;

        int cnum = 0;
        for(int i=indObj.abnIndex; i<tAry.length; i++)
        {
            if(!tAry[i].trim().equalsIgnoreCase("na"))
            {
                cnum += 1;
            }
            else
            {
                if(i<indObj.refIndex){
                    indObj.refIndex -= 1;
                }
            }
        }
        indObj.totLen = (param.add_Ref<0) ? (cnum+1) : (cnum+2);
        indObj.plexNum = (param.add_Ref<0) ? (cnum) : (cnum+1);
        indObj.refIndex = (param.add_Ref<0) ? indObj.refIndex : (indObj.abnIndex+cnum);

        if(param.geneflag){
            indObj.genecIndex = indObj.proteincIndex;
        }

        //region check the existence of columns
        if((!refexit) && (param.add_Ref<0)){
            System.out.println("TMT-Integrator can't find the reference channel matching \"" + param.refTag + "\" from column " + String.join(",", tAry) + ".");
            System.out.println("Please check if the reference tag is correctly defined in the parameter file.");
            System.exit(1);
        }
        if((RefNum>1)&& (param.add_Ref<0))
        {
            System.out.println("There are more than one \""+ param.refTag + "\" in the column " + String.join(",", tAry) + ".");
            System.out.println("Repeated reference tag at column: " + ref_error + ".");
            System.out.println("Please make sure the reference tag is unique among all the column names.");
            System.exit(1);
        }
        if(indObj.pepcIndex<0){
            System.out.println("TMT-Integrator can't find the 'Peptide' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(indObj.pepProbcIndex<0){
            System.out.println("TMT-Integrator can't find the 'PeptideProphet Probability' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(indObj.assignedModcIndex<0){
            System.out.println("TMT-Integrator can't find the 'Assigned Modifications' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if((indObj.ptmLocalcIndex<0) && (param.minSiteProb>0)) {
            System.out.println("TMT-Integrator can't find the corresponding site localization column. Please check if the column: "+param.columntag+" is in the psm tables.");
            System.exit(1);
        }
        if(indObj.proteinIDcIndex<0){
            System.out.println("TMT-Integrator can't find the 'Protein' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(indObj.genecIndex<0){
            System.out.println("TMT-Integrator can't find the 'Gene' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(indObj.isUniquecIndex<0){
            System.out.println("TMT-Integrator can't find the 'Is Unique' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(indObj.rtIndex<0){
            System.out.println("TMT-Integrator can't find the 'Retention' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(indObj.ms1IntIndex<0){
            System.out.println("TMT-Integrator can't find the 'Intensity' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(indObj.purityIndex<0){
            System.out.println("TMT-Integrator can't find the 'Purity' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(indObj.peptideIndex<0){
            System.out.println("TMT-Integrator can't find the 'Peptide' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(indObj.chargeIndex<0){
            System.out.println("TMT-Integrator can't find the 'Charge' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(indObj.observedMzIndex<0){
            System.out.println("TMT-Integrator can't find the 'Observed M/Z' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(indObj.pepMassIndex<0){
            System.out.println("TMT-Integrator can't find the 'Calculated Peptide Mass' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(indObj.modifiedPeptideIndex<0){
            System.out.println("TMT-Integrator can't find the 'Modified Peptide' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }

        //endregion
    }

    private static void UpdateColumns(File PsmF, boolean best_psm)  throws IOException
    {
        List<String> AllPsmLi = new ArrayList<String>();
        TreeMap<String, List<String>> PsmMap = new TreeMap<String, List<String>>();
        PsmMap.put("NotUsed", new ArrayList<String>());

        //region Read PsmF
        BufferedReader br = new BufferedReader(new FileReader(PsmF.getAbsolutePath()));
        String title = br.readLine();
        GetColumnIndex(title, PsmF.getAbsolutePath());
        ds_Index indObj = param.indMap.get(PsmF.getAbsolutePath());
        List<Double> TmtIntLi = new ArrayList<Double>();
        String line = "";
        while ((line = br.readLine()) != null)
        {
            String[] strAry = line.split("\t");

            String NewPsm = ""; //Update isUsed = false
            for(int i=0; i<indObj.abnIndex;i++)
            {
                NewPsm+=strAry[i]+"\t";
            }
            NewPsm+="false\t";
            for(int i=indObj.abnIndex; i<strAry.length;i++){
                NewPsm+=strAry[i]+"\t";
            }

            double SumTmtInt = 0; //Sum Tmt Int
            for(int i=indObj.abnIndex; i<indObj.flength; i++){
            	try{
            		SumTmtInt+=Double.parseDouble(strAry[i]);
            	}
            	catch (Exception ex) {
                    System.err.println(ex.getMessage() + ": " + strAry[i] + " is not a number.");
                    System.err.println(String.join("\t", strAry));
                    ex.printStackTrace();
                    System.exit(1);
              }                
            }
            TmtIntLi.add(SumTmtInt);
            NewPsm+=SumTmtInt;
            AllPsmLi.add(NewPsm);
        }
        br.close();
        Collections.sort(TmtIntLi);

        int TmtThresIndex = (int) Math.floor(AllPsmLi.size()*param.minPercent);
        double TmtThres = TmtIntLi.isEmpty() ? 0 : TmtIntLi.get(TmtThresIndex);

        indObj.isUsedIndex = indObj.abnIndex;
        indObj.abnIndex+=1;
        indObj.refIndex+=1;
        indObj.flength+=1;

        List<String> NewModTagLi = new ArrayList<String>();

        //region Collapse Psm
        for(String Psm : AllPsmLi)
        {
            String[] strAry = Psm.split("\t");
            String fn = strAry[0].substring(strAry[0].lastIndexOf('_')+1, strAry[0].indexOf('.'));
            String PepSeq = strAry[indObj.pepcIndex];
            double purity = Double.parseDouble(strAry[indObj.purityIndex]);
            double tmtInt = Double.parseDouble(strAry[strAry.length-1]);
            double pepProb = Double.parseDouble(strAry[indObj.pepProbcIndex]);
            boolean isUnique = Boolean.parseBoolean(strAry[indObj.isUniquecIndex]);
            String assignedMod = strAry[indObj.assignedModcIndex];
            String protein = strAry[indObj.proteincIndex];
            String proteinID = strAry[indObj.proteinIDcIndex];
            String gene = strAry[indObj.genecIndex].length()>0?strAry[indObj.genecIndex]:proteinID;
            String mapGenes = (indObj.mapGeneIndex>=0) ? strAry[indObj.mapGeneIndex]: "";
            double refInt = (param.add_Ref<0) ? Double.parseDouble(strAry[indObj.refIndex]) : 10000 ; //set up a random value to pass the criteria
            int ntt = Integer.parseInt(strAry[indObj.numEnzyTermi]);
            double psm_glycoQval;
            try {
                psm_glycoQval = param.glycoQval >= 0 ? TryParseDouble(strAry[indObj.glycoQvalIndex]) : -1;    // only parse this if glycan FDR checking requested
            } catch (IndexOutOfBoundsException ex) {
                // non-glyco search
                System.out.println("Glycan FDR control requested but no such column found. No Glycan FDR applied.");
                param.glycoQval = -1;
                psm_glycoQval = -1;
            }
            boolean isAllowed = true;
            //region allow overlabeled
            if((!param.allow_overlabel) && (assignedMod.contains("S(229."))){
                isAllowed = false;
            }
            //endregion

            boolean labelflag=false;
            String firstAAtag = PepSeq.charAt(0)+"(229.";
            //region check if peptides with the labels
            if(param.allow_unlabeled)
            {
                labelflag=true;
            }
            else if (!param.allow_unlabeled && (assignedMod.contains("n(42.") || assignedMod.contains("n(229.") || assignedMod.contains(", 1S(229.")
                || assignedMod.indexOf("1S(229.")==0 || assignedMod.contains("N-term(42.") || assignedMod.contains("N-term(229.")
                || assignedMod.contains("N-term(144.1") || assignedMod.contains("9K(304.2") || assignedMod.contains("N-term(304.2")))
            {
                labelflag=true;
            }
            //endregion

            boolean modflag = false;
            //region  Use peptide with specified modification
            if (param.modTagLi.get(0).trim().equalsIgnoreCase("none"))
            {
                modflag = true;
            }
            else
            {
                if (assignedMod.length() == 0) {
                    // skip PSMs with nothing in assigned mods column if mods requested (shouldn't happen, but can if two mods are reported on same site and Philosopher can't separate them)
                    modflag = false;
                    continue;
                }
                for(String term : param.modTagLi){
                    term = (!term.equalsIgnoreCase("N-glyco") && !term.equalsIgnoreCase("O-glyco"))?term.substring(0, term.indexOf(".")+3): term;
                    if(term.equalsIgnoreCase("N-glyco")){
                        // if skipping PSMs that failed glycan FDR, skip if psm glycan q-val is larger than threshold provided in param
                        if (param.glycoQval >= 0 && psm_glycoQval > param.glycoQval) {
                            continue;
                        }
                        param.modAA="N";
                        String[] assignedModAry=assignedMod.split(",");
                        for(String aMod : assignedModAry)
                        {
                            if(aMod.contains("N(")){
                                double mass = Double.valueOf(aMod.substring(aMod.indexOf("N(")+2,aMod.indexOf(")")));
                                modflag=mass>=100?true:false;

                                String mod = getAssignedModIndex(aMod, strAry[indObj.observedModIndex], param.useGlycoComposition);
                                if(!NewModTagLi.contains(mod) && modflag){
                                    NewModTagLi.add(mod);
                                }
                            }
                        }
                    }
                    else if(term.equalsIgnoreCase("O-glyco")){
                        // if skipping PSMs that failed glycan FDR, skip if psm glycan q-val is larger than threshold provided in param
                        if (param.glycoQval >= 0 && psm_glycoQval > param.glycoQval) {
                            continue;
                        }
                        param.modAA="S|T";
                        String[] assignedModAry=assignedMod.split(",");
                        for(String aMod : assignedModAry)
                        {
                            boolean tflag = false;
                            if(aMod.contains("S(")){
                                double mass = Double.valueOf(aMod.substring(aMod.indexOf("S(")+2,aMod.indexOf(")")));
                                tflag=mass>=100?true:false;
                            }
                            else if(aMod.contains("T(")){
                                double mass = Double.valueOf(aMod.substring(aMod.indexOf("T(")+2,aMod.indexOf(")")));
                                tflag=mass>=100?true:false;
                            }

                            String mod = getAssignedModIndex(aMod, strAry[indObj.observedModIndex], param.useGlycoComposition);
                            if(tflag && !NewModTagLi.contains(mod)){
                                NewModTagLi.add(mod);
                            }

                            if(tflag){
                                modflag=true;
                            }
                        }
                    }
                    else if(assignedMod.contains(term)){
                        modflag= true;
                        break;
                    }
                }
            }
            //endregion

            boolean uniqueflag = false;
            //region Use unique peptide or not
            if((isUnique && param.uniquePep) || (!param.uniquePep)){
                uniqueflag = true;
            }
            //endregion

            boolean peflag = true;
            //region check excluded  proteins
            if(protein.contains("contam")){
                peflag = false;
            }
            if((param.protExcludeAry != null) && (!param.protExcludeAry[0].contains("none")))
            {
                for(String term : param.protExcludeAry){
                    if(protein.contains(term)){
                        peflag = false;
                        break;
                    }
                }
            }
            //endregion

            int gene_category = 0;
            //region gene category

            if(mapGenes!="")
            {
                String[] gAry = mapGenes.trim().split(",");
                if(gAry.length<=1)
                {
                    if((gAry[0].trim().equalsIgnoreCase(""))||(gAry[0].trim().equalsIgnoreCase(gene.trim())))
                    {
                        gene_category = 2;
                    }
                    else
                    {
                        boolean IsInUnion = false;
                        if(param.AllGeneLi.contains(gAry[0].trim()))
                        {
                            IsInUnion=true;
                        }
                        if(IsInUnion)
                        {
                            gene_category=0;
                        }
                        else
                        {
                            gene_category=1;
                        }
                    }
                }
                else
                {
                    boolean IsInUnion = false;
                    for(String str  : gAry)
                    {
                        if(!str.trim().equalsIgnoreCase(gene.trim()))
                        {
                            if(param.AllGeneLi.contains(str.trim()))
                            {
                                IsInUnion=true;
                                break;
                            }
                        }
                    }
                    if(IsInUnion)
                    {
                        gene_category=0;
                    }
                    else
                    {
                        gene_category=1;
                    }
                }
            }

            //endregion

            if((purity>=param.minPurity) && (pepProb>=param.minPepProb) && (tmtInt>=TmtThres) &&  isAllowed && labelflag && modflag && uniqueflag && peflag && (gene_category>=param.uniqueGene) && (refInt > 0) && (ntt>=param.min_ntt)){
                String NewPsm = ""; //Update isUsed = true
                strAry[indObj.isUsedIndex] = "true";
                for(int i=0; i<strAry.length;i++){
                    NewPsm+=strAry[i]+"\t";
                }
                NewPsm+=gene_category;

                String Key = fn+"_"+strAry[indObj.peptideIndex]+"_"+strAry[indObj.chargeIndex]+"_"+strAry[indObj.pepMassIndex];
                if(PsmMap.containsKey(Key)){
                    List<String> PsmLi = PsmMap.get(Key);
                    PsmLi.add(NewPsm);
                }
                else{
                    List<String> PsmLi = new ArrayList<String>();
                    PsmLi.add(NewPsm);
                    PsmMap.put(Key, PsmLi);
                }
            }
            else{
                List<String> PsmLi = PsmMap.get("NotUsed");
                PsmLi.add(Psm +"\t"+gene_category);
            }
        }
        //endregion

        param.modTagLi.addAll(NewModTagLi);

        if(best_psm){
            //region Select the best Psm
            for(String key : PsmMap.keySet()){
                if(!key.equals("NotUsed"))
                {
                    List<String> PsmLi = PsmMap.get(key);
                    double maxInt = 0;
                    int bestPsmIndex = -1;
                    for(int i=0; i<PsmLi.size();i++){
                        String[] strAry = PsmLi.get(i).split("\t");
                        double TmtInt = Double.parseDouble(strAry[strAry.length-2]);
                        if(TmtInt>maxInt){
                            maxInt = TmtInt;
                            bestPsmIndex = i;
                        }
                    }
                    //Update isUsed
                    List<String> NewPsmLi = new ArrayList<String>();
                    for(int i=0; i<PsmLi.size();i++){
                        String[] strAry = PsmLi.get(i).split("\t");
                        if(i==bestPsmIndex){
                            strAry[indObj.isUsedIndex] = "true";
                        }
                        else{
                            strAry[indObj.isUsedIndex] = "false";
                        }
                        String NewPsm = "";
                        for(int j=0; j<strAry.length;j++){
                            NewPsm+=strAry[j]+"\t";
                        }
                        NewPsmLi.add(NewPsm);
                    }
                    PsmLi.clear();
                    PsmLi.addAll(NewPsmLi);
                }
            }
            //endregion
        }

        //region Print PsmF
        String NewPath = PsmF.getAbsolutePath().replace(".tsv",".ti");
        int PrintNum = indObj.abnIndex + param.channelNum;
        BufferedWriter wr = new BufferedWriter(new FileWriter(NewPath));
        String[] tAry = title.split("\t");
        String ntStr = "";
        for(int i = 0; i < (indObj.abnIndex-1); i++){
            wr.write(tAry[i] + "\t");
            ntStr += tAry[i]+"\t";
        }
        wr.write("Is Used (TMT-I)\t");
        ntStr += "Is Used (TMT-I)\t";
        for(int i=indObj.abnIndex-1; i<(PrintNum-1);i++)
        {
            ntStr += tAry[i] + "\t";
            if(!tAry[i].trim().equalsIgnoreCase("na"))
            {
                wr.write(tAry[i] + "\t");
            }
        }
        if(param.add_Ref>=0)
        {
            ntStr += "Virtual_Reference_"+PsmF.getParentFile().getName()+"\t";
            wr.write("Virtual_Reference_"+PsmF.getParentFile().getName()+"\t");
            param.refTag = "Virtual_Reference";
        }
        wr.newLine();
        String[] ntAry = ntStr.split("\t");
        for(List<String> PsmLi : PsmMap.values()){
            for(String Psm : PsmLi){
                String[] pAry = Psm.split("\t");
                for(int i=0;i<PrintNum; i++){
                    if(!ntAry[i].trim().equalsIgnoreCase("na"))
                    {
                        wr.write(pAry[i]+"\t");
                    }
                }
                if(param.add_Ref>=0)
                {
                    double rAbn = 0;
                    if(param.add_Ref==0) //summation
                    {
                        for(int i=indObj.abnIndex; i<PrintNum; i++)
                        {
                            if(!ntAry[i].trim().equalsIgnoreCase("na"))
                            {
                                double value = TryParseDouble(pAry[i]);
                                rAbn += (value>=0) ? value:0;
                            }
                        }
                    }
                    else if(param.add_Ref==1) //average
                    {
                        int count=0;
                        for(int i=indObj.abnIndex; i<PrintNum; i++)
                        {
                            if(!ntAry[i].trim().equalsIgnoreCase("na"))
                            {
                                double value = TryParseDouble(pAry[i]);
                                rAbn += (value >= 0) ? value:0;
                                count += (value > 0) ? 1:0;
                            }
                        }
                        rAbn = (count>0) ? (rAbn/count) : 0;
                    }
                    else if(param.add_Ref==2) //median
                    {
                        List<Double> rAbnLi = new ArrayList<Double>();
                        for(int i=indObj.abnIndex; i<PrintNum; i++)
                        {
                            if(!ntAry[i].trim().equalsIgnoreCase("na"))
                            {
                                double value = TryParseDouble(pAry[i]);
                                if(value>0)
                                {
                                    rAbnLi.add(value);
                                }
                            }
                        }
                        rAbn = (rAbnLi.size()>0) ? TakeMedian(rAbnLi) : 0;
                    }
                    wr.write(rAbn + "\t");
                }
                wr.newLine();
            }
        }
        wr.close();
        //endregion
    }

    /**
     * Helper method for getting the part of a modification to use as the index. Supports using the glycan composition
     * instead of mass if specified.
     * Note: if observed mods is empty, default to using the mass value instead
     * @param inputMod single Assigned modification string
     * @param observedMods contents of the corresponding observed mods column (only needed for glyco mode)
     * @param useGlycanComposition whether to use the glycan composition or mass as the index
     * @return mod string to use as index
     */
    public static String getAssignedModIndex(String inputMod, String observedMods, boolean useGlycanComposition) {
        String mod;
        if (useGlycanComposition && observedMods.length() > 0) {
            // if using composition for index, read it from observed mods column. Still get AA site from assigned mods
            mod = inputMod.substring(inputMod.indexOf("(")-1, inputMod.indexOf("("));
            mod = String.format("%s(%s)", mod, observedMods);
        } else {
            // read mass from assigned mod, as would do for any other mod
            mod = inputMod.substring(inputMod.indexOf("(")-1);
        }
        return mod;
    }

    private static double TryParseDouble(String str)
    {
        double value = 0;
        try
        {
            value = Double.parseDouble(str);
        }
        catch (NumberFormatException e)
        {
            value = -1;
        }
        return value;
    }

    private static double TakeMedian(List<Double> tmpLi)
    {
        Collections.sort(tmpLi);
        double mValue;
        if(tmpLi.size() == 0)
        {
            mValue = -9999;
        }
        else if(tmpLi.size() == 1)
        {
            mValue = tmpLi.get(0);
        }
        else if(tmpLi.size() == 2)
        {
            mValue = (tmpLi.get(0)+tmpLi.get(1))/2;
        }
        else
        {
            int mod = tmpLi.size()%2;
            int i1= tmpLi.size()/2;
            if(mod==0) //even
            {
                int i2 = i1-1;
                mValue = (tmpLi.get(i1) + tmpLi.get(i2))/2 ;
            }
            else //odd
            {
                mValue = tmpLi.get(i1);
            }
        }

        return mValue;
    }
}