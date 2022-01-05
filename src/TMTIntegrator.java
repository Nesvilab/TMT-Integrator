import java.io.*;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.lang.*;

public class TMTIntegrator
{
    private static ds_Parameters param = new ds_Parameters();

    public static void main(String[] args) throws IOException
    {
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
            LoadParam(YamlFile); //Load & check parameters
            LoadFasta(); //Load fast file

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

            long start = System.currentTimeMillis();
            NumberFormat formatter = new DecimalFormat("#0.00000");

            GetAllGenes(param.FileLi); //Get All Genes

            for(File PsmF : param.FileLi){UpdateColumns(PsmF, param.bestPsm); }

            long end1 = System.currentTimeMillis();
            System.out.println("UpdateColumns--- " + formatter.format((end1 - start) / (1000d * 60)) + " min.");

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
                if((param.abn_type==1) && (param.protNorm>1 || param.protNorm<3)){
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
                if(line.contains(":"))
                {
                    String header =  line.substring(0, line.indexOf(":")).trim();
                    String value = line.substring(line.indexOf(":")+1, line.length()).trim();
                    value = value.contains("#")? value.substring(0, value.indexOf("#")).trim() : value.trim();
                    if(header.equals("protein_database"))
                    {
                        param.fastaF= new File(value);
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
                            param.columntag = targetmass!=""?(param.modAA.replace("|","")+":"+targetmass):"";
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
                }
            }
            catch(Exception e)
            {
                System.out.println("Error at: " + line);
                System.exit(1);
            }
        }
        br.close();
    }

    private static void LoadFasta() throws IOException
    {
        BufferedReader br = new BufferedReader(new FileReader(param.fastaF));
        String line = br.readLine();
        String KeyStr = GenerateKey(line);
        String ValueStr = "";
        while ((line = br.readLine()) != null)
        {
            if(line.contains(">"))
            {
                param.fastaMap.put(KeyStr, ValueStr);

                KeyStr = GenerateKey(line);
                ValueStr = "";
            }
            else
            {
                ValueStr += line;
            }
        }
        br.close();
        param.fastaMap.put(KeyStr, ValueStr);
    }

    private static String GenerateKey(String line)
    {
        line = line.replace(">", "");
        String[] KeyAry = null;
        String KeyStr = "";

        if(line.contains("sp|") || line.contains("rev_sp|") || line.contains("tr|") || line.contains("rev_tr|")){
            KeyAry = line.split("\\|");
            if(line.contains("rev_")){
                KeyStr = "rev_"+KeyAry[1];
            }
            else{
                KeyStr = KeyAry[1];
            }
        }
        else if(line.contains("|")){
            KeyAry = line.split("\\|");
            KeyStr = KeyAry[0];
        }
        else{
            KeyAry = line.split(" ");
            KeyStr = KeyAry[0];
        }

        return KeyStr;
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
                System.exit(0);
            }
        }
    }

    private static void GetColumnIndex(String title, String fName)
    {
        String[] tAry = title.split("\t");
        boolean refexit = false;

        int RefNum = 0;
        for (String str : tAry) {
            if (str.contains(param.refTag)) {
                RefNum += 1;
            }
        }

        ds_Index indObj = param.indMap.get(fName);
        for(int i = 0; i < tAry.length; i++){
            if(param.add_Ref>=0)
            {
                indObj.refIndex = tAry.length;
                refexit = true;
            }
            else if(tAry[i].contains(param.refTag))
            {
                indObj.refIndex = i;
                refexit = true;
            }

            if(tAry[i].equals("PeptideProphet Probability")){
                indObj.pepProbcIndex = i;
            }
            if(tAry[i].equals("Peptide")){
                indObj.pepcIndex = i;
            }
            if(tAry[i].equals("Assigned Modifications")){
                indObj.assignedModcIndex = i;
            }
            if(tAry[i].equals("Phospho Site Localization") || tAry[i].equals(param.columntag)){
                indObj.ptmLocalcIndex = i;
            }
            if(tAry[i].equals("Protein ID")){
                indObj.proteinIDcIndex = i;
            }
            if(tAry[i].equals("Protein")){
                indObj.proteincIndex = i;
            }
            if(tAry[i].equals("Gene")){
                indObj.genecIndex = i;
            }
            if(tAry[i].equals("Is Unique")){
                indObj.isUniquecIndex = i;
            }
            if(tAry[i].equals("Retention")){
                indObj.rtIndex = i;
            }
            if(tAry[i].equals("Intensity")){
                indObj.ms1IntIndex = i;
            }
            if(tAry[i].equals("Purity")){
                indObj.purityIndex = i;
            }
            if(tAry[i].equals("Peptide")){
                indObj.peptideIndex = i;
            }
            if(tAry[i].equals("Charge")){
                indObj.chargeIndex = i;
            }
            if(tAry[i].equals("Observed M/Z")){
                indObj.observedMzIndex = i;
            }
            if(tAry[i].equals("Calculated Peptide Mass") || tAry[i].equals("Peptide Mass")){
                indObj.pepMassIndex = i;
            }
            if(tAry[i].equals("Mapped Genes")){
                indObj.mapGeneIndex = i;
            }
            if(tAry[i].equals("Modified Peptide")){
                indObj.modifiedPeptideIndex = i;
            }
            if(tAry[i].equals("Number of Enzymatic Termini")){
                indObj.numEnzyTermi = i;
            }
            if(tAry[i].equals("Glycan q-value")){
                indObj.glycoQvalIndex = i;
            }
        }

        indObj.abnIndex = tAry.length-param.channelNum;
        indObj.flength = tAry.length;

        //region check the existence of columns
        if((!refexit) && (param.add_Ref<0)){
            System.out.println("TMT-Integrator can't find the reference channel. Please check if the reference tag is correctly defined in the parameter file.");
            System.exit(1);
        }
        if((RefNum>1)&& (param.add_Ref<0))
        {
            System.out.println("There are more than one reference tag in the column names. Please check if the reference tag is unique among all the column names.");
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
            System.out.print(1);
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
                SumTmtInt+=Double.parseDouble(strAry[i]);
            }
            TmtIntLi.add(SumTmtInt);
            NewPsm+=SumTmtInt;
            AllPsmLi.add(NewPsm);
        }
        br.close();
        Collections.sort(TmtIntLi);

        int TmtThresIndex = (int) Math.floor(AllPsmLi.size()*param.minPercent);
        double TmtThres = TmtIntLi.get(TmtThresIndex);

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
            String gene = strAry[indObj.genecIndex];
            String proteinID = strAry[indObj.proteinIDcIndex];
            String protein = strAry[indObj.proteincIndex];
            String mapGenes = (indObj.mapGeneIndex>=0) ? strAry[indObj.mapGeneIndex]: "";
            double refInt = (param.add_Ref<0) ? Double.parseDouble(strAry[indObj.refIndex]) : 10000 ; //set up a random value to pass the criteria
            int ntt = Integer.parseInt(strAry[indObj.numEnzyTermi]);
            double psm_glycoQval = param.glycoQval >=0 ? TryParseDouble(strAry[indObj.glycoQvalIndex]) : -1;    // only parse this if glycan FDR checking requested

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

                                String mod = aMod.substring(aMod.indexOf("(")-1);
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

                            String mod = aMod.substring(aMod.indexOf("(")-1);
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

            if((purity>=param.minPurity) && (pepProb>=param.minPepProb) && (tmtInt>=TmtThres) && (gene.length()>0) &&  isAllowed && labelflag && modflag && uniqueflag && peflag && (gene_category>=param.uniqueGene) && (refInt > 0) && (ntt>=param.min_ntt)){
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
        for(int i = 0; i < (indObj.abnIndex-1); i++){
            wr.write(tAry[i] + "\t");
        }
        wr.write("Is Used (TMT-I)\t");
        for(int i=indObj.abnIndex-1; i<(PrintNum-1);i++)
        {
            wr.write(tAry[i] + "\t");
        }
        if(param.add_Ref>=0)
        {
            wr.write("Virtual_Reference_"+PsmF.getParentFile().getName()+"\t");
            param.refTag = "Virtual_Reference";
        }
        wr.newLine();
        for(List<String> PsmLi : PsmMap.values()){
            for(String Psm : PsmLi){
                String[] pAry = Psm.split("\t");
                for(int i=0;i<PrintNum; i++){
                    wr.write(pAry[i]+"\t");
                }
                if(param.add_Ref>=0)
                {
                    double rAbn = 0;
                    if(param.add_Ref==0) //summation
                    {
                        for(int i=indObj.abnIndex; i<PrintNum; i++)
                        {
                            double value = TryParseDouble(pAry[i]);
                            rAbn += (value>=0) ? value:0;
                        }
                    }
                    else if(param.add_Ref==1) //average
                    {
                        int count=0;
                        for(int i=indObj.abnIndex; i<PrintNum; i++)
                        {
                            double value = TryParseDouble(pAry[i]);
                            rAbn += (value >= 0) ? value:0;
                            count += (value > 0) ? 1:0;
                        }
                        rAbn = (count>0) ? (rAbn/count) : 0;
                    }
                    else if(param.add_Ref==2) //median
                    {
                        List<Double> rAbnLi = new ArrayList<Double>();
                        for(int i=indObj.abnIndex; i<PrintNum; i++)
                        {
                            double value = TryParseDouble(pAry[i]);
                            if(value>0)
                            {
                                rAbnLi.add(value);
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