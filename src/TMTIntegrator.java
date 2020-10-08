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
            for (int i = 1 ; i < args.length ; i++){
                param.FileLi.add(new File(args[i]));
            }

            for(int i=0;i<param.FileLi.size();i++)
            {
                param.fNameLi.add(param.FileLi.get(i).getAbsolutePath());
            }
            Collections.sort(param.fNameLi);
            LoadParam(YamlFile); //Load parameters

            LoadFasta(); //Load fast file

            long start = System.currentTimeMillis();
            NumberFormat formatter = new DecimalFormat("#0.00000");

            GetAllGenes(param.FileLi); //Get All Genes

            for(File PsmF : param.FileLi){UpdateColumns(PsmF, param.bestPsm); }

            long end1 = System.currentTimeMillis();
            System.out.println("UpdateColumns--- " + formatter.format((end1 - start) / (1000d * 60)) + " min.");

            int gop = 4; //group options
            if(param.minSiteProb<0){//not ptm data
                gop = 2;
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
                    for(int i=0;i<=gop;i++) {//groupby
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
                for(int i=0; i<=gop; i++){//groupby
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
                        param.protExcludeAry = value.split(",");
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
                        String[] strAry = value.split(",");
                        for(int i=0;i<strAry.length;i++)
                        {
                            param.modTagLi.add(strAry[i].trim());
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
        String[] KeyAry = line.split(" ");
        String KeyStr = KeyAry[0];
        line = line.replace(">", "");
        String ValueStr = "";
        while ((line = br.readLine()) != null)
        {
            if(line.contains(">"))
            {
                KeyStr = KeyStr.replace(">", "");
                param.fastaMap.put(KeyStr, ValueStr);
                line = line.replace(">", "");

                KeyAry = line.split(" ");
                KeyStr = KeyAry[0];
                ValueStr = "";
            }
            else
            {
                ValueStr += line;
            }
        }
        br.close();
        KeyStr = KeyStr.replace(">", "");
        param.fastaMap.put(KeyStr, ValueStr);
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
                String line = "";
                while ((line = br.readLine()) != null)
                {
                    String[] strAry = line.split("\t");
                    String gene = strAry[param.ci.genecIndex].trim();
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

        for(int i = 0; i < tAry.length; i++){
            if(tAry[i].contains(param.refTag)){
                param.ci.refIndexMap.put(fName,i);
                refexit = true;
            }
            else if(param.add_Ref>=0)
            {
                param.ci.refIndexMap.put(fName,tAry.length);
                refexit = true;
            }

            if(tAry[i].equals("PeptideProphet Probability")){
                param.ci.pepProbcIndex = i;
            }
            if(tAry[i].equals("Peptide")){
                param.ci.pepcIndex = i;
            }
            if(tAry[i].equals("Assigned Modifications")){
                param.ci.assignedModcIndex = i;
            }
            if(tAry[i].equals("Phospho Site Localization")){
                param.ci.phosphoLocalcIndex = i;
            }
            if(tAry[i].equals("Protein")){
                param.ci.proteinIDcIndex = i;
            }
            if(tAry[i].equals("Gene")){
                param.ci.genecIndex = i;
            }
            if(tAry[i].equals("Is Unique")){
                param.ci.isUniquecIndex = i;
            }
            if(tAry[i].equals("Is Used")){
                param.ci.isUsedcIndex = i;
            }
            if(tAry[i].equals("Retention")){
                param.ci.rtIndex = i;
            }
            if(tAry[i].equals("Intensity")){
                param.ci.ms1IntIndex = i;
            }
            if(tAry[i].equals("Purity")){
                param.ci.purityIndex = i;
            }
            if(tAry[i].equals("Peptide")){
                param.ci.peptideIndex = i;
            }
            if(tAry[i].equals("Charge")){
                param.ci.chargeIndex = i;
            }
            if(tAry[i].equals("Observed M/Z")){
                param.ci.observedMzIndex = i;
            }
            if(tAry[i].equals("Calculated Peptide Mass") || tAry[i].equals("Peptide Mass")){
                param.ci.pepMassIndex = i;
            }
            if(tAry[i].equals("Number of Phospho Sites")){
                param.ci.numPhosphoSiteIndex = i;
            }
            if(tAry[i].equals("Mapped Genes")){
                param.ci.mapGeneIndex = i;
            }
            if(tAry[i].equals("Modified Peptide")){
                param.ci.modifiedPeptideIndex = i;
            }
            if(tAry[i].equals("Number of Enzymatic Termini")){
                param.ci.numEnzyTermi = i;
            }
        }

        param.ci.abncIndex = tAry.length-param.channelNum;

        //region check the existence of columns
        if((!refexit) && (param.add_Ref<0)){
            System.out.println("TMT-Integrator can't find the reference channel. Please check if the reference tag is correctly defined in the parameter file.");
            System.exit(1);
        }
        if(param.ci.pepcIndex<0){
            System.out.println("TMT-Integrator can't find the 'Peptide' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.pepProbcIndex<0){
            System.out.println("TMT-Integrator can't find the 'PeptideProphet Probability' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.assignedModcIndex<0){
            System.out.println("TMT-Integrator can't find the 'Assigned Modifications' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if((param.ci.phosphoLocalcIndex<0) && (param.minSiteProb>0)) {
            System.out.println("TMT-Integrator can't find the 'Phospho Site Localization' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if((param.ci.numPhosphoSiteIndex<0) && (param.minSiteProb>0)) {
            System.out.println("TMT-Integrator can't find the 'Number of Phospho Sites' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.proteinIDcIndex<0){
            System.out.println("TMT-Integrator can't find the 'Protein' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.genecIndex<0){
            System.out.println("TMT-Integrator can't find the 'Gene' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.isUniquecIndex<0){
            System.out.println("TMT-Integrator can't find the 'Is Unique' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.isUsedcIndex<0){
            System.out.println("TMT-Integrator can't find the 'Is Used' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.rtIndex<0){
            System.out.println("TMT-Integrator can't find the 'Retention' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.abncIndex<0){
            System.out.println("TMT-Integrator can't find the 'Abundance' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.ms1IntIndex<0){
            System.out.println("TMT-Integrator can't find the 'Intensity' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.purityIndex<0){
            System.out.println("TMT-Integrator can't find the 'Purity' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.peptideIndex<0){
            System.out.println("TMT-Integrator can't find the 'Peptide' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.chargeIndex<0){
            System.out.println("TMT-Integrator can't find the 'Charge' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.observedMzIndex<0){
            System.out.println("TMT-Integrator can't find the 'Observed M/Z' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.pepMassIndex<0){
            System.out.println("TMT-Integrator can't find the 'Calculated Peptide Mass' column. Please check if the column is in the psm tables.");
            System.exit(1);
        }
        if(param.ci.modifiedPeptideIndex<0){
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
        List<Double> TmtIntLi = new ArrayList<Double>();
        String line = "";
        while ((line = br.readLine()) != null)
        {
            String[] strAry = line.split("\t");

            String NewPsm = ""; //Update isUsed = false
            strAry[param.ci.isUsedcIndex] = "false";
            for(int i=0; i<strAry.length;i++){
                NewPsm+=strAry[i]+"\t";
            }

            double SumTmtInt = 0; //Sum Tmt Int
            int length = param.ci.abncIndex+param.channelNum;
            for(int i=param.ci.abncIndex; i<length; i++){
                SumTmtInt+=Double.parseDouble(strAry[i]);
            }
            TmtIntLi.add(SumTmtInt);
            NewPsm+=SumTmtInt;
            AllPsmLi.add(NewPsm);
        }
        br.close();
        Collections.sort(TmtIntLi);
        //endregion

        int TmtThresIndex = (int) Math.floor(AllPsmLi.size()*param.minPercent);
        double TmtThres = TmtIntLi.get(TmtThresIndex);

        //region Collapse Psm
        for(String Psm : AllPsmLi)
        {
            String[] strAry = Psm.split("\t");
            String fn = strAry[0].substring(strAry[0].lastIndexOf('_')+1, strAry[0].indexOf('.'));
            String PepSeq = strAry[param.ci.pepcIndex];
            String modPepSeq = strAry[param.ci.modifiedPeptideIndex];
            double purity = Double.parseDouble(strAry[param.ci.purityIndex]);
            double tmtInt = Double.parseDouble(strAry[strAry.length-1]);
            double pepProb = Double.parseDouble(strAry[param.ci.pepProbcIndex]);
            boolean isUnique = Boolean.parseBoolean(strAry[param.ci.isUniquecIndex]);
            int numPhosphoSite = param.ci.numPhosphoSiteIndex>=0? Integer.parseInt(strAry[param.ci.numPhosphoSiteIndex]) : 0;
            String assignedMod = strAry[param.ci.assignedModcIndex];
            String gene = strAry[param.ci.genecIndex];
            String proteinID = strAry[param.ci.proteinIDcIndex];
            String mapGenes = (param.ci.mapGeneIndex>=0) ? strAry[param.ci.mapGeneIndex]: "";
            double refInt = (param.add_Ref<0) ? Double.parseDouble(strAry[param.ci.refIndexMap.get(PsmF.getAbsolutePath())]) : 10000 ; //set up a random value to pass the criteria
            int ntt = Integer.parseInt(strAry[param.ci.numEnzyTermi]);

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
                for(String term : param.modTagLi){
                    if(modPepSeq.contains(term)){
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
                    if(proteinID.contains(term)){
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
                strAry[param.ci.isUsedcIndex] = "true";
                for(int i=0; i<strAry.length;i++){
                    NewPsm+=strAry[i]+"\t";
                }
                NewPsm+=gene_category; //X

                String Key = fn+"_"+strAry[param.ci.peptideIndex]+"_"+strAry[param.ci.chargeIndex]+"_"+strAry[param.ci.pepMassIndex];
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
                //PsmLi.add(Psm);
                PsmLi.add(Psm +"\t"+gene_category); //X
            }
        }
        //endregion

        if(best_psm){
            //region Select the best Psm
            for(String key : PsmMap.keySet()){
                if(!key.equals("NotUsed")){
                    List<String> PsmLi = PsmMap.get(key);
                    double maxInt = 0;
                    int bestPsmIndex = -1;
                    for(int i=0; i<PsmLi.size();i++){
                        String[] strAry = PsmLi.get(i).split("\t");
                        //double TmtInt = Double.parseDouble(strAry[strAry.length-1]);
                        double TmtInt = Double.parseDouble(strAry[strAry.length-2]); //X
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
                            strAry[param.ci.isUsedcIndex] = "true";
                        }
                        else{
                            strAry[param.ci.isUsedcIndex] = "false";
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

        //region Write PsmF
        String NewPath = PsmF.getAbsolutePath().replace(".tsv",".ti");
        //PsmF.delete();
        int PrintNum = param.ci.abncIndex + param.channelNum;
        BufferedWriter wr = new BufferedWriter(new FileWriter(NewPath));
        String[] tAry = title.split("\t");
        for(int i = 0; i < PrintNum; i++){
            wr.write(tAry[i] + "\t");
        }
        if(param.add_Ref>=0)
        {
            wr.write("Virtual_Reference_"+PsmF.getParentFile().getName()+"\t");
            //wr.write("Virtual_Reference_"+PsmF.getName()+"\t");
            //param.ci.refcIndex = PrintNum;
            param.refTag = "Virtual_Reference";
        }
        //wr.write("\tMapped Genes\tTMT Sum. Abund. \tGene Category"); //X
        wr.newLine();
        for(List<String> PsmLi : PsmMap.values()){
            for(String Psm : PsmLi){
                String[] pAry = Psm.split("\t");
                for(int i=0;i<PrintNum; i++){
                    //for(int i=0;i<pAry.length; i++){  //X
                    wr.write(pAry[i]+"\t");
                }
                if(param.add_Ref>=0)
                {
                    double rAbn = 0;
                    if(param.add_Ref==0) //summation
                    {
                        for(int i=param.ci.abncIndex; i<PrintNum; i++)
                        {
                            double value = TryParseDouble(pAry[i]);
                            rAbn += (value>=0) ? value:0;
                        }
                    }
                    else if(param.add_Ref==1) //average
                    {
                        int count=0;
                        for(int i=param.ci.abncIndex; i<PrintNum; i++)
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
                        for(int i=param.ci.abncIndex; i<PrintNum; i++)
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