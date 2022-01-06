
import java.io.*;
import java.lang.reflect.Array;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

public class integrate
{
    private static ds_Parameters param = new ds_Parameters();
    private static int groupBy = -1;
    private static int protNorm = -1;
    private static int TotLen = -1;
    private static int PlexNum = -1;
    private static TreeMap<String, TreeMap<String, ds_PsmInfo>> gPsmMap = new TreeMap<String, TreeMap<String, ds_PsmInfo>>();
    private static TreeMap<String, TreeMap<String, double[]>> gAbnMap = new TreeMap<String, TreeMap<String, double[]>>();

    public static void run(ds_Parameters parameter, int gb, int pn) throws IOException
    {
        param = parameter;
        groupBy = gb;
        protNorm = pn;
        TotLen = (param.add_Ref<0) ? (param.channelNum+1) : (param.channelNum+2);
        PlexNum = (param.add_Ref<0) ? (param.channelNum) : (param.channelNum+1);

        Boolean SecondProcess = false;
        if(gb==4){
            groupBy = 3;
            SecondProcess = true;
        }

        long end1 = System.currentTimeMillis();
        NumberFormat formatter = new DecimalFormat("#0.00000");

        //Load all psm.tsvs in FileMap
        TreeMap<String, List<String>> FileMap = LoadPsms(param.FileLi); //Store the psms, key: fName; value: psm list (0 is the title)

        long  end2 = System.currentTimeMillis();
        System.out.println("LoadPsms--- " + formatter.format((end2 - end1) / (1000d * 60)) + " min.");

        if(param.abn_type==0)
        {
            for(String fName : FileMap.keySet()){//Take log and normalize
                List<String> PsmLi = FileMap.get(fName);
                ds_Index indObj = param.indMap.get(fName);
                List<String> NewPsmLi = LogNor(PsmLi, indObj);
                PsmLi.clear();
                PsmLi.addAll(NewPsmLi);
            }
        }

        long  end3 = System.currentTimeMillis();
        System.out.println("Take log and normalize--- " + formatter.format((end3 - end2) / (1000d * 60)) + " min.");

        if(param.psmNorm){ //PSM normalization
            for(String fName : FileMap.keySet()){
                List<String> PsmLi = FileMap.get(fName);
                ds_Index indObj = param.indMap.get(fName);
                List<String> NewPsmLi = RtNorm(PsmLi, indObj);
                PsmLi.clear();
                PsmLi.addAll(NewPsmLi);
            }
        }

        long  end4 = System.currentTimeMillis();
        System.out.println("PSM normalization--- " + formatter.format((end4 - end3) / (1000d * 60)) + " min.");

        GroupPsm(FileMap);
        FileMap.clear();

        if(param.outlierRemoval){RemoveOutlier();}

        long  end5 = System.currentTimeMillis();
        System.out.println("outlierRemoval--- " + formatter.format((end5 - end4) / (1000d * 60)) + " min.");

        Collapse();

        long  end6 = System.currentTimeMillis();
        System.out.println("Collapse--- " + formatter.format((end6 - end5) / (1000d * 60)) + " min.");

        if(protNorm > 0){ProtNor();}

        long  end7 = System.currentTimeMillis();
        System.out.println("protNorm--- " + formatter.format((end7 - end6) / (1000d * 60)) + " min.");

        if(SecondProcess){
            GenerateSingleSite();
            groupBy = 4;
        }

        if(param.abn_type==0)
        {
            if(protNorm>=0 && protNorm<=2){
                Report("Ratio"); //report ratios
                Report("RatioAbn"); //report abundances
            }
            else{
                Report("RatioAbn"); //report abundances
            }
        }
        else
        {
            Report("RawAbn"); //report abundances
        }

        gPsmMap.clear();
        gAbnMap.clear();

        long end8 = System.currentTimeMillis();
        System.out.println("Report--- " + formatter.format((end8 - end7) / (1000d * 60)) + " min.");
    }

    private static TreeMap<String, List<String>> LoadPsms(List<File> FileLi) throws IOException
    {
        TreeMap<String, List<String>> FileMap = new TreeMap<String, List<String>>();
        for(File f : FileLi)
        {
            try
            {
                String fName = f.getAbsolutePath();
                List<String> PsmLi = new ArrayList<String>();

                String NewPath = f.getAbsolutePath().replace(".tsv",".ti");
                BufferedReader br = new BufferedReader(new FileReader(NewPath));
                String title = br.readLine();
                param.TitleMap.put(fName, title);
                ds_Index indObj = param.indMap.get(fName);

                String line = "";
                while ((line = br.readLine()) != null)
                {
                    String[] strAry = line.split("\t");
                    boolean isUsed = Boolean.parseBoolean(strAry[indObj.isUsedIndex]);
                    double refInt = Double.parseDouble(strAry[indObj.refIndex]);
                    String peptide = strAry[indObj.pepcIndex];
                    String assignedMod = strAry[indObj.assignedModcIndex];
                    String ptmLocal = indObj.ptmLocalcIndex>=0? strAry[indObj.ptmLocalcIndex]: "";
                    String proteinID = strAry[indObj.proteinIDcIndex];
                    String gene = strAry[indObj.genecIndex];
                    double MS1Int = Double.parseDouble(strAry[indObj.ms1IntIndex]);

                    if(isUsed)
                    {
                        //region Get num. of phospho site with probablity > threshold ; and the start, end position

                        int slCount = 0;
                        String slpStr = ""; //site localized position
                        TreeMap<Integer, List<String>> slmMap = new TreeMap<Integer, List<String>>(); // key: position, value: site localized mass
                        List<String> slpLi = new ArrayList<String>();
                        int sCount = 0;
                        TreeMap<Integer, Double> probMap = new TreeMap<Integer, Double>();

                        proteinID = proteinID.contains("\"") ? proteinID.replaceAll("\"","") :  proteinID;
                        String ProtSeq = param.fastaMap.get(proteinID);
                        if (ProtSeq == null) {
                            System.out.printf("Error: could not find protein ID %s in database. Stopping analysis\n", proteinID);
                            System.exit(0);
                        }
                        int pepIndex = (ProtSeq.indexOf(peptide)>0)?ProtSeq.indexOf(peptide):0;
                        int sIndex = pepIndex;
                        int eIndex = pepIndex+peptide.length();

                        if(param.minSiteProb>0)
                        {
                            //region Get num. of possible modified site
                            String masstag = "("+param.columntag.substring(param.columntag.indexOf(":")+1)+")";
                            int index = assignedMod.indexOf(masstag);
                            while(index>=0)
                            {
                                sCount += 1;
                                index = assignedMod.indexOf(masstag, index+1);
                            }
                            //endregion

                            //region Get positions
                            int lindex = ptmLocal.indexOf("(");
                            int rindex = ptmLocal.indexOf(")");
                            int gap = 0;
                            boolean isFirst = true;
                            while(lindex>=0)
                            {
                                double prob = Double.parseDouble(ptmLocal.substring(lindex+1, rindex));
                                int key = -1;
                                if(isFirst)
                                {
                                    key = lindex;
                                    isFirst = false;
                                }
                                else
                                {
                                    gap += (rindex - lindex + 1);
                                    key = lindex - gap;
                                }
                                probMap.put(key, prob);
                                lindex = ptmLocal.indexOf("(", lindex+1);
                                rindex = ptmLocal.indexOf(")", rindex+1);
                            }
                            boolean isStart = false;
                            for(int key : probMap.keySet())
                            {
                                if(probMap.get(key) > param.minSiteProb)
                                {
                                    slCount += 1;
                                    String str = peptide.charAt(key-1) + "" + (pepIndex+key);
                                    slpStr  += str;
                                    slpLi.add(str);
                                }
                                if(!isStart)
                                {
                                    sIndex = key;
                                    isStart = true;
                                }
                                eIndex = key;
                            }

                            sIndex = (sIndex < 0) ? 0 : sIndex + pepIndex;
                            eIndex = (eIndex < 0) ? 0 : eIndex + pepIndex;
                            //endregion
                        }
                        else if(param.minSiteProb==0)
                        {
                            //region Get site position

                            List<Integer> positionLi = new ArrayList<Integer>();
                            String[] aModAry = assignedMod.split(",");
                            for(String aMod : aModAry){
                                String mod = TMTIntegrator.getAssignedModIndex(aMod, strAry[indObj.observedModIndex], param.useGlycoComposition);
                                if(param.modTagLi.contains(mod))
                                {
                                    int pos = Integer.valueOf(aMod.substring(0,aMod.indexOf("(")-1).trim());
                                    positionLi.add(pos);
                                    String tmp=mod.replace("(","").replace(")","");
                                    tmp = tmp.substring(1);     // remove first character (AA code) from index

                                    if(slmMap.containsKey(pos)){
                                        List<String> slmLi = slmMap.get(pos);
                                        slmLi.add(tmp);
                                    }
                                    else{
                                        List<String> slmLi = new ArrayList<String>();
                                        slmLi.add(tmp);
                                        slmMap.put(pos,slmLi);
                                    }
                                }
                            }
                            Collections.sort(positionLi);

                            for(int i=0;i<positionLi.size();i++)
                            {
                                String str = peptide.charAt(positionLi.get(i)-1) + "" + (pepIndex+positionLi.get(i));
                                slpStr  += str;
                                slpLi.add(str);
                            }
                            slCount = positionLi.size();
                            sCount = positionLi.size();
                            sIndex = positionLi.get(0) + pepIndex;
                            eIndex = positionLi.get(positionLi.size()-1) + pepIndex;

                            //endregion
                        }

                        //endregion

                        String groupkey = "";
                        //region  Set up key

                        if(groupBy==1){ //1: ProteinID
                            groupkey = proteinID;
                        }
                        else if (groupBy==2) { //2: Peptide
                            if((param.minSiteProb>=0)&&(sCount >0)){
                                groupkey = proteinID + "%" + sIndex + "%" + eIndex;
                            }
                            else if(param.minSiteProb<0)
                            {
                                groupkey = proteinID + "%" + peptide;
                            }
                            else{
                                groupkey = "null";
                            }
                        }
                        else if (groupBy==3){ //3. Multiple Phospho Sites
                            if(sCount >0)
                            {
                                groupkey = proteinID + "%" + sIndex + "%" + eIndex + "%" + sCount + "%" + slCount;
                                groupkey = (slpStr.length() > 0) ? (groupkey + "%" + slpStr) : groupkey;
                                groupkey = (peptide.split("-").length>1)? (groupkey+"%"+peptide.split("-")[1]): groupkey;

                                //Mark localized site in peptide sequence
                                String NewPepSeq = "";
                                for(int i= 0; i<peptide.length(); i++){
                                    NewPepSeq += (probMap.containsKey(i+1) && (probMap.get(i+1) > param.minSiteProb)) ?
                                            Character.toLowerCase(peptide.charAt(i)) : peptide.charAt(i);
                                }
                                peptide = NewPepSeq;
                            }
                            else{
                                groupkey = "null";
                            }
                        }
                        else if (groupBy==4) { //4. Single Phospho Sites
                            if (slpLi.size() > 0) {
                                String NewPepSeq = ""; //Mark localized site in peptide sequence
                                for(int i= 0; i<peptide.length(); i++){
                                    NewPepSeq += (probMap.containsKey(i+1) && (probMap.get(i+1) >= param.minSiteProb)) ?
                                            Character.toLowerCase(peptide.charAt(i)) : peptide.charAt(i);
                                }
                                peptide = NewPepSeq;
                            }
                            else{
                                groupkey = "null";
                            }
                        }
                        else if (groupBy==5){ //5. Multi-mass for glycosylation
                            if(sCount >0)
                            {
                                String slmStr="";
                                for(List<String> slmLi : slmMap.values()){
                                    for(String slm : slmLi){
                                        slmStr += slm+"_";
                                    }
                                }
                                slmStr = slmStr.substring(0,slmStr.length()-1);

                                groupkey = proteinID + "%" + sIndex + "%" + eIndex + "%" + sCount + "%" + slCount;
                                groupkey = (slpStr.length() > 0) ? (groupkey + "%" + slpStr) : groupkey;
                                groupkey += "%"+slmStr;

                                //Mark localized site in peptide sequence
                                String NewPepSeq = "";
                                for(int i= 0; i<peptide.length(); i++){
                                    NewPepSeq += (probMap.containsKey(i+1) && (probMap.get(i+1) > param.minSiteProb)) ?
                                            Character.toLowerCase(peptide.charAt(i)) : peptide.charAt(i);
                                }
                                peptide = NewPepSeq;
                            }
                            else{
                                groupkey = "null";
                            }
                        }
                        else if (groupBy==0) { //0. Gene
                            groupkey = gene;
                        }

                        //endregion

                        //region Use MS1 intensity
                        if(param.ms1Int){
                            double Sum = 0;
                            for(int x=indObj.abnIndex; x<strAry.length; x++){
                                Sum += Double.parseDouble(strAry[x]);
                            }
                            refInt = MS1Int * (refInt/Sum);
                        }
                        //endregion

                        if(groupkey != "null")
                        {
                            if((groupBy == 4)&&(slpLi.size() > 0))
                            { //4. Single Phospho Site
                                for(int i=0; i<slpLi.size(); i++)
                                {
                                    String NewPsm = proteinID +"%" + slpLi.get(i) + "#" + refInt + "#" + peptide + "#" + line;
                                    PsmLi.add(NewPsm);
                                }
                            }
                            else
                            {
                                String NewPsm = groupkey + "#" + refInt + "#" + peptide + "#" + line;
                                PsmLi.add(NewPsm);
                            }
                        }
                    }
                }
                br.close();
                FileMap.put(fName, PsmLi);
            }
            catch (Exception e)
            {
                System.out.println("Error at: "+f.getAbsolutePath()+" "+e);
                System.exit(0);
            }
        }

        return FileMap;
    }

    private static Boolean TryParseInt(String value)
    {
        try {
            Integer.parseInt(value);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    private static List<String> RtNorm(List<String> PsmLi, ds_Index indObj)
    {
        List<String> NewPsmLi = new ArrayList<String>();

        double minRt = 999999;
        double maxRt = 0;
        //region Get min. and max. Rt

        for(String psm : PsmLi){
            String[] strAry = psm.split("\t");
            double Rt = Double.parseDouble(strAry[indObj.rtIndex]);
            if(minRt>Rt){
                minRt = Rt;
            }
            if(maxRt < Rt){
                maxRt = Rt;
            }
        }

        //endregion

        int bin = 10;
        TreeMap<Double, List<String>> binMap = new TreeMap<Double, List<String>>();
        //region Generate bins
        double gap = (maxRt - minRt)/bin;
        double sRt = minRt;
        double eRt = maxRt + gap;
        while(sRt<eRt){
            binMap.put(sRt, new ArrayList<String>());
            sRt += gap;
        }
        binMap.put((maxRt+gap), new ArrayList<String>());
        //endregion

        List<Double> keyLi = new ArrayList<Double>(binMap.keySet());
        Collections.sort(keyLi);
        //region Assign psms to bins
        for(String psm : PsmLi)
        {
            String[] strAry = psm.split("\t");
            double Rt = Double.parseDouble(strAry[indObj.rtIndex]);
            for(int i=0; i<keyLi.size()-1;i++)
            {
                if((Rt>=keyLi.get(i)) && (Rt<keyLi.get(i+1))){
                    List<String> strLi = binMap.get(keyLi.get(i));
                    strLi.add(psm);
                    break;
                }
            }
        }
        //endregion

        //region Normalize each channel
        for(double key : binMap.keySet())
        {
            //System.out.print(key+"\t");
            List<String> SubPsmLi = binMap.get(key);
            NewPsmLi.addAll(PsmNor(SubPsmLi, indObj));
        }
        //endregion

        return NewPsmLi;
    }

    private static List<String> PsmNor(List<String> SubPsmLi, ds_Index indObj)
    {
        List<String> NewPsmLi = new ArrayList<String>();

        double[][] ratio2DAry = new double[SubPsmLi.size()][PlexNum];
        //region Convert to 2D array
        for(int i=0; i<SubPsmLi.size(); i++){
            String[] strAry = SubPsmLi.get(i).split("\t");
            for(int j=indObj.abnIndex; j<strAry.length; j++){
                ratio2DAry[i][j-indObj.abnIndex] = Double.parseDouble(strAry[j]);
            }
        }
        //endregion

        double[] MedAry = new double[PlexNum];
        //region Get median ratio for each channel
        for(int j=0; j<PlexNum; j++){
            List<Double> tmpLi = new ArrayList<Double>();
            for(int i=0; i<SubPsmLi.size(); i++){
                if(ratio2DAry[i][j] != -9999){
                    tmpLi.add(ratio2DAry[i][j]);
                }
            }
            MedAry[j] = TakeMedian(tmpLi);
        }
        //endregion

        //region Subtract median ratio from each channel
        for(int i=0; i<ratio2DAry.length; i++){
            for(int j=0; j<ratio2DAry[i].length; j++){
                ratio2DAry[i][j] = (ratio2DAry[i][j] != -9999) ? (ratio2DAry[i][j]-MedAry[j]) : ratio2DAry[i][j];
            }
        }
        //endregion

        //region Update new ratios
        for(int i=0; i<SubPsmLi.size(); i++) {
            String[] strAry = SubPsmLi.get(i).split("\t");
            String NewPsm = "";
            for(int j=0; j<strAry.length; j++){
                NewPsm += (j<indObj.abnIndex) ? (strAry[j] + "\t") : (ratio2DAry[i][j-indObj.abnIndex] + "\t");
            }
            NewPsmLi.add(NewPsm);
        }
        //endregion

        return NewPsmLi;
    }

    private static List<String> LogNor(List<String> PsmLi, ds_Index indObj)
    {
        List<String> NewPsmLi = new ArrayList<String>();
        for(int i = 0; i < PsmLi.size(); i++)
        {
            String[] strAry = PsmLi.get(i).split("\t");
            double refValue = Double.parseDouble(strAry[indObj.refIndex]) > 0 ? Log2(Double.parseDouble(strAry[indObj.refIndex])) : 0;
            double[] ratioAry = new double[PlexNum];
            for(int j = indObj.abnIndex; j < strAry.length; j++) {
                ratioAry[j-indObj.abnIndex] = (Double.parseDouble(strAry[j]) > 0) ?
                        (Log2(Double.parseDouble(strAry[j])) - refValue) : -9999;
            }
            String NewPsm = "";
            for (int j=0; j<strAry.length; j++)
            {
                NewPsm += (j<indObj.abnIndex) ? (strAry[j] + "\t") : (ratioAry[j-indObj.abnIndex] + "\t");
            }
            NewPsmLi.add(NewPsm);
        }
        return NewPsmLi;
    }

    private static void GroupPsm(TreeMap<String, List<String>> FileMap) throws IOException
    {
        for(String fName : FileMap.keySet())
        {
            ds_Index indObj = param.indMap.get(fName);
            List<String> PsmLi = FileMap.get(fName);
            for(String Psm : PsmLi){
                String[] strAry = Psm.split("\t");
                String[] strAry1 = strAry[0].split("#");
                String groupkey = strAry1[0];
                String NewPepSeq = strAry1[2];
                String gene = strAry[indObj.genecIndex];

                if(gPsmMap.containsKey(groupkey))
                {
                    TreeMap<String, ds_PsmInfo> fMap = gPsmMap.get(groupkey);
                    if(fMap.containsKey(fName))
                    {
                        ds_PsmInfo pi = fMap.get(fName);
                        pi.gene = gene;
                        pi.peptide = (pi.peptide.length() < NewPepSeq.length()) ? NewPepSeq : pi.peptide;
                        pi.PsmLi.add(Psm);
                    }
                    else
                    {
                        ds_PsmInfo pi = new ds_PsmInfo();
                        pi.gene = gene;
                        pi.peptide = NewPepSeq;
                        pi.PsmLi.add(Psm);
                        fMap.put(fName, pi);
                    }
                }
                else
                {
                    TreeMap<String, ds_PsmInfo> fMap = new TreeMap<String, ds_PsmInfo>();
                    ds_PsmInfo pi = new ds_PsmInfo();
                    pi.gene = gene;
                    pi.peptide = NewPepSeq;
                    pi.PsmLi.add(Psm);

                    fMap.put(fName, pi);
                    gPsmMap.put(groupkey, fMap);
                }
            }
        }

        //Compute the total reference intensity
        for(String groupkey : gPsmMap.keySet())
        {
            TreeMap<String, ds_PsmInfo> fMap = gPsmMap.get(groupkey);
            for(ds_PsmInfo pi : fMap.values())
            {
                List<Double> refIntLi = new ArrayList<Double>();
                for(String Psm : pi.PsmLi)
                {
                    String[] strAry = Psm.split("\t");
                    String[] strAry1 = strAry[0].split("#");
                    double refInt = Double.parseDouble(strAry1[1]);
                    refIntLi.add(refInt);
                    pi.tolrefInt += refInt;  //sum all peptides
                }
                Collections.sort(refIntLi);

                //Top3 intensive peptides
                if(param.top3Pep){
                    if(pi.PsmLi.size()>=3){
                        pi.tolrefInt = 0;
                        for(int i=refIntLi.size()-1; i>=refIntLi.size()-3; i--){
                            pi.tolrefInt += refIntLi.get(i);
                        }
                    }
                }
            }
        }
    }

    private static void RemoveOutlier()
    {
        for(String groupkey : gPsmMap.keySet())
        {
            TreeMap<String, ds_PsmInfo> fMap = gPsmMap.get(groupkey);
            for (String fName : fMap.keySet())
            {
                ds_PsmInfo pi = fMap.get(fName);
                ds_Index indObj = param.indMap.get(fName);
                if(pi.PsmLi.size()>3) { //only process pi with at least 4 psms, for speeding up
                    List<String> NewPsmLi = IqrPsm(pi.PsmLi, indObj);
                    pi.PsmLi.clear();
                    pi.PsmLi.addAll(NewPsmLi);
                }
            }
        }
    }

    private static List<String> IqrPsm(List<String> PsmLi, ds_Index indObj)
    {
        List<String> NewPsmLi = new ArrayList<String>();

        //region Get ratios
        double[][] Ratio2DAry = new double[PsmLi.size()][PlexNum];
        for(int i=0; i<PsmLi.size(); i++){
            String[] strAry = PsmLi.get(i).split("\t");
            for(int j=indObj.abnIndex; j<strAry.length; j++){
                Ratio2DAry[i][j-indObj.abnIndex] = Double.parseDouble(strAry[j]);
            }
        }
        //endregion

        //region Remove outlier from each channel
        for(int j =0; j < PlexNum; j++)
        {
            List<Double> tmpLi = new ArrayList<Double>();
            for (int i = 0; i< PsmLi.size(); i++){
                if(Ratio2DAry[i][j] != -9999)
                    tmpLi.add(Ratio2DAry[i][j]);
            }
            if(tmpLi.size()>3){ //use IQR for outlier removal
                double[] IqrAry = IqrRatio(tmpLi);
                double LowerIqr = IqrAry[0];
                double UpperIqr = IqrAry[1];

                for(int i=0; i<PsmLi.size(); i++){
                    if((Ratio2DAry[i][j] < LowerIqr) || (Ratio2DAry[i][j] > UpperIqr)){
                        Ratio2DAry[i][j] = -9999;
                    }
                }
            }
        }
        //endregion

        //region Update the ratios
        for(int i=0; i<PsmLi.size(); i++){
            String[] strAry = PsmLi.get(i).split("\t");
            String NewPsm = "";
            for(int j=0; j<strAry.length; j++){
                NewPsm += (j<indObj.abnIndex) ? (strAry[j] + "\t") : (Ratio2DAry[i][j-indObj.abnIndex] + "\t");
            }
            NewPsmLi.add(NewPsm);
        }
        //endregion

        return NewPsmLi;
    }

    private static double[] IqrRatio(List<Double> RatioLi) //Interquartile range method for outlier removal
    {
        Collections.sort(RatioLi);
        double[] IqrAry = new double[2]; //0: LowerIqr; 1: UpperIqr

        //region Get Q1 and Q3
        int n = RatioLi.size()/4;
        List<Double> Q1Li = new ArrayList<Double>();
        List<Double> Q3Li = new ArrayList<Double>();
        for(int i = 0; i <= n; i++){
            Q1Li.add(RatioLi.get(i));
        }
        for(int i = RatioLi.size()-n-1; i<RatioLi.size(); i++)
        {
            Q3Li.add(RatioLi.get(i));
        }
        //endregion

        double Q1 = TakeMedian(Q1Li);
        double Q3 = TakeMedian(Q3Li);
        IqrAry[0] = Q1 - 1.5*(Q3 - Q1);
        IqrAry[1] = Q3 + 1.5*(Q3 - Q1);

        return IqrAry;
    }

    private static void Collapse()
    {
        for(String groupkey : gPsmMap.keySet())
        {
            TreeMap<String, ds_PsmInfo> fMap = gPsmMap.get(groupkey);
            TreeMap<String, double[]> fAbnMap = new TreeMap<String, double[]>();
            List<String> ProtIdLi = new ArrayList<String>();
            int NumPsm = 0;
            double MaxPepProb = 0;
            List<String> geneLi = new ArrayList<String>();
            List<String> pepLi = new ArrayList<String>();
            for(String fName : fMap.keySet())
            {
                ds_PsmInfo pi = fMap.get(fName);
                ds_Index indObj = param.indMap.get(fName);
                double[] mAry = new double[TotLen]; //index=10 for SumAbn
                double[][] Ratio2DAry = new double[pi.PsmLi.size()][PlexNum];
                ds_Ratio[][] Obj2DAry = new ds_Ratio[pi.PsmLi.size()][PlexNum];

                //region Get  Max. peptide probability and all peptide sequences
                for(int i = 0; i < pi.PsmLi.size(); i++){
                    String[] strAry = pi.PsmLi.get(i).split("\t");
                    double PepProb = Double.parseDouble(strAry[indObj.pepProbcIndex]);
                    if(PepProb > MaxPepProb){
                        MaxPepProb = PepProb;
                    }
                    if(!pepLi.contains(pi.peptide)){
                        pepLi.add(pi.peptide);
                    }
                    if(!geneLi.contains(pi.gene)){
                        geneLi.add(pi.gene);
                    }
                }
                //endregion

                //region Get ratios
                if(param.aggregation_method==0)
                {
                    for(int i = 0; i < pi.PsmLi.size(); i++)
                    {
                        String[] strAry = pi.PsmLi.get(i).split("\t");
                        for(int j=indObj.abnIndex; j<strAry.length; j++){
                            Ratio2DAry[i][j-indObj.abnIndex] = Double.parseDouble(strAry[j]);
                        }
                    }
                }
                else if(param.aggregation_method==1)
                {
                    for(int i = 0; i < pi.PsmLi.size(); i++)
                    {
                        String[] strAry = pi.PsmLi.get(i).split("\t");
                        for(int j=indObj.abnIndex; j<strAry.length; j++)
                        {
                            ds_Ratio rObj = new ds_Ratio();
                            rObj.preInt = Double.parseDouble(strAry[indObj.ms1IntIndex]);
                            rObj.ratio = Double.parseDouble(strAry[j]);
                            Obj2DAry[i][j-indObj.abnIndex] = rObj;
                        }
                    }
                }

                //endregion

                //region Take median ratios
                for(int j =0; j < PlexNum; j++)
                {
                    if(param.aggregation_method==0)
                    {
                        List<Double> tmpLi = new ArrayList<Double>();
                        for (int i = 0; i< pi.PsmLi.size(); i++){
                            if(Ratio2DAry[i][j] != -9999)
                                tmpLi.add(Ratio2DAry[i][j]);
                        }
                        mAry[j] = TakeMedian(tmpLi);
                    }
                    else if(param.aggregation_method==1)
                    {
                        List<ds_Ratio> rObjLi = new ArrayList<ds_Ratio>();
                        for (int i = 0; i< pi.PsmLi.size(); i++){
                            if(Obj2DAry[i][j].ratio != -9999)
                                rObjLi.add(Obj2DAry[i][j]);
                        }
                        mAry[j] = TakeWeightedMedian(rObjLi);
                    }
                }
                //endregion

                mAry[PlexNum] = pi.tolrefInt;
                fAbnMap.put(fName, mAry);

                for(int i = 0; i < pi.PsmLi.size(); i++){ //Get ProteinID
                    String[] strAry = pi.PsmLi.get(i).split("\t");
                    String ProtId = strAry[indObj.proteinIDcIndex];
                    if(!ProtIdLi.contains(ProtId)){
                        ProtIdLi.add(ProtId);
                    }
                }

                if(groupBy == 1) {
                    NumPsm += pi.PsmLi.size();
                }
                else if(groupBy == 0){
                    NumPsm += pi.PsmLi.size();
                }
                //endregion
            }

            String ggpStr = ""; //global gene/peptide string
            for(String gene : geneLi)
            {
                ggpStr+=gene+";";
            }
            ggpStr=ggpStr.substring(0,ggpStr.lastIndexOf(";"));

            String ProtStr = "";
            for(String ProtId : ProtIdLi){
                ProtStr += ProtId+";";
            }
            ProtStr=ProtStr.substring(0, ProtStr.length()-1);

            if(groupBy==0){
                ggpStr = NumPsm +"\t" + ProtStr;
            }
            else if(groupBy==1){
                ggpStr = NumPsm + "\t" + ggpStr;
            }
            else{
                String pepStr = "";
                for(String pep: pepLi){
                    pepStr+=pep+";";
                }
                pepStr=pepStr.substring(0,pepStr.length()-1);
                ggpStr +="\t"+ProtStr+"\t"+pepStr;
            }
            groupkey = (ggpStr!="") ? (groupkey+"\t"+ggpStr+"\t"+MaxPepProb) : groupkey+"\t"+MaxPepProb;
            gAbnMap.put(groupkey, fAbnMap);
        }
    }

    private static void ProtNor()
    {
        if(protNorm==1 || protNorm==2){
            TreeMap<String, double[]> MedMap = GetProtMedian(false); //get median
            double m0= CalGlobalMed(MedMap);
            SubtractProt(1, MedMap, -1, -1); //subtract protein ratios

            if(protNorm==2){
                TreeMap<String, double[]> AbsMedMap = GetProtMedian(true);
                double Med0 = CalGlobalMed(AbsMedMap);
                SubtractProt(2, AbsMedMap, m0, Med0); //subtract protein ratios
            }
        }
        else if(protNorm==3){
            if(param.abn_type==0){//have to convert ratio2Abn if not the raw-based abundance
                Ratio2Abn();
            }
            SLNorm();
            IRSNorm();
        }
    }

    private static void Ratio2Abn()
    {
        double GloMinRefInt = CalGloMinRefInt();
        for(String groupkey : gAbnMap.keySet()){
            TreeMap<String, double[]> fAbnMap = gAbnMap.get(groupkey);

            double AvgAbn = 0;
            //region M2: Calculate average abundance (using GlobalMinRefInt)
            for(String fName : param.fNameLi)
            {
                double[] mAry = fAbnMap.containsKey(fName) ? fAbnMap.get(fName) : new double[TotLen];
                if(mAry [PlexNum]>0)
                {
                    AvgAbn += mAry [PlexNum];
                }
                else
                {
                    AvgAbn += GloMinRefInt;
                }
            }
            AvgAbn = AvgAbn/param.fNameLi.size();
            //endregion

            //region Convert ratio to abundance
            for(String fName : param.fNameLi)
            {
                double[] mAry = fAbnMap.containsKey(fName)?fAbnMap.get(fName):null;
                ds_Index indObj = param.indMap.get(fName);
                int refIndex = indObj.refIndex - indObj.abnIndex;
                if(mAry!=null){
                    for(int i=0; i<PlexNum; i++){
                        if(mAry[i]!=-9999 && i!=refIndex){
                            mAry[i] =  Log2(Math.pow(2, mAry[i]) * AvgAbn);
                        }
                    }
                }
            }
            //endregion
        }
    }

    private static void SLNorm()
    {
        List<String> fNameLi = new ArrayList<>(param.TitleMap.keySet());

        TreeMap<String, double[] > sumMap = new TreeMap<String,double[] >();
        //region 1. get the summation of all channels

        //Need to check reference channel, need to be excluded from the summation
        for(String fName : fNameLi){
            List<double[]> OrgLi = new ArrayList<double[]>();
            for(TreeMap<String, double[]> fAbnMap : gAbnMap.values()){
                if (fAbnMap.containsKey(fName)){
                    OrgLi.add(fAbnMap.get(fName));
                }
            }
            double[] sumAry = new double[TotLen];
            for(int i = 0; i < sumAry.length; i++) {
                double sum = 0;
                for(double[] mAry : OrgLi){
                    if(mAry[i] != -9999){
                        sum += mAry[i];
                    }
                }
                sumAry[i] = sum;
            }
            sumMap.put(fName, sumAry);
        }
        //endregion

        double sumAvg = 0;
        //region 2. compute sum average
        int count = 0;
        for (double[] sumAry : sumMap.values()){
            for(int i=0;i<sumAry.length;i++)
            {
                sumAvg += sumAry[i]>0?sumAry[i]:0;
                count += sumAry[i]>0? 1:0;
            }
        }
        sumAvg/=count;
        //endregion

        //region 3. adjust intensity using factors
        for(String fName : fNameLi){
            double[] sumAry = sumMap.get(fName);
            for(TreeMap<String, double[]> fAbnMap : gAbnMap.values()) {
                if (fAbnMap.containsKey(fName)) {
                    double[] mAry = fAbnMap.get(fName);
                    for(int i=0;i<TotLen;i++){
                        if(mAry[i]!=-9999){
                            mAry[i] = mAry[i]*(sumAvg/sumAry[i]);
                        }
                    }
                }
            }
        }
        //endregion
    }

    private static void IRSNorm()
    {
        double GloMinRefInt = CalGloMinRefInt();
        List<String> fNameLi = new ArrayList<>(param.TitleMap.keySet());
        for(String groupkey : gAbnMap.keySet()){
            double AvgRef=0;
            int count=0;
            TreeMap<String, double[]> fAbnMap = gAbnMap.get(groupkey);
            for(String fName : fNameLi){
                if(fAbnMap.containsKey(fName)){
                    double[] mAry = fAbnMap.get(fName);
                    AvgRef+=mAry[mAry.length-1]>0? Math.log(mAry[mAry.length-1]):0;
                    count+=mAry[mAry.length-1]>0? 1:0;
                }
            }
            AvgRef=Math.exp(AvgRef/count);
            for(String fName : fNameLi) {
                if (fAbnMap.containsKey(fName)) {
                    double[] mAry = fAbnMap.get(fName);
                    double fac = (mAry[mAry.length-1]>0)?AvgRef/mAry[mAry.length-1]:GloMinRefInt;
                    for(int i=0;i<param.channelNum;i++){
                        mAry[i] = (mAry[i]!=-9999) ? mAry[i]*fac:-9999;
                    }
                }
            }
        }
    }

    private static TreeMap<String, double[]> GetProtMedian(boolean UseAbsValue)
    {
        List<String> fNameLi = new ArrayList<>(param.TitleMap.keySet());
        TreeMap<String, double[]> MedMap = new TreeMap<String, double[]>();
        //region Get the median
        for(String fName : fNameLi)
        {
            //double[] MedianAry = new double[10];
            double[] MedianAry = new double[PlexNum];
            List<double[]> OrgLi = new ArrayList<double[]>();
            for(TreeMap<String, double[]> fAbnMap : gAbnMap.values()){
                if (fAbnMap.containsKey(fName)){
                    OrgLi.add(fAbnMap.get(fName));
                }
            }
            for(int i = 0; i < MedianAry.length; i++) {
                List<Double> tmpLi = new ArrayList<Double>();
                for(double[] mAry : OrgLi){
                    if(mAry[i] != -9999){
                        if(UseAbsValue){
                            tmpLi.add(Math.abs(mAry[i]));
                        }
                        else{
                            tmpLi.add(mAry[i]);
                        }
                    }
                }
                MedianAry[i] = TakeMedian(tmpLi);
            }
            MedMap.put(fName, MedianAry);
        }

        return MedMap;
    }

    private static Double CalGlobalMed(TreeMap<String, double[]> MedMap)
    {
        List<Double> tmpLi = new ArrayList<Double>();
        for(double[] mAry : MedMap.values()){
            for(double m : mAry){
                if(m != -9999){
                    tmpLi.add(m);
                }
            }
        }
        double Med0 = TakeMedian(tmpLi);
        return Med0;
    }

    private static void SubtractProt(int NormIndex, TreeMap<String, double[]> MedMap, double m0, double Med0)
    {
        for(String groupkey : gAbnMap.keySet()){
            TreeMap<String, double[]> fAbnMap = gAbnMap.get(groupkey);
            for(String fName : fAbnMap.keySet()){
                double[] MedianAry = MedMap.get(fName);
                double[] mAry =  fAbnMap.get(fName);
                for(int i = 0 ; i < MedianAry.length; i++){
                    if(mAry[i] != -9999){
                        if(NormIndex==1){
                            mAry[i] = mAry[i] - MedianAry[i];
                        }
                        else if(NormIndex==2){
                            //mAry[i] = (mAry[i]/MedianAry[i])*Med0+m0;
                            mAry[i] = (mAry[i]/MedianAry[i])*Med0;
                        }
                    }
                }
            }
        }
    }

    private static void Report(String flag) throws IOException
    {
        String gTag="";
        //region Change groupBy name
        if(groupBy==0){
            gTag="gene";
        }else if(groupBy==1){
            gTag="protein";
        }else if(groupBy==2){
            gTag="peptide";
        }else if(groupBy==3){
            gTag="multi-site";
        }else if(groupBy==4){
            gTag="single-site";
        }else if(groupBy==5){
            gTag="multi-mass";
        }
        //endregion

        String pnTag="";
        //region Change protNorm name
        if(protNorm==0){
            pnTag="None";
        }else if(protNorm==1){
            //pnTag="MC";
            pnTag="MD";
        }else if(protNorm==2){
            pnTag="GN";
        }else if(protNorm==3){
            pnTag="SL+IRS";
        }
        //endregion

        String path = "";
        if(flag == "Ratio"){
            path = param.reportPath + File.separator +"ratio_" + gTag +"_"+pnTag+".tsv";
        }
        else if(flag == "RatioAbn"){
            path = param.reportPath + File.separator  + "abundance_" + gTag +"_"+pnTag+".tsv";
        }
        else if(flag == "RawAbn"){
            path = param.reportPath + File.separator  + "raw2abundance_" + gTag +"_"+pnTag+".tsv";
        }

        BufferedWriter wr = new BufferedWriter(new FileWriter(path));
        //region Write title
        if(groupBy==1){ //1: Protin ID
            wr.write("Index\tNumberPSM\tGene");
        }
        else if (groupBy==2){ //2: Peptide
            wr.write("Index\tGene\tProteinID\tPeptide");
        }
        else if ((groupBy==3) || (groupBy==4) || (groupBy==5)){ //3. Site
            wr.write("Index\tGene\tProteinID\tPeptide");
        }
        else if (groupBy==0){ //0. Gene
            wr.write("Index\tNumberPSM\tProteinID");
        }
        wr.write("\tMaxPepProb\tReferenceIntensity");
        for(String fName : param.fNameLi)
        {
            ds_Index indObj = param.indMap.get(fName);
            String[] tAry = param.TitleMap.get(fName).split("\t");
            for(int i = indObj.abnIndex; i < tAry.length; i++)
            {
                if(!tAry[i].contains(param.refTag)){
                    wr.write( "\t" + tAry[i].replace(" Abundance",""));
                }
            }
        }
        if(param.print_RefInt){
            for(String fName : param.fNameLi)
            {
                ds_Index indObj = param.indMap.get(fName);
                File f =new File(fName);
                String[] tAry = param.TitleMap.get(fName).split("\t");

                String fldstr=f.getParent().substring(f.getParent().lastIndexOf("\\")+1);
                wr.write("\tRefInt_"+ (param.add_Ref<0?tAry[indObj.refIndex].replace(" Abundance",""):fldstr));
            }
        }
        wr.newLine();

        //endregion

        double GloMinRefInt = CalGloMinRefInt();

        //region write ratios / abundances
        for(String groupkey : gAbnMap.keySet()) {
            boolean isPrint = true;

            TreeMap<String, double[]> fAbnMap = gAbnMap.get(groupkey);
            double AvgAbn = 0;

            //region M2: Calculate average abundance (using GlobalMinRefInt)
            for(String fName : param.fNameLi)
            {
                double[] mAry = fAbnMap.containsKey(fName) ? fAbnMap.get(fName) : new double[TotLen];
                if(mAry [PlexNum]>0)
                {
                    AvgAbn += mAry [PlexNum];
                }
                else
                {
                    AvgAbn += GloMinRefInt;
                }
            }
            AvgAbn = AvgAbn/param.fNameLi.size();
            //endregion

            String[] strAry = groupkey.split("\t");
            double MaxPepProb = Double.parseDouble(strAry[strAry.length-1]);
            if((AvgAbn<=0) || (MaxPepProb<param.max_pep_prob_thres)){
                isPrint=false;
            }

            //region Print front columns

            if(isPrint){
                wr.write(groupkey.replace("%","_"));
            }

            //endregion

            if(isPrint){
                if(protNorm<3){
                    if(flag == "RawAbn"){
                        wr.write("\t" + AvgAbn);
                    }
                    else{
                        wr.write("\t" + Log2(AvgAbn));
                    }
                }
                else{
                    wr.write("\t" + AvgAbn);
                }

                String AbnStr = "";
                for(String fName : param.fNameLi)
                {
                    ds_Index indObj = param.indMap.get(fName);
                    int refIndex = indObj.refIndex - indObj.abnIndex;
                    double[] NullAry = new double[TotLen];
                    for (int x = 0; x < NullAry.length; x++) {
                        NullAry[x] = -9999;
                    }
                    double[] mAry = fAbnMap.containsKey(fName) ? fAbnMap.get(fName) : NullAry;

                    //region Record reference abundances
                    if (mAry[PlexNum] > 0) {
                        if(protNorm<3){
                            if(flag == "RawAbn"){
                                AbnStr += "\t" + mAry[PlexNum];
                            }
                            else{
                                AbnStr += "\t" + Log2(mAry[PlexNum]);
                            }
                        }
                        else{
                            AbnStr += "\t" + mAry[PlexNum];
                        }
                    } else {
                        AbnStr += "\tNA";
                    }
                    //endregion

                    if (protNorm < 3) { //none,MC, GN
                        if (flag == "Ratio") {
                            for (int i = 0; i < mAry.length - 1; i++) {
                                if (i != refIndex) {
                                    if (mAry[i] == -9999) {
                                        wr.write("\tNA");
                                    } else {
                                        wr.write("\t" + mAry[i]);
                                    }
                                }
                            }
                        } else if (flag == "RatioAbn") {
                            for (int i = 0; i < mAry.length - 1; i++) {
                                if (i != refIndex) {
                                    if (mAry[i] == -9999) {
                                        wr.write("\tNA");
                                    } else {
                                        double abn = Log2(Math.pow(2, mAry[i]) * AvgAbn);
                                        wr.write("\t" + abn);
                                    }
                                }
                            }
                        }
                        else if (flag == "RawAbn") {
                            for (int i = 0; i < mAry.length - 1; i++) {
                                if (i != refIndex) {
                                    if (mAry[i] == -9999) {
                                        wr.write("\tNA");
                                    } else {
                                        wr.write("\t" + mAry[i]);
                                    }
                                }
                            }
                        }
                    } else { //SL, IRS, SL+IRS
                        for (int i = 0; i < mAry.length - 1; i++) {
                            if (i != refIndex) {
                                if (mAry[i] == -9999) {
                                    wr.write("\tNA");
                                } else {
                                    wr.write("\t" + mAry[i]);
                                }
                            }
                        }
                    }
                }

                if(param.print_RefInt){
                    wr.write(AbnStr);
                }
                wr.newLine();
            }
        }
        //endregion

        wr.close();
    }

    private static double CalGloMinRefInt()
    {
        double GloMinRefInt = 0;
        boolean isFirst = true;
        for(String groupkey : gAbnMap.keySet()) {
            TreeMap<String, double[]> fAbnMap = gAbnMap.get(groupkey);
            for(double[] mAry : fAbnMap.values()){
                if(isFirst) {//Assign initial intensity
                    GloMinRefInt =mAry[PlexNum];
                    isFirst = false;
                }
                if((mAry[PlexNum] < GloMinRefInt) && (mAry[PlexNum] > 0)) {
                    GloMinRefInt = mAry[PlexNum];
                }
            }
        }
        return GloMinRefInt;
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

    private static double TakeWeightedMedian(List<ds_Ratio> rObjLi)
    {
        double mValue = -9999;
        if(rObjLi.size() == 0)
        {
            mValue = -9999;
        }
        else if(rObjLi.size() == 2)
        {
            mValue = (rObjLi.get(0).ratio + rObjLi.get(1).ratio)/2;
        }
        else
        {
            double sum=0;
            double pow=1;
            //1. Precursor Intensity^pow
            for(int i=0; i<rObjLi.size(); i++)
            {
                rObjLi.get(i).weight = Math.pow(rObjLi.get(i).preInt, pow);
                sum += rObjLi.get(i).weight;
            }

            //2. The normalization of weight
            for(int i=0;i<rObjLi.size();i++)
            {
                rObjLi.get(i).weight = rObjLi.get(i).weight/sum;
            }

            //3. Sort ratios
            Collections.sort(rObjLi, new Comparator<ds_Ratio>() {
                @Override
                public int compare(ds_Ratio r1, ds_Ratio r2) {
                    return Double.compare(r1.ratio, r2.ratio);
                }
            });

            //4. take the weighted median
            int index=1;
            while(index<rObjLi.size())
            {
                double w1 = 0;
                double w2 = 0;
                for(int i =0; i < index; i++)
                {
                    w1 += rObjLi.get(i).weight;
                }
                for(int i = index+1;i< rObjLi.size(); i++)
                {
                    w2 += rObjLi.get(i).weight;
                }
                if((w1<=0.5)&&(w2<=0.5))
                {
                    break;
                }
                index+=1;
            }

            if(rObjLi.size()==index)
            {
                Collections.sort(rObjLi, new Comparator<ds_Ratio>() {
                    @Override
                    public int compare(ds_Ratio r1, ds_Ratio r2) {
                        return Double.compare(r1.weight, r2.weight);
                    }
                });
                mValue = rObjLi.get(rObjLi.size()-1).ratio;
            }
            else
            {
                mValue=rObjLi.get(index).ratio;
            }
        }

        return mValue;
    }

    private static double Log2(double value)
    {
        return (Math.log(value)/Math.log(2));
    }

    private static double TryParseIInt(String str)
    {
        double value = 0;
        try
        {
            value = Integer.parseInt(str);
        }
        catch (NumberFormatException e)
        {
            value = -1;
        }
        return value;
    }

    private static void GenerateSingleSite()
    {
        TreeMap<String, List<String>> keyMap = new TreeMap<String, List<String>>();
        int location=5;
        //region Cluster keys based on the index
        for(String groupkey : gAbnMap.keySet()) {
            String[] strAry = groupkey.split("\t");
            String[] indAry = strAry[0].split("%");
            if(TryParseIInt(indAry[indAry.length-1])<0){
                String[] sAry = indAry[location].split(param.modAA);
                for(int i=1; i<sAry.length; i++){
                    String newgroupkey = indAry[0]+"%"+indAry[location].charAt(indAry[location].indexOf(sAry[i])-1)+sAry[i];
                    if(keyMap.containsKey(newgroupkey)){
                        List<String> strLi = keyMap.get(newgroupkey);
                        if(!strLi.contains(groupkey)){
                            strLi.add(groupkey);
                        }
                    }
                    else{
                        List<String> strLi = new ArrayList<String>();
                        strLi.add(groupkey);
                        keyMap.put(newgroupkey, strLi);
                    }
                }
            }
        }
        //endregion

        //region Remove multiply sites, if singly site exists
        for(String key : keyMap.keySet()){
            Boolean isExist = false;
            List<String> strLi = keyMap.get(key);
            if(strLi.size()>1){
                for(String str : strLi){ //check the existence of singly site
                    String[] strAry = str.split("\t|%");
                    if(strAry[4].equalsIgnoreCase("1")){
                        isExist = true;
                        break;
                    }
                }
                if(isExist){
                    List<String> newStrLi = new ArrayList<String>();
                    for(String str : strLi) { //check the existence of singly site
                        String[] strAry = str.split("\t|%");
                        if(strAry[4].equalsIgnoreCase("1")){
                            newStrLi.add(str);
                        }
                    }
                    strLi.clear();
                    strLi.addAll(newStrLi);
                }
            }
        }
        //endregion

        TreeMap<String, TreeMap<String, double[]>> NewgAbnMap = new TreeMap<String, TreeMap<String, double[]>>();
        //region Calculate the median abundance
        for(String key : keyMap.keySet()){
            List<String> strLi = keyMap.get(key);
            if(strLi.size()>1){
                TreeMap<String, List<double[]>> gfAbnMap = new TreeMap<String, List<double[]>>();
                for(String str : strLi){
                    TreeMap<String, double[]> fAbnMap = gAbnMap.get(str);
                    for(String fName : param.fNameLi) {
                        double[] mAry = fAbnMap.containsKey(fName) ? fAbnMap.get(fName) : null;
                        if(mAry!=null){
                            if (gfAbnMap.containsKey(fName)){
                                List<double[]> mLi = gfAbnMap.get(fName);
                                mLi.add(mAry);
                            }
                            else{
                                List<double[]> mLi = new ArrayList<double[]>();
                                mLi.add(mAry);
                                gfAbnMap.put(fName, mLi);
                            }
                        }
                    }
                }
                TreeMap<String, double[]> upAbnMap = new TreeMap<String, double[]>();
                for(String fName : gfAbnMap.keySet() ){
                    List<double[]> mLi = gfAbnMap.get(fName);
                    if(mLi.size()>1){
                        double[] fmAry = new double[TotLen];
                        for(int i=0;i<TotLen;i++){
                            List<Double> vLi = new ArrayList<Double>();
                            for(int j=0;j<mLi.size();j++){
                                double[] mAry = mLi.get(j);
                                if(mAry[i]!=-9999)
                                {
                                    vLi.add(mAry[i]);
                                }
                            }
                            fmAry[i] = TakeMedian(vLi);
                        }
                        upAbnMap.put(fName, fmAry);
                    }
                    else{
                        upAbnMap.put(fName, mLi.get(0));
                    }
                }

                String gene = "";
                String protID = "";
                String peptides = "";
                double maxProb = -1000;
                for(String str : strLi){
                    String[] sAry = str.split("\t");
                    gene = sAry[1];
                    protID = sAry[2];
                    peptides += sAry[3] + ";";
                    maxProb = maxProb<Double.parseDouble(sAry[4])? Double.parseDouble(sAry[4]) : maxProb;
                }
                peptides=peptides.substring(0,peptides.length()-1);
                String newKey = key+ "\t" + gene + "\t" + protID +"\t"+peptides + "\t" + maxProb;
                NewgAbnMap.put(newKey, upAbnMap);
            }
            else{
                String newKey = key + strLi.get(0).substring(strLi.get(0).indexOf("\t"),strLi.get(0).length()-1);
                TreeMap<String, double[]> fAbnMap = gAbnMap.get(strLi.get(0));
                NewgAbnMap.put(newKey, fAbnMap);
            }
        }
        //endregion

        gAbnMap.clear();
        gAbnMap = NewgAbnMap;
    }
}
