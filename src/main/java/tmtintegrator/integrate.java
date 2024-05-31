package tmtintegrator;

import tmtintegrator.constants.GroupBy;
import tmtintegrator.constants.NormType;
import tmtintegrator.integrator.PsmFileLoader;
import tmtintegrator.integrator.PsmNormalizer;
import tmtintegrator.integrator.PsmProcessor;
import tmtintegrator.pojo.*;
import tmtintegrator.utils.Utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class integrate
{
    private static final Pattern pattern = Pattern.compile("([^%]+)%([ncA-Z])(\\d+)");
    private static final Pattern pattern2 = Pattern.compile("([ncA-Z])(\\d+)");
    private static Parameters param = new Parameters();
    private static int groupBy = -1;
    private static int protNorm = -1;
    private static Map<String, Map<String, PsmInfo>> gPsmMap = new TreeMap<>();
    private static Map<String, Map<String, double[]>> gAbnMap = new TreeMap<>();

    public static void run(Parameters parameter, int gb, int pn) throws IOException
    {
        param = parameter;
        groupBy = gb;
        protNorm = pn;

        Boolean SecondProcess = false;
        if(gb==4){
            // if single-site, need to do multi-site first
            groupBy = 3;
            SecondProcess = true;
        }

        long end1 = System.currentTimeMillis();
        NumberFormat formatter = new DecimalFormat("#0.00000");

        //Load all psm.tsvs in FileMap
        PsmFileLoader loader = new PsmFileLoader(param);
        GroupBy groupByEnum = GroupBy.fromValue(groupBy); // TODO: temp implementation, to be refactored
        Map<String, List<String>> FileMap = loader.loadPsmFiles(groupByEnum); // key: file path, value: list of psm (0 is title)

        long  end2 = System.currentTimeMillis();
        System.out.println("LoadPsms--- " + formatter.format((end2 - end1) / (1000d * 60)) + " min.");

        NormType normTypeEnum = NormType.fromValue(protNorm); // TODO: temp implementation, to be refactored
        PsmNormalizer normalizer = new PsmNormalizer(param, normTypeEnum);

        if(param.abn_type==0)
        {
            normalizer.logNormalizeData(FileMap);
        }

        long  end3 = System.currentTimeMillis();
        System.out.println("Take log and normalize--- " + formatter.format((end3 - end2) / (1000d * 60)) + " min.");

        if(param.psmNorm){ //PSM normalization
            normalizer.rtNormalizeData(FileMap);
        }

        long  end4 = System.currentTimeMillis();
        System.out.println("PSM normalization--- " + formatter.format((end4 - end3) / (1000d * 60)) + " min.");

        PsmProcessor processor = new PsmProcessor(param, GroupBy.fromValue(groupBy)); // TODO: temp implementation, to be refactored
        processor.groupPsm(FileMap);
        gPsmMap = processor.getGroupPsmMap(); // TODO: temp implementation, to be removed.
        FileMap.clear();
        System.out.println("Group psm: " + (System.currentTimeMillis() - end4) + " ms");

        if(param.outlierRemoval){
            processor.removeOutlier();
        }

        long  end5 = System.currentTimeMillis();
        System.out.println("outlierRemoval--- " + formatter.format((end5 - end4) / (1000d * 60)) + " min.");

        processor.collapse();
        gAbnMap = processor.getGroupAbundanceMap(); // TODO: temp implementation, to be removed.
        long  end6 = System.currentTimeMillis();
        System.out.println("Collapse--- " + formatter.format((end6 - end5) / (1000d * 60)) + " min.");

        normalizer.setGroupAbundanceMap(gAbnMap);
        if(protNorm > 0){
            normalizer.proteinNormalize();
        }

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

        Files.createDirectories(Paths.get(param.reportPath));

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
            //wr.write("Index\tGene\tProteinID\tPeptide");
            wr.write("Index\tGene\tProteinID\tPeptide\tSequenceWindow\tStart\tEnd");
        }
        else if ((groupBy==3) || (groupBy==4) || (groupBy==5)){ //3. Site
            //wr.write("Index\tGene\tProteinID\tPeptide\tSequenceWindow");
            wr.write("Index\tGene\tProteinID\tPeptide\tSequenceWindow\tStart\tEnd");
        }
        else if (groupBy==0){ //0. Gene
            wr.write("Index\tNumberPSM\tProteinID");
        }
        wr.write("\tMaxPepProb\tReferenceIntensity");
        for(String fName : param.fNameLi)
        {
            Index indObj = param.indMap.get(fName);
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
                Index indObj = param.indMap.get(fName);
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

            Map<String, double[]> fAbnMap = gAbnMap.get(groupkey);
            double AvgAbn = 0;

            //region M2: Calculate average abundance (using GlobalMinRefInt)
            for(String fName : param.fNameLi)
            {
                Index indObj = param.indMap.get(fName);
                double[] mAry = fAbnMap.containsKey(fName) ? fAbnMap.get(fName) : new double[indObj.totLen];
                if(mAry [indObj.plexNum]>0)
                {
                    AvgAbn += mAry [indObj.plexNum];
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

            if (groupBy == 4) {
                int site = -1;
                String[] ss = groupkey.split("\t");
                Matcher matcher = pattern.matcher(ss[0]);
                if (matcher.matches()) {
                    site = Integer.parseInt(matcher.group(3));
                } else {
                    System.err.println("Error: cannot parse site from " + ss[0]);
                    System.exit(1);
                }

                int offset = site - Integer.parseInt(ss[5]);
                String s = ss[4].replace(".", "");

                int dotIndex = ss[4].indexOf('.');
                if (dotIndex == -1) {
                    System.err.println("There are no find dots in " + s);
                    System.exit(1);
                }

                int startIndex = dotIndex + offset - 7;
                int endIndex = dotIndex + offset + 8;

                if (startIndex < 0) {
                    s = "_".repeat(-startIndex) + s;
                    endIndex -= startIndex;
                    startIndex = 0;
                }
                if (endIndex > s.length()) {
                    s += "_".repeat(endIndex - s.length());
                }

                String sequenceWindow =  s.substring(startIndex, startIndex + 7) + Character.toLowerCase(s.charAt(startIndex + 7)) + s.substring(startIndex + 8, endIndex);

                ss[4] = sequenceWindow;
                groupkey = String.join("\t", ss);
            }

            if(isPrint){
                wr.write(groupkey.replace("%","_"));
            }

            if(isPrint){
                if(protNorm<3){
                    if(flag == "RawAbn"){
                        wr.write("\t" + AvgAbn);
                    }
                    else{
                        wr.write("\t" + (param.log2transformed?Log2(AvgAbn):AvgAbn));
                    }
                }
                else{
                    wr.write("\t" + AvgAbn);
                }

                String AbnStr = "";
                for(String fName : param.fNameLi)
                {
                    Index indObj = param.indMap.get(fName);
                    int refIndex = indObj.refIndex - indObj.abnIndex;
                    double[] NullAry = new double[indObj.totLen];
                    for (int x = 0; x < NullAry.length; x++) {
                        NullAry[x] = -9999;
                    }
                    double[] mAry = fAbnMap.containsKey(fName) ? fAbnMap.get(fName) : NullAry;

                    //region Record reference abundances
                    if (mAry[indObj.plexNum] > 0) {
                        if(protNorm<3){
                            if(flag == "RawAbn"){
                                AbnStr += "\t" + mAry[indObj.plexNum];
                            }
                            else{
                                AbnStr += "\t" + (param.log2transformed?Log2(mAry[indObj.plexNum]): mAry[indObj.plexNum]);
                            }
                        }
                        else{
                            AbnStr += "\t" + mAry[indObj.plexNum];
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
                                        wr.write("\t" + (param.log2transformed?mAry[i]:Math.pow(2, mAry[i])));
                                    }
                                }
                            }
                        } else if (flag == "RatioAbn") {
                            for (int i = 0; i < mAry.length - 1; i++) {
                                if (i != refIndex) {
                                    if (mAry[i] == -9999) {
                                        wr.write("\tNA");
                                    } else {
                                        double abn = param.log2transformed?Log2(Math.pow(2, mAry[i]) * AvgAbn): Math.pow(2, mAry[i]) * AvgAbn;
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
                    }
                    else { //SL, IRS, SL+IRS
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
            Map<String, double[]> fAbnMap = gAbnMap.get(groupkey);
            for(String fName : fAbnMap.keySet()){
                Index indObj = param.indMap.get(fName);
                double[] mAry = fAbnMap.get(fName);
                if(isFirst) {//Assign initial intensity
                    GloMinRefInt =mAry[indObj.plexNum];
                    isFirst = false;
                }
                if((mAry[indObj.plexNum] < GloMinRefInt) && (mAry[indObj.plexNum] > 0)) {
                    GloMinRefInt = mAry[indObj.plexNum];
                }
            }
        }
        return GloMinRefInt;
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
        TreeMap<String, TreeMap<String, Integer>> keyPepMap = new TreeMap<String, TreeMap<String, Integer>>(); //to store the peptide start position
        int location=5;
        //region Cluster keys based on the index
        for(String groupkey : gAbnMap.keySet()) {
            String[] strAry = groupkey.split("\t");
            String[] indAry = strAry[0].split("%");

            if(TryParseIInt(indAry[indAry.length-1])<0){
                String[] sAry = indAry[location].split(param.modAA);

                //Find the peptide start position in the protein sequence
                int firstIndex = -1;
                String pepseq = strAry[3];
                for(int j=0; j<pepseq.length(); j++){
                    if(Character.isLowerCase(pepseq.charAt(j))){
                        firstIndex = j;
                        break;
                    }
                }
                int pepsIndex = Integer.parseInt(sAry[1]) - firstIndex;

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

                    int pepIndex = Integer.parseInt(sAry[i]) - pepsIndex;
                    if(keyPepMap.containsKey(newgroupkey)){
                        TreeMap<String, Integer> indMap = keyPepMap.get(newgroupkey);
                        indMap.put(pepseq, pepIndex);
                    }
                    else{
                        TreeMap<String, Integer> indMap = new TreeMap<String, Integer>();
                        indMap.put(pepseq, pepIndex);
                        keyPepMap.put(newgroupkey, indMap);
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

        Map<String, Map<String, double[]>> NewgAbnMap = new TreeMap<>();
        //region Calculate the median abundance
        for(String key : keyMap.keySet()){
            List<String> strLi = keyMap.get(key);
            if(strLi.size()>1){
                TreeMap<String, List<double[]>> gfAbnMap = new TreeMap<String, List<double[]>>();
                for(String str : strLi){
                    Map<String, double[]> fAbnMap = gAbnMap.get(str);
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
                    Index indObj = param.indMap.get(fName);
                    List<double[]> mLi = gfAbnMap.get(fName);
                    if(mLi.size()>1){
                        double[] fmAry = new double[indObj.totLen];
                        for(int i=0;i<indObj.totLen;i++){
                            List<Double> vLi = new ArrayList<Double>();
                            for(int j=0;j<mLi.size();j++){
                                double[] mAry = mLi.get(j);
                                if(mAry[i]!=-9999)
                                {
                                    vLi.add(mAry[i]);
                                }
                            }
                            fmAry[i] = Utils.takeMedian(vLi);
                        }
                        upAbnMap.put(fName, fmAry);
                    }
                    else{
                        upAbnMap.put(fName, mLi.get(0));
                    }
                }

                String str = strLi.get(0);
                String[] sAry = str.split("\t");
                String gene = sAry[1];
                String proteinID = sAry[2];
                Set<String> pepLi = new TreeSet<>();
                String extpepLi = sAry[4];
                String start = sAry[5];
                String end = sAry[6];
                double maxProb = 0;

                for (String s : strLi) {
                    String[] ss = s.split("\t");
                    pepLi.add(ss[3]);
                    maxProb = Math.max(maxProb, Float.parseFloat(ss[7]));
                }

                String newKey = key+ "\t" + gene + "\t" + proteinID + "\t" + String.join(";", pepLi) + "\t" + extpepLi + "\t" + start + "\t" + end + "\t" + maxProb;
                NewgAbnMap.put(newKey, upAbnMap);
            }
            else{
                String[] sAry = strLi.get(0).split("\t");
                String newKey = key + "\t" + sAry[1] + "\t" + sAry[2] + "\t" + sAry[3] + "\t" + sAry[4] + "\t" + sAry[5] + "\t" + sAry[6] + "\t" + sAry[7];
                Map<String, double[]> fAbnMap = gAbnMap.get(strLi.get(0));
                NewgAbnMap.put(newKey, fAbnMap);
            }
        }
        //endregion

        gAbnMap.clear();
        gAbnMap = NewgAbnMap;
    }
}
