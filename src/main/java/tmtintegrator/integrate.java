package tmtintegrator;

import tmtintegrator.constants.GroupBy;
import tmtintegrator.constants.NormType;
import tmtintegrator.integrator.PsmFileLoader;
import tmtintegrator.integrator.PsmNormalizer;
import tmtintegrator.integrator.PsmProcessor;
import tmtintegrator.integrator.ReportGenerator;
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
    private static Parameters param = new Parameters();
    private static int groupBy = -1;
    private static int protNorm = -1;
//    private static Map<String, Map<String, PsmInfo>> gPsmMap = new TreeMap<>();
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
//        gPsmMap = processor.getGroupPsmMap(); // TODO: temp implementation, to be removed.
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
            processor.generateSingleSite();
            groupBy = 4; // FIXME: not a good practice, should have a flog to indicate the second process
        }

        ReportGenerator reporter = new ReportGenerator(param, GroupBy.fromValue(groupBy), NormType.fromValue(protNorm), gAbnMap); // FIXME: temp implementation, to be refactored
        reporter.generateReport();

        gAbnMap.clear();

        long end8 = System.currentTimeMillis();
        System.out.println("Report--- " + formatter.format((end8 - end7) / (1000d * 60)) + " min.");
    }
}
