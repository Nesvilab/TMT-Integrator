package tmtintegrator.pojo;

import tmtintegrator.constants.GroupBy;
import tmtintegrator.constants.ReportInfoIndex;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Columns for report info part.
 */
public class ReportInfo {

    public static final String HEADER = "Index\tGene\tProteinID\tPeptide\tSequenceWindow\tStart\tEnd\tMaxPepProb\tSpectrum Number\tProtein\tEntry Name\tProtein Description\tMapped Genes\tMapped Proteins";
    public static final String HEADER_PROTEIN = "Index\tNumberPSM\tGene\tMaxPepProb\tProtein\tProtein ID\tEntry Name\tProtein Description\tOrganism\tIndistinguishable Proteins";
    public static final String HEADER_MODIFIED_PEPTIDE = "Index\tGene\tProteinID\tPeptide\tAssigned Modification\tSequenceWindow\tStart\tEnd\tMaxPepProb\tSpectrum Number\tProtein\tEntry Name\tProtein Description\tMapped Genes\tMapped Proteins";
    private final String index;
    private String numPSM;
    private final String gene;
    private String proteinID;
    private String peptide;
    private String sequenceWindow;
    private String start;
    private String end;
    private final String maxPepProb;
    private String assignedModifications;

    public ReportInfo(String[] parts, GroupBy groupBy) {
        if (groupBy == GroupBy.PROTEIN_ID) {
            this.index = parts[ReportInfoIndex.INDEX];
            this.numPSM = parts[ReportInfoIndex.NUMBER_PSM];
            this.gene = parts[ReportInfoIndex.GENE_PROTEIN];
            this.maxPepProb = parts[ReportInfoIndex.MAX_PEP_PROB_PROTEIN];
        } else {
            this.index = parts[ReportInfoIndex.INDEX].replace("%", "_");
            this.gene = parts[ReportInfoIndex.GENE];
            this.proteinID = parts[ReportInfoIndex.PROTEIN_ID];
            this.peptide = parts[ReportInfoIndex.PEPTIDE];
            this.sequenceWindow = parts[ReportInfoIndex.SEQUENCE_WINDOW];
            this.start = parts[ReportInfoIndex.START];
            this.end = parts[ReportInfoIndex.END];
            this.maxPepProb = parts[ReportInfoIndex.MAX_PEP_PROB];
            if (ReportInfoIndex.ASSIGNED_MODIFICATIONS < parts.length) {
                this.assignedModifications = parts[ReportInfoIndex.ASSIGNED_MODIFICATIONS];
            }
        }
    }

    /**
     * Propagate extra info for peptide and site level report from original psm file.
     *
     * @param extraPsmInfoMap <peptide, extra columns>
     */
    public void propagateExtraInfo(Map<String, ExtraPsmInfo> extraPsmInfoMap, BufferedWriter writer, String assignedMods) throws IOException {
        // get all peptide for this index
        String[] peptides = peptide.split(";");
        // deduplicate and group extra info for each peptide in the same index
        Set<String> spectrums = new TreeSet<>();
        Set<String> proteins = new TreeSet<>();
        Set<String> entryNames = new TreeSet<>();
        Set<String> proteinDescs = new TreeSet<>();
        Set<String> mappedGenes = new TreeSet<>();
        Set<String> mappedProteins = new TreeSet<>();

        for (String peptide : peptides) {
            ExtraPsmInfo extraPsmInfo = extraPsmInfoMap.get(peptide.toUpperCase()); // upper case in psm file
            spectrums.addAll(extraPsmInfo.getSpectrums());
            proteins.add(extraPsmInfo.getProtein());
            entryNames.add(extraPsmInfo.getEntryName());
            proteinDescs.add(extraPsmInfo.getProteinDesc());
            if (!extraPsmInfo.getMappedGenes().isEmpty()) {
                mappedGenes.add(extraPsmInfo.getMappedGenes());
            }
            if (!extraPsmInfo.getMappedProteins().isEmpty()) {
                mappedProteins.add(extraPsmInfo.getMappedProteins());
            }
        }

        // generate the final report info
        String result = index + "\t" +
                gene + "\t" +
                proteinID + "\t" +
                peptide + "\t";

        if (!assignedMods.isEmpty()) {
            result += assignedMods + "\t";
        }

        result += sequenceWindow + "\t" +
                start + "\t" +
                end + "\t" +
                maxPepProb + "\t" +
                spectrums.size() + "\t" + // Just for clarity in case of large number of spectrum
                String.join(";", proteins) + "\t" +
                String.join(";", entryNames) + "\t" +
                String.join(";", proteinDescs) + "\t" +
                String.join(";", mappedGenes) + "\t" +
                String.join(";", mappedProteins);

        writer.write(result);
    }

    public void propagateExtraInfoProtein(Map<String, ExtraPsmInfo> proteinMap, BufferedWriter writer) throws IOException {
        // protein level report index is protein ID
        String proteinId = index;
        ExtraPsmInfo extraPsmInfo = proteinMap.get(proteinId);
        String result = index + "\t" +
                numPSM + "\t" +
                gene + "\t" +
                maxPepProb + "\t" +
                extraPsmInfo.getProtein() + "\t" +
                proteinId + "\t" +
                extraPsmInfo.getEntryName() + "\t" +
                extraPsmInfo.getProteinDesc() + "\t" +
                extraPsmInfo.getOrganism() + "\t" +
                extraPsmInfo.getIndistinguishableProteins();

        writer.write(result);
    }
}
