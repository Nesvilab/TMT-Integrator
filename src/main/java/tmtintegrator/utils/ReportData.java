package tmtintegrator.utils;

import tmtintegrator.constants.GroupBy;
import tmtintegrator.pojo.ExtraPsmInfo;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.ProteinIndex;
import tmtintegrator.pojo.ReportInfo;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Maintain all necessary data for the report
 */
public class ReportData {
    private final Map<String, ExtraPsmInfo> extraPsmInfoMap; // <peptide sequence, extra psm info>
    private final Map<String, ExtraPsmInfo> proteinMap; // <protein, extra psm info>

    public ReportData() {
        extraPsmInfoMap = new HashMap<>();
        proteinMap = new HashMap<>();
    }

    /**
     * Information that will be used in the report from psm.tsv
     * - protein
     * - protein ID
     * - Entry Name
     * - Gene
     * - Protein Description
     * - Mapped Genes
     * - Mapped Proteins
     *
     * @param fields one psm line
     * @param index  index of the fields
     */
    public void updateExtraPsmInfo(String[] fields, Index index) {
        updateExtraPsmInfoMap(fields, index);
        updateProteinMap(fields, index);
    }

    public void updateExtraProteinInfo(String[] fields, ProteinIndex index) {
        String proteinId = fields[index.proteinIDIdx];
        ExtraPsmInfo extraPsmInfo = proteinMap.computeIfAbsent(proteinId, k -> new ExtraPsmInfo());

        extraPsmInfo.setOrganism(fields[index.organismIdx].trim());

        // indistinguishable Proteins may be empty in protein.tsv
        if (fields.length > index.indistinguishableProteinsIdx) {
            extraPsmInfo.setIndistinguishableProteins(fields[index.indistinguishableProteinsIdx]);
        } else {
            extraPsmInfo.setIndistinguishableProteins("");
        }
    }

    public void writeReportInfo(BufferedWriter writer, String infoString, GroupBy groupBy) throws IOException {
        String[] parts = infoString.split("\t");
        // extract report info
        ReportInfo reportInfo = new ReportInfo(parts, groupBy);

        // propagate extra info from original psm file
        if (groupBy == GroupBy.PROTEIN_ID) {
            reportInfo.propagateExtraInfoProtein(proteinMap, writer);
        } else {
            reportInfo.propagateExtraInfo(extraPsmInfoMap, writer);
        }
    }

    // region helper methods
    private void updateExtraPsmInfoMap(String[] fields, Index index) {
        String peptide = fields[index.pepcIndex];
        ExtraPsmInfo extraPsmInfo = extraPsmInfoMap.computeIfAbsent(peptide, k -> new ExtraPsmInfo());

        extraPsmInfo.addSpectrum(fields[index.spectrumIndex]);

        if (extraPsmInfo.getSpectrums().size() == 1) {
            // Only update the first time
            extraPsmInfo.setProtein(fields[index.proteincIndex]);
            extraPsmInfo.setProteinId(fields[index.proteinIDcIndex]);
            extraPsmInfo.setEntryName(fields[index.entryNameIndex]);
            extraPsmInfo.setGene(fields[index.genecIndex]);
            extraPsmInfo.setProteinDesc(fields[index.proteinDescIndex]);
            extraPsmInfo.setMappedGenes(fields[index.mapGeneIndex]);
            extraPsmInfo.setMappedProteins(fields[index.mappedProteinsIndex]);
        }
    }

    private void updateProteinMap(String[] fields, Index index) {
        String proteinId = fields[index.proteinIDcIndex];
        ExtraPsmInfo extraPsmInfo = proteinMap.computeIfAbsent(proteinId, k -> new ExtraPsmInfo());

        // Assume these fields are unique for each protein ID
        extraPsmInfo.setProtein(fields[index.proteincIndex]);
        extraPsmInfo.setEntryName(fields[index.entryNameIndex]);
        extraPsmInfo.setGene(fields[index.genecIndex]);
        extraPsmInfo.setProteinDesc(fields[index.proteinDescIndex]);
    }
    // endregion
}
