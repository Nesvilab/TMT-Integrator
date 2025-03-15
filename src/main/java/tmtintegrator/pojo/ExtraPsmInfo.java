package tmtintegrator.pojo;

import java.util.*;

public class ExtraPsmInfo {
    // region psm.tsv columns
    private final Set<String> spectrums;
    private String protein;
    private String entryName;
    private String gene;
    private String proteinDesc;
    private String mappedGenes;
    private String mappedProteins;
    // endregion

    // region protein.tsv columns
    private String organism;
    private String indistinguishableProteins;
    // endregion

    public ExtraPsmInfo() {
        this.spectrums = new TreeSet<>();
    }

    // region getters and setters
    public Set<String> getSpectrums() {
        return spectrums;
    }

    public void addSpectrum(String spectrum) {
        this.spectrums.add(spectrum);
    }

    public String getProtein() {
        return protein;
    }

    public void setProtein(String protein) {
        this.protein = protein;
    }

    public String getEntryName() {
        return entryName;
    }

    public void setEntryName(String entryName) {
        this.entryName = entryName;
    }

    public String getGene() {
        return gene;
    }

    public void setGene(String gene) {
        this.gene = gene;
    }

    public String getProteinDesc() {
        return proteinDesc;
    }

    public void setProteinDesc(String proteinDesc) {
        this.proteinDesc = proteinDesc;
    }

    public String getMappedGenes() {
        return mappedGenes;
    }

    public void setMappedGenes(String mappedGenes) {
        this.mappedGenes = mappedGenes;
    }

    public String getMappedProteins() {
        return mappedProteins;
    }

    public void setMappedProteins(String mappedProteins) {
        this.mappedProteins = mappedProteins;
    }

    public String getOrganism() {
        return organism;
    }

    public void setOrganism(String organism) {
        this.organism = organism;
    }

    public String getIndistinguishableProteins() {
        return indistinguishableProteins;
    }

    public void setIndistinguishableProteins(String indistinguishableProteins) {
        this.indistinguishableProteins = indistinguishableProteins;
    }

    // endregion
}
