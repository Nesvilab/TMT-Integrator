package tmtintegrator.constants;

/**
 * Level of data summarization.
 *
 * @author rogerli on 05/2024
 */
public enum GroupBy {

    GENE(0), // PSM aggregation to the gene level
    PROTEIN_ID(1), // protein level
    PEPTIDE(2), // peptide sequence level
    MULTI_PHOSPHO_SITE(3), // multiple PTM sites
    SINGLE_PHOSPHO_SITE(4), // single PTM site
    MULTI_MASS_GLYCO(5), // multi-mass for glycolysation
    MODIFIED_PEPTIDE(6); // modified peptide report

    private final int value;

    GroupBy(int value) {
        this.value = value;
    }

    public int getValue() {
        return value;
    }

    public static GroupBy fromValue(int value) {
        for (GroupBy groupBy : GroupBy.values()) {
            if (groupBy.getValue() == value) {
                return groupBy;
            }
        }
        throw new IllegalArgumentException("Invalid groupBy value: " + value);
    }
}
