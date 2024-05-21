package tmtintegrator.pojo;

public enum GroupBy {
    GENE(0),
    PROTEIN_ID(1),
    PEPTIDE(2),
    MULTI_PHOSPHO_SITE(3),
    SINGLE_PHOSPHO_SITE(4),
    MULTI_MASS_GLYCO(5);

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
