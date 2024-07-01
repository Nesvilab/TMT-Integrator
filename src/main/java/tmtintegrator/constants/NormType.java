package tmtintegrator.constants;

/**
 * Normalization types for PSM data.
 */
public enum NormType {
    ALL_NORM(-1),
    NONE(0),
    MC(1), // Median centering
    GN(2), // median centering + variance scaling
    SL_IRS(3); // sample loading and internal reference scaling

    private final int value;

    NormType(int value) {
        this.value = value;
    }

    public int getValue() {
        return value;
    }

    public static NormType fromValue(int value) {
        for (NormType normType : NormType.values()) {
            if (normType.getValue() == value) {
                return normType;
            }
        }
        throw new IllegalArgumentException("Invalid protNorm value: " + value);
    }
}
