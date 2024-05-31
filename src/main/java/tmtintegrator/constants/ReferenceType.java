package tmtintegrator.constants;

/**
 * Reference channel type
 *
 * @author rogerli on 05/2024
 */
public enum ReferenceType {
    RAW_ABUNDANCE(-2), // export raw abundance
    NONE(-1), // don't add reference
    SUMMATION(0), // use summation as reference
    AVERAGE(1), // use average as reference
    MEDIAN(2); // use median as reference

    private final int value;

    ReferenceType(int value) {
        this.value = value;
    }

    public int getValue() {
        return value;
    }

    public static ReferenceType fromValue(int value) {
        for (ReferenceType referenceType : ReferenceType.values()) {
            if (referenceType.getValue() == value) {
                return referenceType;
            }
        }
        throw new IllegalArgumentException("Invalid referenceType value: " + value);
    }
}
