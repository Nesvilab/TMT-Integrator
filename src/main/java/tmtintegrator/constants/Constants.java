package tmtintegrator.constants;

/**
 * Constants class for common arithmetic operations, and psm processing.
 *
 * @author rogerli on 05/2024
 */
public final class Constants {

    private Constants() {
        // private constructor to prevent instantiation
        throw new AssertionError("The Constants class cannot be instantiated");
    }

    public static final int BIN_NUM = 10; // number of bins for retention time
    public static final int PSM_NUM_THRESHOLD = 4; // threshold for outlier removal
}
