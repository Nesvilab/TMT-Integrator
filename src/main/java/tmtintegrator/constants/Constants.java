package tmtintegrator.constants;

import java.util.regex.Pattern;

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

    // regex to match mod tag with format: A(123.456), group 1: A, group 2: 123.456
    public static final Pattern MOD_TAG_PATTERN = Pattern.compile("([A-Z])\\((\\d+\\.\\d+)\\)");
    // regex to match glyco mod with format: N(203.079), group 1: N, group 2: 203.079
    public static final Pattern GLYCO_MOD_PATTERN = Pattern.compile("([A-Z])\\((\\d+\\.\\d+)\\)");
    public static final int AA_GROUP = 1;
    public static final int MASS_GROUP = 2;
    public static final Pattern KEY_PATTERN = Pattern.compile("([^%]+)%([ncA-Z])(\\d+)"); // Example: A0A096%S13
    public static final int BIN_NUM = 10; // number of bins for retention time
    public static final int PSM_NUM_THRESHOLD = 4; // threshold for outlier removal
}
