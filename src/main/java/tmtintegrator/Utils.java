package tmtintegrator;

import java.util.Collections;
import java.util.List;

public final class Utils {

    private Utils() {
        // private constructor to prevent instantiation
        throw new AssertionError("The Utils class cannot be instantiated");
    }

    public static double tryParseDouble(String input) {
        try {
            return Double.parseDouble(input);
        } catch (NumberFormatException e) {
            return Double.NaN; // TODO: need test
        }
    }

    public static double takeMedian(List<Double> numbers) {
        if (numbers == null || numbers.isEmpty()) {
            return Double.NaN;
        }

        Collections.sort(numbers);
        int size = numbers.size();
        int middleIndex = size / 2;

        if (size % 2 == 1) { // Odd number of elements
            return numbers.get(middleIndex);
        } else if (size == 2) { // Special case for two elements only
            return (numbers.get(0) + numbers.get(1)) / 2.0;
        } else if (size > 2) { // Even number of elements
            return (numbers.get(middleIndex - 1) + numbers.get(middleIndex)) / 2.0;
        }

        return Double.NaN; // Fallback for any other unexpected cases
    }

}
