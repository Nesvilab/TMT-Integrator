/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
