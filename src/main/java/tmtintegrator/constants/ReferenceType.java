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
 * Reference channel type
 */
public enum ReferenceType {
    RAW_ABUNDANCE(-2), // export raw abundance
    REAL(-1), // use real reference
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
