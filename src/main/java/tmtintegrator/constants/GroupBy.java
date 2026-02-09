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
 * Level of data summarization.
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
