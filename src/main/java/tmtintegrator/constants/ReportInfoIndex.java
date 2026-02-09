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
 * Column index for report info columns.
 */
public class ReportInfoIndex {

    private ReportInfoIndex() {
        throw new AssertionError("ReportInfoIndex is a utility class");
    }

    /**
     * Column format for report info part for peptide and site level:
     * <Index> <Gene> <ProteinID> <Peptide> <SequenceWindow> <Start> <End> <MaxPepProb>
     */
    public static final int INDEX = 0;
    public static final int GENE = 1;
    public static final int PROTEIN_ID = 2;
    public static final int PEPTIDE = 3;
    public static final int SEQUENCE_WINDOW = 4;
    public static final int START = 5;
    public static final int END = 6;
    public static final int MAX_PEP_PROB = 7;
    public static final int ASSIGNED_MODIFICATIONS = 8;

    /**
     * Column format for report info part for protein level:
     * <Index> <NumberPSM> <Gene> <MaxPepProb>
     */
    public static final int NUMBER_PSM = 1;
    public static final int GENE_PROTEIN = 2;
    public static final int MAX_PEP_PROB_PROTEIN = 3;
}
