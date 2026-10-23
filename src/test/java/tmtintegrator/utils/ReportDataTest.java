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

package tmtintegrator.utils;

import org.junit.jupiter.api.Test;

import tmtintegrator.constants.GroupBy;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.ReportInfo;

import java.io.BufferedWriter;
import java.io.StringWriter;

import static org.junit.jupiter.api.Assertions.assertEquals;

class ReportDataTest {

    // group key info part: Index, Gene, ProteinID, Peptide, SequenceWindow, Start, End, MaxPepProb
    private static final String KEY_PREFIX = "P12345%n[304.2071]PEPTIDE\tGENE1\tP12345\tPEPTIDE\tABCDEFG.PEPTIDE.HIJKLMN\t10\t16\t0.95";

    private static String writeReportInfo(String groupKey, GroupBy groupBy) throws Exception {
        Index index = new Index();
        index.spectrumIndex = 0;
        index.pepcIndex = 1;
        index.proteincIndex = 2;
        index.genecIndex = 3;
        index.proteinIDcIndex = 4;

        ReportData reportData = new ReportData();
        reportData.updateExtraPsmInfo(new String[]{"run1.00001.00001.2", "PEPTIDE", "sp|P12345|PROT_HUMAN", "GENE1", "P12345"}, index);

        StringWriter stringWriter = new StringWriter();
        try (BufferedWriter writer = new BufferedWriter(stringWriter)) {
            reportData.writeReportInfo(writer, groupKey, groupBy);
        }
        return stringWriter.toString();
    }

    @Test
    void modifiedPeptideWithAssignedModsKeepsColumns() throws Exception {
        String[] columns = writeReportInfo(KEY_PREFIX + "\tN-term(304.2071)", GroupBy.MODIFIED_PEPTIDE).split("\t", -1);
        assertEquals(ReportInfo.HEADER_MODIFIED_PEPTIDE.split("\t").length, columns.length);
        assertEquals("N-term(304.2071)", columns[4]);
        assertEquals("0.95", columns[8]);
    }

    @Test
    void modifiedPeptideWithEmptyAssignedModsKeepsColumns() throws Exception {
        // An empty Assigned Modifications is the last field of the key and must not shift the other columns
        String[] columns = writeReportInfo(KEY_PREFIX + "\t", GroupBy.MODIFIED_PEPTIDE).split("\t", -1);
        assertEquals(ReportInfo.HEADER_MODIFIED_PEPTIDE.split("\t").length, columns.length);
        assertEquals("", columns[4]);
        assertEquals("ABCDEFG.PEPTIDE.HIJKLMN", columns[5]);
        assertEquals("0.95", columns[8]);
    }

    @Test
    void peptideReportHasNoAssignedModsColumn() throws Exception {
        String[] columns = writeReportInfo(KEY_PREFIX, GroupBy.PEPTIDE).split("\t", -1);
        assertEquals(ReportInfo.HEADER.split("\t").length, columns.length);
        assertEquals("ABCDEFG.PEPTIDE.HIJKLMN", columns[4]);
        assertEquals("0.95", columns[7]);
    }
}
