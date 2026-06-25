package tmtintegrator.pojo.psm;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Nested;

import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.Subplex;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class PsmRecordTest {

    // region reflection helpers

    @SuppressWarnings("unchecked")
    private static <T> T getField(PsmRecord record, String fieldName) throws Exception {
        Field field = PsmRecord.class.getDeclaredField(fieldName);
        field.setAccessible(true);
        return (T) field.get(record);
    }

    // endregion

    @Nested
    class StripBracketModsTest {

        @Test
        void removesSimpleBracketMod() {
            assertEquals("PEPTIDE", PsmRecord.stripBracketMods("PEP[79.9663]TIDE"));
        }

        @Test
        void removesNTermMod() {
            assertEquals("PEPTIDE", PsmRecord.stripBracketMods("n[42.0106]PEPTIDE"));
        }

        @Test
        void removesCTermMod() {
            assertEquals("PEPTIDE", PsmRecord.stripBracketMods("PEPTIDEc[57.0215]"));
        }

        @Test
        void removesMultipleMods() {
            assertEquals("PEPTIDE", PsmRecord.stripBracketMods("n[42.0106]PEP[79.9663]TIDEc[57.0215]"));
        }

        @Test
        void removesMultipleInternalMods() {
            assertEquals("ABCDEF", PsmRecord.stripBracketMods("A[100.0]BC[200.0]D[300.0]EF"));
        }

        @Test
        void handlesNoMods() {
            assertEquals("PEPTIDE", PsmRecord.stripBracketMods("PEPTIDE"));
        }

        @Test
        void handlesNull() {
            assertNull(PsmRecord.stripBracketMods(null));
        }

        @Test
        void handlesEmpty() {
            assertEquals("", PsmRecord.stripBracketMods(""));
        }

        @Test
        void handlesOnlyNTermMod() {
            assertEquals("A", PsmRecord.stripBracketMods("n[229.1629]A"));
        }

        @Test
        void handlesComplexMasses() {
            assertEquals("SAMPLER", PsmRecord.stripBracketMods("S[167.0]A[304.2072]M[147.0354]PLER"));
        }
    }

    @Nested
    class ParseChargeFromSpectrumTest {

        @Test
        void parsesTypicalSpectrum() {
            assertEquals(3, PsmRecord.parseChargeFromSpectrum("file.00123.00123.3"));
        }

        @Test
        void parsesCharge2() {
            assertEquals(2, PsmRecord.parseChargeFromSpectrum("run1.00456.00456.2"));
        }

        @Test
        void parsesHighCharge() {
            assertEquals(5, PsmRecord.parseChargeFromSpectrum("experiment.99999.99999.5"));
        }

        @Test
        void throwsForNull() {
            assertThrows(IllegalArgumentException.class, () -> PsmRecord.parseChargeFromSpectrum(null));
        }

        @Test
        void throwsForEmpty() {
            assertThrows(IllegalArgumentException.class, () -> PsmRecord.parseChargeFromSpectrum(""));
        }

        @Test
        void throwsForNoDot() {
            assertThrows(IllegalArgumentException.class, () -> PsmRecord.parseChargeFromSpectrum("nodotshere"));
        }

        @Test
        void throwsWhenDotIsLastChar() {
            assertThrows(IllegalArgumentException.class, () -> PsmRecord.parseChargeFromSpectrum("file.123."));
        }

        @Test
        void throwsForNonNumericCharge() {
            assertThrows(IllegalArgumentException.class, () -> PsmRecord.parseChargeFromSpectrum("file.123.abc"));
        }

        @Test
        void parsesSingleDotSpectrum() {
            assertEquals(4, PsmRecord.parseChargeFromSpectrum("something.4"));
        }
    }

    @Nested
    class GenerateModifiedSequenceTest {

        private String invokeGenerateModifiedSequence(String sequence, String assignMods) throws Exception {
            Parameters params = new Parameters();
            params.subplexes.add(new Subplex("ref", 6));
            Index index = new Index();
            PsmRecord record = new PsmRecord(params, index);
            Method method = PsmRecord.class.getDeclaredMethod("generateModifiedSequence", String.class, String.class);
            method.setAccessible(true);
            return (String) method.invoke(record, sequence, assignMods);
        }

        @Test
        void noModsReturnsSequenceUnchanged() throws Exception {
            assertEquals("PEPTIDE", invokeGenerateModifiedSequence("PEPTIDE", ""));
        }

        @Test
        void singleInternalMod() throws Exception {
            String result = invokeGenerateModifiedSequence("PEPTIDE", "3P(79.9663)");
            assertEquals("PEP[79.9663]TIDE", result);
        }

        @Test
        void nTermMod() throws Exception {
            String result = invokeGenerateModifiedSequence("PEPTIDE", "N-term(42.0106)");
            assertEquals("n[42.0106]PEPTIDE", result);
        }

        @Test
        void cTermMod() throws Exception {
            // This tests the bug fix: c-term was previously using modMap.get(0) instead of modMap.get(seqLen+1)
            String result = invokeGenerateModifiedSequence("PEPTIDE", "C-term(57.0215)");
            assertEquals("PEPTIDEc[57.0215]", result);
        }

        @Test
        void bothTerminalMods() throws Exception {
            String result = invokeGenerateModifiedSequence("PEP", "N-term(42.0106), C-term(57.0215)");
            assertEquals("n[42.0106]PEPc[57.0215]", result);
        }

        @Test
        void multipleMods() throws Exception {
            String result = invokeGenerateModifiedSequence("PEPTIDE", "N-term(42.0106), 3P(79.9663), C-term(57.0215)");
            assertEquals("n[42.0106]PEP[79.9663]TIDEc[57.0215]", result);
        }

        @Test
        void multipleInternalMods() throws Exception {
            String result = invokeGenerateModifiedSequence("ABCDEF", "1A(100.0), 3C(200.0), 5E(300.0)");
            assertEquals("A[100.0]BC[200.0]DE[300.0]F", result);
        }
    }

    // region shared test helpers

    private static Parameters buildStandardParams() {
        Parameters params = new Parameters();
        params.isHyperplexing = false;
        params.addIsobaricFilter = false;
        params.add_Ref = 0; // summation virtual reference
        params.glycoQval = -1;
        params.allow_unlabeled = true;
        params.subplexes.add(new Subplex("ref_126", 6));
        return params;
    }

    private static Index buildStandardIndex() {
        Index index = new Index();
        index.spectrumIndex = 0;
        index.peptideIndex = 1;
        index.extpepIndex = 2;
        index.chargeIndex = 3;
        index.rtIndex = 4;
        index.pepMassIndex = 5;
        index.pepProbcIndex = 6;
        index.numEnzyTermi = 7;
        index.protsIndex = 8;
        index.proteIndex = 9;
        index.labelMs1IntIndices[0] = 10;
        index.assignedModcIndex = 11;
        index.purityIndex = 12;
        index.isUniquecIndex = 13;
        index.proteincIndex = 14;
        index.proteinIDcIndex = 15;
        index.genecIndex = 16;

        Index.SubplexIndex si = new Index.SubplexIndex();
        si.abnIndex = 17;
        si.refIndex = 17;
        index.subplexIndices.add(si);

        return index;
    }

    private static String[] buildStandardFields() {
        return new String[]{
            "run1.00100.00100.2",   // 0: spectrum
            "PEPTIDE",              // 1: peptide
            "K.PEPTIDE.R",          // 2: extended peptide
            "2",                    // 3: charge
            "30.5",                 // 4: retention time
            "800.4",               // 5: calculated peptide mass
            "0.99",                 // 6: probability
            "2",                    // 7: num enzymatic termini
            "10",                   // 8: protein start
            "17",                   // 9: protein end
            "50000",               // 10: MS1 intensity
            "",                     // 11: assigned modifications (empty)
            "0.95",                 // 12: purity
            "true",                 // 13: is unique
            "sp|P12345|PROT_HUMAN", // 14: protein
            "P12345",              // 15: protein ID
            "GENE1",               // 16: gene
            "1000",                // 17: channel 1 (ref)
            "2000",                // 18: channel 2
            "3000",                // 19: channel 3
            "4000",                // 20: channel 4
            "5000",                // 21: channel 5
            "6000",                // 22: channel 6
        };
    }

    private static PsmRecord parseStandardLine(String[] fields) {
        Parameters params = buildStandardParams();
        Index index = buildStandardIndex();
        PsmRecord record = new PsmRecord(params, index);
        record.parsePsmRecord(String.join("\t", fields), new ArrayList<>(), new ArrayList<>());
        return record;
    }

    private static Parameters buildHyperplexingParams() {
        Parameters params = new Parameters();
        params.isHyperplexing = true;
        params.addIsobaricFilter = false;
        params.add_Ref = 0; // summation virtual reference
        params.glycoQval = -1;
        params.allow_unlabeled = true;
        params.subplexes.add(new Subplex("ref_126", 3));
        params.subplexes.add(new Subplex("ref_126", 3));
        params.hyperplexLabelNames.add("Light");
        params.hyperplexLabelNames.add("Heavy");
        params.hyperplexLabelSubplexCounts.add(1);
        params.hyperplexLabelSubplexCounts.add(1);
        return params;
    }

    private static Index buildHyperplexingIndex() {
        Index index = new Index();
        index.labelFreeModifiedPeptideIndex = 0;
        index.labelSpectrumIndices[0] = 1; // Light Spectrum
        index.labelSpectrumIndices[1] = -1; // Medium (not present)
        index.labelSpectrumIndices[2] = 2; // Heavy Spectrum
        index.labelAssignedModIndices[0] = 3;
        index.labelMs1IntIndices[0] = 4; // Light Intensity
        index.labelMs1IntIndices[2] = 5; // Heavy Intensity (index 2 = Heavy)
        index.pepcIndex = 0;
        index.spectrumIndex = 1;
        index.rtIndex = 6;
        index.pepMassIndex = 7;
        index.pepProbcIndex = 8;
        index.purityIndex = 9;
        index.isUniquecIndex = 10;
        index.proteincIndex = 11;
        index.proteinIDcIndex = 12;
        index.genecIndex = 13;
        index.observedMzIndex = 14;

        Index.SubplexIndex lightSi = new Index.SubplexIndex();
        lightSi.abnIndex = 15;
        lightSi.refIndex = 15;
        index.subplexIndices.add(lightSi);

        Index.SubplexIndex heavySi = new Index.SubplexIndex();
        heavySi.abnIndex = 18;
        heavySi.refIndex = 18;
        index.subplexIndices.add(heavySi);

        return index;
    }

    private static String[] buildHyperplexingFields() {
        return new String[]{
            "n[229.16]PEPTIDE",      // 0: Label Free Modified Peptide
            "run1.00100.00100.3",     // 1: Light Spectrum
            "run1.00100.00100.3",     // 2: Heavy Spectrum
            "1K(229.1629)",           // 3: Light Assigned Modifications
            "20000",                  // 4: Light Intensity
            "30000",                  // 5: Heavy Intensity
            "25.5",                   // 6: RT
            "850.4",                  // 7: Calc Pep Mass
            "0.98",                   // 8: Probability
            "0.92",                   // 9: Purity
            "true",                   // 10: Is Unique
            "sp|P12345|PROT_HUMAN",   // 11: Protein
            "P12345",                 // 12: Protein ID
            "GENE1",                  // 13: Gene
            "425.2",                  // 14: Observed M/Z
            "1000", "2000", "3000",   // 15-17: Light channels
            "4000", "5000", "6000",   // 18-20: Heavy channels
        };
    }

    private static PsmRecord parseHyperplexingLine(String[] fields) {
        Parameters params = buildHyperplexingParams();
        Index index = buildHyperplexingIndex();
        PsmRecord record = new PsmRecord(params, index);
        record.parsePsmRecord(String.join("\t", fields), new ArrayList<>(), new ArrayList<>());
        return record;
    }

    // endregion

    @Nested
    class ParseStandardFieldsTest {

        @Test
        void parsesAllPublicFields() {
            PsmRecord record = parseStandardLine(buildStandardFields());

            assertEquals("PEPTIDE", record.getPeptide());
            assertEquals("PEPTIDE", record.getModifiedPeptide());
            assertEquals("K.PEPTIDE.R", record.getExtendedPeptide());
            assertEquals(0.99, record.getProbability(), 1e-6);
            assertEquals(30.5, record.getRetention(), 1e-6);
            assertEquals("GENE1", record.getGene());
            assertEquals("P12345", record.getProteinId());
            assertEquals("", record.getAssignedModifications());
            assertEquals(50000.0, record.getMs1Intensity(), 1e-6);
            assertEquals(21000.0, record.getSumTmtIntensity(), 1e-6);
            assertEquals(9, record.getPepsIndex()); // proteinStart(10) - 1
        }

        @Test
        void parsesAllPrivateFields() throws Exception {
            PsmRecord record = parseStandardLine(buildStandardFields());

            assertEquals("run1.00100.00100.2", getField(record, "spectrum"));
            assertEquals(2, (int) getField(record, "charge"));
            assertEquals(800.4, (double) getField(record, "calculatedPepMass"), 1e-6);
            assertEquals(2, (int) getField(record, "numEnzymaticTermini"));
            assertEquals(10, (int) getField(record, "proteinStart"));
            assertEquals(17, (int) getField(record, "proteinEnd"));
            assertEquals(0.95, (double) getField(record, "purity"), 1e-6);
            assertTrue((boolean) getField(record, "isUnique"));
            assertEquals("sp|P12345|PROT_HUMAN", getField(record, "protein"));
            assertEquals("", getField(record, "ptmLocalization"));
            assertEquals("", getField(record, "mappedGenes"));
        }
    }

    @Nested
    class ParseHyperplexingFieldsTest {

        @Test
        void parsesAllPublicFields() {
            PsmRecord record = parseHyperplexingLine(buildHyperplexingFields());

            assertEquals("PEPTIDE", record.getPeptide());
            assertEquals("n[229.16]PEPTIDE", record.getModifiedPeptide());
            assertEquals("1K(229.1629)", record.getAssignedModifications());
            assertEquals("P12345", record.getProteinId());
            assertEquals("GENE1", record.getGene());
            assertEquals(0.98, record.getProbability(), 1e-6);
            assertEquals(25.5, record.getRetention(), 1e-6);
            assertEquals(50000.0, record.getMs1Intensity(), 1e-6);
            assertEquals(21000.0, record.getSumTmtIntensity(), 1e-6);
            assertEquals(-1, record.getPepsIndex());
        }

        @Test
        void parsesAllPrivateFields() throws Exception {
            PsmRecord record = parseHyperplexingLine(buildHyperplexingFields());

            assertEquals("run1.00100.00100.3", getField(record, "spectrum"));
            assertEquals(3, (int) getField(record, "charge"));
            assertEquals(850.4, (double) getField(record, "calculatedPepMass"), 1e-6);
            assertEquals(0.92, (double) getField(record, "purity"), 1e-6);
            assertTrue((boolean) getField(record, "isUnique"));
            assertEquals("sp|P12345|PROT_HUMAN", getField(record, "protein"));
            assertEquals("", getField(record, "ptmLocalization"));
            assertEquals("", getField(record, "mappedGenes"));
            assertEquals(2, (int) getField(record, "numEnzymaticTermini"));
            assertEquals(0, (int) getField(record, "proteinStart"));
            assertEquals(0, (int) getField(record, "proteinEnd"));
            assertEquals("", getField(record, "extendedPeptide"));
        }

        @Test
        void usesFirstNonEmptyLabelSpectrum() throws Exception {
            String[] fields = buildHyperplexingFields();
            fields[1] = ""; // Light Spectrum is empty

            PsmRecord record = parseHyperplexingLine(fields);

            assertEquals("run1.00100.00100.3", getField(record, "spectrum"));
        }
    }

    @Nested
    class ParseChannelsTest {

        @Test
        void parsesStandardChannels() throws Exception {
            PsmRecord record = parseStandardLine(buildStandardFields());

            List<List<Double>> subplexChannels = getField(record, "subplexChannels");
            assertEquals(1, subplexChannels.size());
            assertEquals(Arrays.asList(1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0), subplexChannels.get(0));
        }

        @Test
        void skipsNaChannelsInStandard() throws Exception {
            Parameters params = buildStandardParams();
            Index index = buildStandardIndex();
            PsmRecord record = new PsmRecord(params, index);
            List<Integer> naChannels = new ArrayList<>(List.of(19)); // skip channel at col 19

            record.parsePsmRecord(String.join("\t", buildStandardFields()), new ArrayList<>(), naChannels);

            List<List<Double>> subplexChannels = getField(record, "subplexChannels");
            assertEquals(5, subplexChannels.get(0).size());
            assertEquals(18000.0, record.getSumTmtIntensity(), 1e-6); // 21000 - 3000
        }

        @Test
        void parsesHyperplexingChannelsPerSubplex() throws Exception {
            PsmRecord record = parseHyperplexingLine(buildHyperplexingFields());

            List<List<Double>> subplexChannels = getField(record, "subplexChannels");
            assertEquals(2, subplexChannels.size());
            assertEquals(Arrays.asList(1000.0, 2000.0, 3000.0), subplexChannels.get(0));
            assertEquals(Arrays.asList(4000.0, 5000.0, 6000.0), subplexChannels.get(1));
        }

        @Test
        void handlesEmptyChannelValues() throws Exception {
            String[] fields = buildHyperplexingFields();
            fields[16] = ""; // one Light channel is empty

            PsmRecord record = parseHyperplexingLine(fields);

            List<List<Double>> subplexChannels = getField(record, "subplexChannels");
            assertEquals(Arrays.asList(1000.0, 0.0, 3000.0), subplexChannels.get(0));
        }
    }

    @Nested
    class ApportionMs1IntensitiesTest {

        @Test
        void standardSingleSubplexGetsFullMs1() throws Exception {
            PsmRecord record = parseStandardLine(buildStandardFields());

            double[] subplexMs1 = getField(record, "subplexMs1Intensities");
            assertEquals(1, subplexMs1.length);
            assertEquals(50000.0, subplexMs1[0], 1e-6);
        }

        @Test
        void standardMultiSubplexApportionsByTmtShare() throws Exception {
            // Two subplexes of 3 channels each
            Parameters params = new Parameters();
            params.isHyperplexing = false;
            params.addIsobaricFilter = false;
            params.add_Ref = 0;
            params.glycoQval = -1;
            params.allow_unlabeled = true;
            params.subplexes.add(new Subplex("ref_126", 3));
            params.subplexes.add(new Subplex("ref_131", 3));

            Index index = buildStandardIndex();
            // Second subplex channels at col 20-22
            Index.SubplexIndex si2 = new Index.SubplexIndex();
            si2.abnIndex = 20;
            si2.refIndex = 20;
            index.subplexIndices.add(si2);

            PsmRecord record = new PsmRecord(params, index);
            record.parsePsmRecord(String.join("\t", buildStandardFields()), new ArrayList<>(), new ArrayList<>());

            double[] subplexMs1 = getField(record, "subplexMs1Intensities");
            assertEquals(2, subplexMs1.length);
            // Plex1: channels 1000+2000+3000=6000, Plex2: channels 4000+5000+6000=15000
            // Total TMT = 21000, MS1 = 50000
            assertEquals(50000.0 * 6000.0 / 21000.0, subplexMs1[0], 1e-6);
            assertEquals(50000.0 * 15000.0 / 21000.0, subplexMs1[1], 1e-6);
        }

        @Test
        void hyperplexingUsesPerLabelMs1Directly() throws Exception {
            PsmRecord record = parseHyperplexingLine(buildHyperplexingFields());

            double[] subplexMs1 = getField(record, "subplexMs1Intensities");
            assertEquals(2, subplexMs1.length);
            assertEquals(20000.0, subplexMs1[0], 1e-6); // Light
            assertEquals(30000.0, subplexMs1[1], 1e-6); // Heavy
            assertEquals(50000.0, record.getMs1Intensity(), 1e-6); // total
        }

        @Test
        void hyperplexingHandlesEmptyMs1() throws Exception {
            String[] fields = buildHyperplexingFields();
            fields[4] = ""; // Light Intensity is empty

            PsmRecord record = parseHyperplexingLine(fields);

            double[] subplexMs1 = getField(record, "subplexMs1Intensities");
            assertEquals(0.0, subplexMs1[0], 1e-6);
            assertEquals(30000.0, subplexMs1[1], 1e-6);
            assertEquals(30000.0, record.getMs1Intensity(), 1e-6);
        }

        @Test
        void hyperplexingMultiSubplexPerLabelApportionsByTmtShare() throws Exception {
            // Two subplexes per label (e.g. TMT-35 with 2+1 channels for test simplicity)
            Parameters params = new Parameters();
            params.isHyperplexing = true;
            params.addIsobaricFilter = false;
            params.add_Ref = 0;
            params.glycoQval = -1;
            params.allow_unlabeled = true;
            // Light: 2 subplexes (2 channels + 1 channel)
            params.subplexes.add(new Subplex("ref_126", 2));
            params.subplexes.add(new Subplex("ref_128", 1));
            // Heavy: 2 subplexes (2 channels + 1 channel)
            params.subplexes.add(new Subplex("ref_126", 2));
            params.subplexes.add(new Subplex("ref_128", 1));
            params.hyperplexLabelNames.add("Light");
            params.hyperplexLabelNames.add("Heavy");
            params.hyperplexLabelSubplexCounts.add(2);
            params.hyperplexLabelSubplexCounts.add(2);

            Index index = new Index();
            index.labelFreeModifiedPeptideIndex = 0;
            index.labelSpectrumIndices[0] = 1;
            index.labelSpectrumIndices[2] = 2;
            index.labelAssignedModIndices[0] = 3;
            index.labelMs1IntIndices[0] = 4; // Light Intensity = 20000
            index.labelMs1IntIndices[2] = 5; // Heavy Intensity = 30000
            index.pepcIndex = 0;
            index.spectrumIndex = 1;
            index.rtIndex = 6;
            index.pepMassIndex = 7;
            index.pepProbcIndex = 8;
            index.purityIndex = 9;
            index.isUniquecIndex = 10;
            index.proteincIndex = 11;
            index.proteinIDcIndex = 12;
            index.genecIndex = 13;
            index.observedMzIndex = 14;

            // Light sub1: abnIndex=15 (2 channels), Light sub2: abnIndex=17 (1 channel)
            Index.SubplexIndex lightSi1 = new Index.SubplexIndex();
            lightSi1.abnIndex = 15;
            lightSi1.refIndex = 15;
            index.subplexIndices.add(lightSi1);

            Index.SubplexIndex lightSi2 = new Index.SubplexIndex();
            lightSi2.abnIndex = 17;
            lightSi2.refIndex = 17;
            index.subplexIndices.add(lightSi2);

            // Heavy sub1: abnIndex=18 (2 channels), Heavy sub2: abnIndex=20 (1 channel)
            Index.SubplexIndex heavySi1 = new Index.SubplexIndex();
            heavySi1.abnIndex = 18;
            heavySi1.refIndex = 18;
            index.subplexIndices.add(heavySi1);

            Index.SubplexIndex heavySi2 = new Index.SubplexIndex();
            heavySi2.abnIndex = 20;
            heavySi2.refIndex = 20;
            index.subplexIndices.add(heavySi2);

            String[] fields = {
                "n[229.16]PEPTIDE",      // 0
                "run1.00100.00100.3",     // 1: Light Spectrum
                "run1.00100.00100.3",     // 2: Heavy Spectrum
                "1K(229.1629)",           // 3: Light Assigned Mods
                "20000",                  // 4: Light Intensity
                "30000",                  // 5: Heavy Intensity
                "25.5",                   // 6: RT
                "850.4",                  // 7: Calc Pep Mass
                "0.98",                   // 8: Probability
                "0.92",                   // 9: Purity
                "true",                   // 10: Is Unique
                "sp|P12345|PROT_HUMAN",   // 11: Protein
                "P12345",                 // 12: Protein ID
                "GENE1",                  // 13: Gene
                "425.2",                  // 14: Observed M/Z
                "1000", "2000",           // 15-16: Light sub1 channels (sum=3000)
                "3000",                   // 17: Light sub2 channel (sum=3000)
                "4000", "5000",           // 18-19: Heavy sub1 channels (sum=9000)
                "6000",                   // 20: Heavy sub2 channel (sum=6000)
            };

            PsmRecord record = new PsmRecord(params, index);
            record.parsePsmRecord(String.join("\t", fields), new ArrayList<>(), new ArrayList<>());

            double[] subplexMs1 = getField(record, "subplexMs1Intensities");
            assertEquals(4, subplexMs1.length);
            // Light total MS1 = 20000, Light TMT sums = [3000, 3000], total = 6000
            assertEquals(20000.0 * 3000.0 / 6000.0, subplexMs1[0], 1e-6); // Light sub1
            assertEquals(20000.0 * 3000.0 / 6000.0, subplexMs1[1], 1e-6); // Light sub2
            // Heavy total MS1 = 30000, Heavy TMT sums = [9000, 6000], total = 15000
            assertEquals(30000.0 * 9000.0 / 15000.0, subplexMs1[2], 1e-6); // Heavy sub1
            assertEquals(30000.0 * 6000.0 / 15000.0, subplexMs1[3], 1e-6); // Heavy sub2
            assertEquals(50000.0, record.getMs1Intensity(), 1e-6); // total
        }
    }

    @Nested
    class UseMS1IntensityTest {

        @Test
        void scaledRefIntensityForStandardSingleSubplex() throws Exception {
            PsmRecord record = parseStandardLine(buildStandardFields());

            // Before: refIntensity = sum of channels = 21000 (summation mode)
            // subplexMs1[0] = 50000, subplexTmtSums[0] = 21000
            record.useMS1Intensity();

            // After: refIntensity = ms1 * (oldRef / tmtSum) = 50000 * (21000/21000) = 50000
            List<Double> refs = getField(record, "subplexRefIntensities");
            assertEquals(50000.0, refs.get(0), 1e-6);
        }

        @Test
        void scaledRefIntensityForHyperplexingSubplexes() throws Exception {
            PsmRecord record = parseHyperplexingLine(buildHyperplexingFields());

            // Before: refs = [6000, 15000] (summation of each subplex's channels)
            // subplexMs1 = [20000, 30000], subplexTmtSums = [6000, 15000]
            record.useMS1Intensity();

            // Light: 20000 * (6000/6000) = 20000
            // Heavy: 30000 * (15000/15000) = 30000
            List<Double> refs = getField(record, "subplexRefIntensities");
            assertEquals(20000.0, refs.get(0), 1e-6);
            assertEquals(30000.0, refs.get(1), 1e-6);
        }

        @Test
        void skipsSubplexWithZeroTmtSum() throws Exception {
            String[] fields = buildHyperplexingFields();
            // Set all Light channels to 0
            fields[15] = "0";
            fields[16] = "0";
            fields[17] = "0";

            PsmRecord record = parseHyperplexingLine(fields);
            record.useMS1Intensity();

            // Light tmtSum=0, so ref should remain unchanged (0.0 from summation of zeros)
            List<Double> refs = getField(record, "subplexRefIntensities");
            assertEquals(0.0, refs.get(0), 1e-6);
            // Heavy should still be scaled
            assertEquals(30000.0, refs.get(1), 1e-6);
        }
    }

    @Nested
    class SetActivePlexTest {

        @Test
        void switchesBetweenSubplexes() {
            PsmRecord record = parseHyperplexingLine(buildHyperplexingFields());

            record.setActivePlex(0);
            assertEquals(Arrays.asList(1000.0, 2000.0, 3000.0), record.getChannels());
            assertEquals(6000.0, record.getRefIntensity(), 1e-6);

            record.setActivePlex(1);
            assertEquals(Arrays.asList(4000.0, 5000.0, 6000.0), record.getChannels());
            assertEquals(15000.0, record.getRefIntensity(), 1e-6);
        }
    }

    @Nested
    class PassResolutionSnrTest {

        @Test
        void failsWhenSnrSumBelowThreshold() throws Exception {
            Parameters params = buildStandardParams();
            params.addIsobaricFilter = true;
            params.minSNR = 100.0f;
            params.minResolution = 45000;

            Index index = buildStandardIndex();
            index.subplexIndices.get(0).resOffset = 23;
            index.subplexIndices.get(0).snrOffset = 29;

            PsmRecord record = new PsmRecord(params, index);

            // Append 6 Resolution + 6 SNR columns (all low SNR)
            String[] base = buildStandardFields();
            String[] fields = new String[base.length + 12];
            System.arraycopy(base, 0, fields, 0, base.length);
            // Resolution columns (col 23-28)
            for (int i = 0; i < 6; i++) fields[23 + i] = "50000";
            // SNR columns (col 29-34) — all low, sum = 6 < 100
            for (int i = 0; i < 6; i++) fields[29 + i] = "1";

            record.parsePsmRecord(String.join("\t", fields), new ArrayList<>(), new ArrayList<>());

            assertFalse((boolean) getField(record, "passResolutionSnr"));
        }

        @Test
        void failsWhenResolutionBelowThreshold() throws Exception {
            Parameters params = buildStandardParams();
            params.addIsobaricFilter = true;
            params.minSNR = 1.0f;
            params.minResolution = 45000;

            Index index = buildStandardIndex();
            index.subplexIndices.get(0).resOffset = 23;
            index.subplexIndices.get(0).snrOffset = 29;

            PsmRecord record = new PsmRecord(params, index);

            String[] base = buildStandardFields();
            String[] fields = new String[base.length + 12];
            System.arraycopy(base, 0, fields, 0, base.length);
            // Resolution columns — one below threshold
            for (int i = 0; i < 6; i++) fields[23 + i] = "50000";
            fields[23] = "10000"; // below minResolution
            // SNR columns — high enough to pass sum but >=1 triggers resolution check
            for (int i = 0; i < 6; i++) fields[29 + i] = "5";

            record.parsePsmRecord(String.join("\t", fields), new ArrayList<>(), new ArrayList<>());

            assertFalse((boolean) getField(record, "passResolutionSnr"));
        }

        @Test
        void passesWhenAllAboveThresholds() throws Exception {
            Parameters params = buildStandardParams();
            params.addIsobaricFilter = true;
            params.minSNR = 1.0f;
            params.minResolution = 45000;

            Index index = buildStandardIndex();
            index.subplexIndices.get(0).resOffset = 23;
            index.subplexIndices.get(0).snrOffset = 29;

            PsmRecord record = new PsmRecord(params, index);

            String[] base = buildStandardFields();
            String[] fields = new String[base.length + 12];
            System.arraycopy(base, 0, fields, 0, base.length);
            for (int i = 0; i < 6; i++) fields[23 + i] = "50000";
            for (int i = 0; i < 6; i++) fields[29 + i] = "5";

            record.parsePsmRecord(String.join("\t", fields), new ArrayList<>(), new ArrayList<>());

            assertTrue((boolean) getField(record, "passResolutionSnr"));
        }
    }

    @Nested
    class IsPassCriteriaTest {

        /**
         * Build a hyperplexing record (Light + Heavy plexes, summation reference) that passes every
         * non-reference criterion, with the given channel intensities so each plex's summation
         * reference can be controlled (a plex with all-zero channels has a zero reference).
         */
        private PsmRecord buildHyperplexRecord(String light1, String light2, String light3,
                                               String heavy1, String heavy2, String heavy3) {
            Parameters params = buildHyperplexingParams();
            params.modTagSet.add("none"); // make modificationFlag pass
            Index index = buildHyperplexingIndex();
            PsmRecord record = new PsmRecord(params, index);
            String[] fields = buildHyperplexingFields();
            fields[15] = light1;
            fields[16] = light2;
            fields[17] = light3;
            fields[18] = heavy1;
            fields[19] = heavy2;
            fields[20] = heavy3;
            record.parsePsmRecord(String.join("\t", fields), new ArrayList<>(), new ArrayList<>());
            record.updateFlags(new HashSet<>());
            return record;
        }

        private PsmRecord buildStandardRecord(String... channelValues) {
            Parameters params = buildStandardParams();
            params.modTagSet.add("none");
            Index index = buildStandardIndex();
            PsmRecord record = new PsmRecord(params, index);
            String[] fields = buildStandardFields();
            for (int i = 0; i < channelValues.length; i++) {
                fields[17 + i] = channelValues[i];
            }
            record.parsePsmRecord(String.join("\t", fields), new ArrayList<>(), new ArrayList<>());
            record.updateFlags(new HashSet<>());
            return record;
        }

        @Test
        void keepsHyperplexPsmWhenOnlyOnePlexReferenceMissing() {
            // Light plex reference (summation) = 0, Heavy reference > 0. A peptide quantified only in
            // the Heavy label must be retained — it is recoverable in the Heavy plex.
            PsmRecord record = buildHyperplexRecord("0", "0", "0", "4000", "5000", "6000");
            assertTrue(record.isPassCriteria(0));
        }

        @Test
        void dropsHyperplexPsmWhenAllPlexReferencesMissing() {
            // No plex has a usable reference -> nothing to quantify anywhere -> drop the PSM.
            PsmRecord record = buildHyperplexRecord("0", "0", "0", "0", "0", "0");
            assertFalse(record.isPassCriteria(0));
        }

        @Test
        void keepsHyperplexPsmWhenAllPlexReferencesPresent() {
            PsmRecord record = buildHyperplexRecord("1000", "2000", "3000", "4000", "5000", "6000");
            assertTrue(record.isPassCriteria(0));
        }

        @Test
        void dropsStandardPsmWhenReferenceMissing() {
            // Standard single-plex: a zero reference must still drop the PSM (unchanged behavior).
            PsmRecord record = buildStandardRecord("0", "0", "0", "0", "0", "0");
            assertFalse(record.isPassCriteria(0));
        }

        @Test
        void keepsStandardPsmWhenReferencePresent() {
            PsmRecord record = buildStandardRecord("1000", "2000", "3000", "4000", "5000", "6000");
            assertTrue(record.isPassCriteria(0));
        }

        /**
         * Build a non-hyperplexing record with two sub-plexes (e.g. TMT-35: an 18-channel and a
         * 17-channel set, here simplified to 3 + 3). Each sub-plex's summation reference is the sum
         * of its own channels, so a sub-plex with all-zero channels has a zero reference.
         */
        private PsmRecord buildTwoSubplexStandardRecord(String s1c1, String s1c2, String s1c3,
                                                        String s2c1, String s2c2, String s2c3) {
            Parameters params = new Parameters();
            params.isHyperplexing = false;
            params.addIsobaricFilter = false;
            params.add_Ref = 0; // summation virtual reference
            params.glycoQval = -1;
            params.allow_unlabeled = true;
            params.modTagSet.add("none");
            params.subplexes.add(new Subplex("ref_126", 3));
            params.subplexes.add(new Subplex("ref_131", 3));

            Index index = buildStandardIndex();
            Index.SubplexIndex si2 = new Index.SubplexIndex();
            si2.abnIndex = 20;
            si2.refIndex = 20;
            index.subplexIndices.add(si2);

            PsmRecord record = new PsmRecord(params, index);
            String[] fields = buildStandardFields();
            fields[17] = s1c1;
            fields[18] = s1c2;
            fields[19] = s1c3;
            fields[20] = s2c1;
            fields[21] = s2c2;
            fields[22] = s2c3;
            record.parsePsmRecord(String.join("\t", fields), new ArrayList<>(), new ArrayList<>());
            record.updateFlags(new HashSet<>());
            return record;
        }

        @Test
        void dropsStandardMultiSubplexPsmWhenOneSubplexReferenceMissing() {
            // Non-hyperplexing multi-sub-plex (e.g. TMT-35): the sub-plexes share one pooled
            // reference, so a peptide must be referenced in EVERY sub-plex. A zero reference in one
            // sub-plex drops the PSM. This preserves the pre-existing TMT-35 behavior; the
            // any-reference relaxation applies to hyperplexing (SILAC) only.
            PsmRecord record = buildTwoSubplexStandardRecord("0", "0", "0", "4000", "5000", "6000");
            assertFalse(record.isPassCriteria(0));
        }

        @Test
        void keepsStandardMultiSubplexPsmWhenAllSubplexReferencesPresent() {
            PsmRecord record = buildTwoSubplexStandardRecord("1000", "2000", "3000", "4000", "5000", "6000");
            assertTrue(record.isPassCriteria(0));
        }
    }

    @Nested
    class HyperplexMultiSubplexReferenceTest {

        /**
         * Build a hyperplexing record with two labels (Light, Heavy), each itself a TMT-35-style pair
         * of sub-plexes (a 2-channel and a 1-channel sub-plex). Each boolean selects whether that
         * sub-plex carries signal; with summation references a sub-plex with all-zero channels has a
         * zero reference. Sub-plex order is label-major: [Light-1, Light-2, Heavy-1, Heavy-2].
         */
        private PsmRecord buildHyperplexTmt35Record(boolean light1, boolean light2,
                                                    boolean heavy1, boolean heavy2) {
            Parameters params = new Parameters();
            params.isHyperplexing = true;
            params.addIsobaricFilter = false;
            params.add_Ref = 0; // summation virtual reference
            params.glycoQval = -1;
            params.allow_unlabeled = true;
            params.modTagSet.add("none");
            params.subplexes.add(new Subplex("ref_126", 2));
            params.subplexes.add(new Subplex("ref_128", 1));
            params.subplexes.add(new Subplex("ref_126", 2));
            params.subplexes.add(new Subplex("ref_128", 1));
            params.hyperplexLabelNames.add("Light");
            params.hyperplexLabelNames.add("Heavy");
            params.hyperplexLabelSubplexCounts.add(2);
            params.hyperplexLabelSubplexCounts.add(2);

            Index index = new Index();
            index.labelFreeModifiedPeptideIndex = 0;
            index.labelSpectrumIndices[0] = 1;
            index.labelSpectrumIndices[2] = 2;
            index.labelAssignedModIndices[0] = 3;
            index.labelMs1IntIndices[0] = 4;
            index.labelMs1IntIndices[2] = 5;
            index.pepcIndex = 0;
            index.spectrumIndex = 1;
            index.rtIndex = 6;
            index.pepMassIndex = 7;
            index.pepProbcIndex = 8;
            index.purityIndex = 9;
            index.isUniquecIndex = 10;
            index.proteincIndex = 11;
            index.proteinIDcIndex = 12;
            index.genecIndex = 13;
            index.observedMzIndex = 14;

            Index.SubplexIndex l1 = new Index.SubplexIndex();
            l1.abnIndex = 15;
            l1.refIndex = 15;
            index.subplexIndices.add(l1);
            Index.SubplexIndex l2 = new Index.SubplexIndex();
            l2.abnIndex = 17;
            l2.refIndex = 17;
            index.subplexIndices.add(l2);
            Index.SubplexIndex h1 = new Index.SubplexIndex();
            h1.abnIndex = 18;
            h1.refIndex = 18;
            index.subplexIndices.add(h1);
            Index.SubplexIndex h2 = new Index.SubplexIndex();
            h2.abnIndex = 20;
            h2.refIndex = 20;
            index.subplexIndices.add(h2);

            String[] fields = {
                "n[229.16]PEPTIDE", "run1.00100.00100.3", "run1.00100.00100.3", "1K(229.1629)",
                "20000", "30000", "25.5", "850.4", "0.98", "0.92", "true",
                "sp|P12345|PROT_HUMAN", "P12345", "GENE1", "425.2",
                light1 ? "1000" : "0", light1 ? "2000" : "0", // 15,16 Light-1
                light2 ? "3000" : "0",                        // 17    Light-2
                heavy1 ? "4000" : "0", heavy1 ? "5000" : "0", // 18,19 Heavy-1
                heavy2 ? "6000" : "0",                        // 20    Heavy-2
            };

            PsmRecord record = new PsmRecord(params, index);
            record.parsePsmRecord(String.join("\t", fields), new ArrayList<>(), new ArrayList<>());
            record.updateFlags(new HashSet<>());
            return record;
        }

        @Test
        void keepsWhenAtLeastOneLabelFullyReferenced() {
            // Light fully referenced (both sub-plexes), Heavy only partially -> keep (Light is usable).
            PsmRecord record = buildHyperplexTmt35Record(true, true, true, false);
            assertTrue(record.isPassCriteria(0));
        }

        @Test
        void dropsWhenNoLabelFullyReferenced() {
            // Each label is missing one of its two sub-plex references, so no label is quantifiable
            // as a TMT-35 unit -> drop, even though individual sub-plexes carry a reference.
            PsmRecord record = buildHyperplexTmt35Record(true, false, true, false);
            assertFalse(record.isPassCriteria(0));
        }

        @Test
        void masksReferencesOfPartiallyReferencedLabel() throws Exception {
            // Light fully referenced, Heavy partial (Heavy-2 missing). The whole Heavy label must be
            // treated as unquantifiable: every Heavy reference is zeroed, including Heavy-1 which had
            // signal. Light references are preserved. Order: [Light-1, Light-2, Heavy-1, Heavy-2].
            PsmRecord record = buildHyperplexTmt35Record(true, true, true, false);
            List<Double> refs = getField(record, "subplexRefIntensities");
            assertTrue(refs.get(0) > 0, "Light-1 reference should be preserved");
            assertTrue(refs.get(1) > 0, "Light-2 reference should be preserved");
            assertEquals(0.0, refs.get(2), 1e-9, "Heavy-1 reference should be masked (Heavy incomplete)");
            assertEquals(0.0, refs.get(3), 1e-9, "Heavy-2 reference is absent");
        }
    }
}
