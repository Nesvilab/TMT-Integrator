package tmtintegrator.pojo.psm;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Nested;

import static org.junit.jupiter.api.Assertions.*;

class PsmTest {

    @Nested
    class IsNaChannelTest {

        @Test
        void detectsStandardNa() {
            assertTrue(Psm.isNaChannel("Intensity NA"));
        }

        @Test
        void detectsStandardNaCaseInsensitive() {
            assertTrue(Psm.isNaChannel("intensity na"));
            assertTrue(Psm.isNaChannel("INTENSITY NA"));
        }

        @Test
        void detectsLightNa() {
            assertTrue(Psm.isNaChannel("Light Intensity NA"));
        }

        @Test
        void detectsMediumNa() {
            assertTrue(Psm.isNaChannel("Medium Intensity NA"));
        }

        @Test
        void detectsHeavyNa() {
            assertTrue(Psm.isNaChannel("Heavy Intensity NA"));
        }

        @Test
        void detectsNaWithWhitespace() {
            assertTrue(Psm.isNaChannel("  Intensity NA  "));
            assertTrue(Psm.isNaChannel("  Light Intensity NA  "));
        }

        @Test
        void rejectsRegularIntensityColumn() {
            assertFalse(Psm.isNaChannel("Intensity 126"));
        }

        @Test
        void rejectsLightIntensityColumn() {
            assertFalse(Psm.isNaChannel("Light Intensity run_126"));
        }

        @Test
        void rejectsEmptyString() {
            assertFalse(Psm.isNaChannel(""));
        }

        @Test
        void rejectsPartialMatch() {
            assertFalse(Psm.isNaChannel("Intensity"));
            assertFalse(Psm.isNaChannel("NA"));
        }
    }
}
