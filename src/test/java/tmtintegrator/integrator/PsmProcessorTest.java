package tmtintegrator.integrator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static tmtintegrator.integrator.PsmProcessor.ff;

import org.junit.jupiter.api.Test;

class PsmProcessorTest {

  @Test
  void testFt() {
    String[] ss = ff("ABCDEFGH", "IJKLMNOPQRST", "UVWXYZ", 104, 111, 100, 113);
    assertEquals("EFGHIJK", ss[0]);
    assertEquals("LMNOPQRS", ss[1]);
    assertEquals("TUVWXYZ", ss[2]);

    ss = ff("GH", "IJKLMNOPQRST", "U", 104, 111, 100, 113);
    assertEquals("__GHIJK", ss[0]);
    assertEquals("LMNOPQRS", ss[1]);
    assertEquals("TU_____", ss[2]);

    ss = ff("ABCDEFGHIJK", "LMNOPQRST", "UVWXYZ", 104, 111, 103, 113);
    assertEquals("EFGHIJK", ss[0]);
    assertEquals("LMNOPQRS", ss[1]);
    assertEquals("TUVWXYZ", ss[2]);

    ss = ff("ABCDEFGHIJ", "KLMNOPQRST", "UVWXYZ", 104, 111, 102, 113);
    assertEquals("EFGHIJK", ss[0]);
    assertEquals("LMNOPQRS", ss[1]);
    assertEquals("TUVWXYZ", ss[2]);

    ss = ff("ABCDE", "FGHIJKLMNOPQRST", "UVWXYZ", 104, 111, 97, 113);
    assertEquals("EFGHIJK", ss[0]);
    assertEquals("LMNOPQRS", ss[1]);
    assertEquals("TUVWXYZ", ss[2]);

    ss = ff("ABCD", "EFGHIJKLMNOPQRST", "UVWXYZ", 104, 111, 96, 113);
    assertEquals("EFGHIJK", ss[0]);
    assertEquals("LMNOPQRS", ss[1]);
    assertEquals("TUVWXYZ", ss[2]);

    ss = ff("ABCDEFGH", "IJKLMNOPQRSTUVWXY", "Z", 104, 111, 100, 118);
    assertEquals("EFGHIJK", ss[0]);
    assertEquals("LMNOPQRS", ss[1]);
    assertEquals("TUVWXYZ", ss[2]);

    ss = ff("ABCDEFGH", "IJKLMNOPQRSTUVWXYZ", "", 104, 111, 100, 119);
    assertEquals("EFGHIJK", ss[0]);
    assertEquals("LMNOPQRS", ss[1]);
    assertEquals("TUVWXYZ", ss[2]);
  }
}