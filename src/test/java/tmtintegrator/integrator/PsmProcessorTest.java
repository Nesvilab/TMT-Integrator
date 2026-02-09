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