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

import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.Subplex;
import tmtintegrator.pojo.psm.PsmRecord;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class UtilsTest {

    @Nested
    class LogNormalizePsmTest {

        // region reflection helpers

        private void invokeLogNormalizePsm(PsmRecord record) throws Exception {
            Method m = Utils.class.getDeclaredMethod("logNormalizePsm", PsmRecord.class);
            m.setAccessible(true);
            m.invoke(null, record);
        }

        /**
         * Build a single-plex PsmRecord with the given active-plex reference intensity and channels,
         * bypassing file parsing so the reference can be set independently of the channel values
         * (as happens with a real reference channel that is empty while sample channels carry signal).
         */
        private PsmRecord buildRecord(double refIntensity, List<Double> channels) throws Exception {
            Parameters params = new Parameters();
            params.subplexes.add(new Subplex("ref", channels.size()));
            Index index = new Index();
            PsmRecord record = new PsmRecord(params, index);

            Field subplexChannelsField = PsmRecord.class.getDeclaredField("subplexChannels");
            subplexChannelsField.setAccessible(true);
            List<List<Double>> subplexChannels = new ArrayList<>();
            subplexChannels.add(new ArrayList<>(channels));
            subplexChannelsField.set(record, subplexChannels);

            Field subplexRefField = PsmRecord.class.getDeclaredField("subplexRefIntensities");
            subplexRefField.setAccessible(true);
            List<Double> subplexRefs = new ArrayList<>();
            subplexRefs.add(refIntensity);
            subplexRefField.set(record, subplexRefs);

            record.setActivePlex(0);
            return record;
        }

        // endregion

        @Test
        void missingReferenceProducesAllNaNChannels() throws Exception {
            // Real reference channel empty (0) for this plex, but sample channels carry signal.
            // Ratios are undefined without a reference, so every channel must become NaN
            // rather than an un-normalized absolute log value.
            PsmRecord record = buildRecord(0.0, List.of(1000.0, 2000.0, 3000.0));

            invokeLogNormalizePsm(record);

            for (double channel : record.getChannels()) {
                assertTrue(Double.isNaN(channel), "channel should be NaN when the plex reference is missing");
            }
        }

        @Test
        void positiveReferenceProducesLogRatios() throws Exception {
            PsmRecord record = buildRecord(6000.0, List.of(1000.0, 2000.0, 3000.0));

            invokeLogNormalizePsm(record);

            List<Double> channels = record.getChannels();
            assertEquals(Utils.log2(1000.0) - Utils.log2(6000.0), channels.get(0), 1e-9);
            assertEquals(Utils.log2(2000.0) - Utils.log2(6000.0), channels.get(1), 1e-9);
            assertEquals(Utils.log2(3000.0) - Utils.log2(6000.0), channels.get(2), 1e-9);
        }

        @Test
        void zeroIntensityChannelIsNaNWithValidReference() throws Exception {
            PsmRecord record = buildRecord(6000.0, List.of(1000.0, 0.0, 3000.0));

            invokeLogNormalizePsm(record);

            List<Double> channels = record.getChannels();
            assertEquals(Utils.log2(1000.0) - Utils.log2(6000.0), channels.get(0), 1e-9);
            assertTrue(Double.isNaN(channels.get(1)));
            assertEquals(Utils.log2(3000.0) - Utils.log2(6000.0), channels.get(2), 1e-9);
        }
    }
}
