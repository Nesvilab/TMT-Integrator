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

package tmtintegrator.pojo.psm;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class PhosphoSiteData {
    int siteLocalCount;
    final StringBuilder siteLocalPos; // site localized position
    final Map<Integer, List<String>> siteLocalMassMap; // Key: position; Value: site localized mass
    int siteCount;
    final Map<Integer, Double> probMap; // Key: position; Value: site probability
    int startIndex;
    int endIndex;

    PhosphoSiteData() {
        siteLocalCount = 0;
        siteLocalPos = new StringBuilder();
        siteLocalMassMap = new TreeMap<>();
        siteCount = 0;
        probMap = new TreeMap<>();
    }
}
