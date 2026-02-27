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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import tmtintegrator.constants.GroupBy;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;

/**
 * PSM class represents a PSM.tsv file
 */
public class Psm {

    private final File psmFile;
    private final Parameters parameters;
    private final Index index;
    private final List<PsmRecord> psmRecords;
    private boolean isFirstProcess;
    private final List<Integer> naChannels;

    public List<PsmRecord> getPsmRecords() {
        return psmRecords;
    }

    public File getPsmFile() {
        return psmFile;
    }

    public Index getIndex() {
        return index;
    }

    public Psm(Parameters parameters, File psmFile) {
        this.psmFile = psmFile;
        this.parameters = parameters;
        this.index = parameters.indMap.get(psmFile.getAbsolutePath());
        this.psmRecords = new ArrayList<>();
        this.isFirstProcess = true;
        this.naChannels = new ArrayList<>();
    }

    /**
     * Read PSM file and record title and PSM records
     */
    public void readPsmFile(File psmFile, List<Double> tmtIntensities) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(psmFile))) {
            String title = reader.readLine();
            title = removeNaChannels(title);
            parameters.titleMap.put(psmFile.getAbsolutePath(), title);

            String line;
            while ((line = reader.readLine()) != null) {
                PsmRecord psmRecord = new PsmRecord(parameters, index);
                psmRecord.parsePsmRecord(line, tmtIntensities, naChannels);
                psmRecords.add(psmRecord);
            }
        }
    }

    /**
     * Adjust plex channel offsets for NA channels in preceding plexes.
     * Each plex's NA channels shift all subsequent plex offsets left.
     */
    public void adjustPlexChannelOffsets() {
        int cumulativeGap = 0;
        for (int i = 0; i < index.subplexIndices.size(); i++) {
            Index.SubplexIndex si = index.subplexIndices.get(i);
            si.abnIndex -= cumulativeGap;
            si.refIndex -= cumulativeGap;
            int expectedChannels = parameters.subplexes.get(i).channelCount;
            cumulativeGap += expectedChannels - si.usedChannelNum;
        }
    }

    /**
     * Filter out PSMs that are not used in psmMap
     */
    public void filterUnUsedPsm(Map<String, List<PsmRecord>> psmMap) {
        psmRecords.clear();
        for (List<PsmRecord> psmList : psmMap.values()) {
            for (PsmRecord psmRecord : psmList) {
                if (psmRecord.isUsed()) {
                    psmRecords.add(psmRecord);
                }
            }
        }
    }

    /**
     * Backup PSM records for each groupBy
     */
    public void backup() {
        if (isFirstProcess) {
            for (PsmRecord psmRecord : psmRecords) {
                psmRecord.backup();
            }
            isFirstProcess = false;
        }
    }

    /**
     * Analyze phospho sites
     */
    public void analyzeByGroup(GroupBy groupBy) {
        for (PsmRecord psmRecord : psmRecords) {
            psmRecord.analyzeByGroup(groupBy);
        }
        psmRecords.removeIf(PsmRecord::isExcluded);
    }

    /**
     * Use MS1 intensity as reference
     */
    public void useMS1Intensity() {
        for (PsmRecord psmRecord : psmRecords) {
            psmRecord.useMS1Intensity();
        }
    }

    /**
     * Set the active plex for processing. Sets index fields and PsmRecord active channels.
     */
    public void setActivePlex(int plexIdx) {
        index.setActivePlex(plexIdx);
        for (PsmRecord psmRecord : psmRecords) {
            psmRecord.setActivePlex(plexIdx);
        }
    }

    /**
     * Reset modified fields for next groupBy iteration
     */
    public void reset() {
        for (PsmRecord psmRecord : psmRecords) {
            psmRecord.reset();
        }
    }

    // region helper methods
    private String removeNaChannels(String title) {
        String[] titleArr = title.split("\t");
        StringBuilder newTitle = new StringBuilder();
        for (int i = 0; i < titleArr.length; i++) {
            if (!titleArr[i].trim().equalsIgnoreCase("Intensity NA")) {
                newTitle.append(titleArr[i]).append("\t");
            } else {
                naChannels.add(i);
            }
        }
        return newTitle.toString().trim();
    }
    // endregion
}
