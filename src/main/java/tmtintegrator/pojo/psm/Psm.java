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
     *
     * @param psmFile        PSM file
     * @param tmtIntensities TMT intensities
     * @throws IOException if an I/O error occurs
     */
    public void readPsmFile(File psmFile, List<Double> tmtIntensities) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(psmFile))) {
            // remove NA channels from title and record it
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
     * Adjust deuterium channel offset for NA channels in non-deuterium channels
     */
    public void adjustDChannelOffset() {
        if (parameters.isTmt35) {
            int gap = 18 - index.usedChannelNum;
            index.abnDIndex -= gap;
            index.refDIndex -= gap;
        }
    }

    /**
     * Filter out PSMs that are not used in psmMap
     *
     * @param psmMap map of PSMs
     */
    public void filterUnUsedPsm(Map<String, List<PsmRecord>> psmMap) {
        // Filter out PSMs that are not used in psmMap
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
     * Reset PSM records for each groupBy
     */
    public void resetPsmRecords() {
        if (isFirstProcess) {
            for (PsmRecord psmRecord : psmRecords) {
                psmRecord.backup();
            }
            if (parameters.isTmt35) {
                index.copyUsedChannelNum = index.usedChannelNum;
                index.allChannelOffset = index.abnIndex;
                index.copyRefIndex = index.refIndex;
            }
            isFirstProcess = false;
        } else {
            for (PsmRecord psmRecord : psmRecords) {
                psmRecord.reset();
            }
        }
    }

    /**
     * Analyze phospho sites
     *
     * @param groupBy group by option
     */
    public void analyzeByGroup(GroupBy groupBy) {
        for (PsmRecord psmRecord : psmRecords) {
            psmRecord.analyzeByGroup(groupBy);
        }
        // exclude PSMs with isExcluded set
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
     * Set deuterium channels (2nd round of TMT-35)
     */
    public void setDChannels() {
        index.usedChannelNum = index.usedDChannelNum;
        index.abnIndex = index.abnDIndex;
        index.refIndex = index.refDIndex;
        for (PsmRecord psmRecord : psmRecords) {
            psmRecord.setDChannels();
        }
    }

    /**
     * Reset deuterium channels (after 2nd round of TMT-35)
     */
    public void resetDChannels() {
        if (parameters.isTmt35) {
            index.usedChannelNum = index.copyUsedChannelNum;
            index.abnIndex = index.allChannelOffset;
            index.refIndex = index.copyRefIndex;
        }
    }

    // region helper methods
    private String removeNaChannels(String title) {
        String[] titleArr = title.split("\t");
        StringBuilder newTitle = new StringBuilder();
        // remove NA channels from title and record NA index
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
