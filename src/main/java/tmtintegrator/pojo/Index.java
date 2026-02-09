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

package tmtintegrator.pojo;

public class Index {
    public int spectrumIndex = -1;
    public int pepcIndex = -1;
    public int peptideIndex = -1;
    public int modifiedPeptideIndex = -1;
    public int extpepIndex = -1;
    public int chargeIndex = -1;
    public int rtIndex = -1;
    public int observedMzIndex = -1;
    public int pepMassIndex = -1;
    public int pepProbcIndex = -1;
    public int numEnzyTermi = -1;
    public int protsIndex = -1;
    public int proteIndex = -1;
    public int ms1IntIndex = -1;
    public int assignedModcIndex = -1;
    public int observedModIndex = -1;
    public int purityIndex = -1;
    public int isUniquecIndex = -1;
    public int proteincIndex = -1;
    public int proteinIDcIndex = -1;
    public int entryNameIndex = -1;
    public int genecIndex = -1;
    public int proteinDescIndex = -1;
    public int mapGeneIndex = -1;
    public int mappedProteinsIndex = -1;
    public int ptmLocalcIndex = -1;
    public int glycoCompositionIndex = -1;
    public int glycoQvalIndex = -1;
    public int usedChannelNum = -1;
    public int usedDChannelNum = -1; // number of used deuterium channels
    public int copyUsedChannelNum = -1; // for 2 rounds TMT-35 resetting
    public int refIndex = -1;
    public int refDIndex = -1; // index of deuterium reference channel
    public int copyRefIndex = -1; // for 2 rounds TMT-35 resetting
    public int abnIndex = -1; // offset for channel intensity
    public int abnDIndex = -1; // offset for deuterium channel intensity
    public int allChannelOffset = -1; // offset for all channel intensity (for TMT-35 resetting)
    public int resOffset = -1; // offset for channel resolution
    public int snrOffset = -1; // offset for channel SNR
}
