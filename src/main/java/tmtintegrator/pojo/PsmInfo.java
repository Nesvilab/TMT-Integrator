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

import tmtintegrator.pojo.psm.PsmRecord;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

public class PsmInfo {

    public String gene = "";
    public double totalRefInt = 0;
    public List<PsmRecord> psmRecords = new ArrayList<>();
    private final Set<P> pp = new TreeSet<>();

    public void addP(String peptide, String l, String r, int pepsIndex) {
        P p1 = new P(peptide, l, r, pepsIndex);
        pp.add(p1);
    }

    public Set<String> getPeptideIndex() {
        Set<String> s = new HashSet<>();
        for (P p : pp) {
            s.add(p.toString());
        }
        return s;
    }


    public static class P implements Comparable<P> {

        public final String peptide, l, r;
        public final int pepsIndex;

        public P(String peptide, String l, String r, int pepsIndex) {
            this.peptide = peptide;
            this.l = l;
            this.r = r;
            this.pepsIndex = pepsIndex;
        }

        @Override
        public int compareTo(P o) {
            return peptide.compareTo(o.peptide);
        }

        @Override
        public boolean equals(Object o) {
            if (o instanceof P) {
                return peptide.compareTo(((P) o).peptide) == 0;
            }
            return false;
        }

        @Override
        public int hashCode() {
            return peptide.hashCode();
        }

        @Override
        public String toString() {
            return peptide + "@" + pepsIndex + "@" + l + "." + peptide + "." + r;
        }
    }
}
