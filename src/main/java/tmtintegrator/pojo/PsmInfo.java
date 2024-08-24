package tmtintegrator.pojo;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

public class PsmInfo {

    public String gene = "";
    public double totalRefInt = 0;
    public List<String> psmList = new ArrayList<>();
    private Set<P> pp = new TreeSet<>();

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
