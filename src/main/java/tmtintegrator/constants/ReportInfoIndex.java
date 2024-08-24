package tmtintegrator.constants;

/**
 * Column index for report info columns.
 */
public class ReportInfoIndex {

    private ReportInfoIndex() {
        throw new AssertionError("ReportInfoIndex is a utility class");
    }

    /**
     * Column format for report info part for peptide and site level:
     * <Index> <Gene> <ProteinID> <Peptide> <SequenceWindow> <Start> <End> <MaxPepProb>
     */
    public static final int INDEX = 0;
    public static final int GENE = 1;
    public static final int PROTEIN_ID = 2;
    public static final int PEPTIDE = 3;
    public static final int SEQUENCE_WINDOW = 4;
    public static final int START = 5;
    public static final int END = 6;
    public static final int MAX_PEP_PROB = 7;

    /**
     * Column format for report info part for protein level:
     * <Index> <NumberPSM> <Gene> <MaxPepProb>
     */
    public static final int NUMBER_PSM = 1;
    public static final int GENE_PROTEIN = 2;
    public static final int MAX_PEP_PROB_PROTEIN = 3;
}
