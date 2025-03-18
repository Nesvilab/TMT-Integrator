package tmtintegrator.config;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.stream.Collectors;

import tmtintegrator.constants.Constants;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;
import tmtintegrator.pojo.ProteinIndex;

/**
 * Load configuration parameters from a YAML file and command line arguments.
 *
 * @author rogerli on 05/2024
 */
public class ConfigLoader {
    private final Parameters parameters;

    public ConfigLoader() {
        this.parameters = new Parameters();
    }

    public Parameters getParameters() {
        return parameters;
    }

    /**
     * Parse a YAML file line by line
     *
     * @param yamlFile YAML file to parse
     * @throws IOException if an I/O error occurs
     */
    public void loadParameters(File yamlFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(yamlFile))) {
            String line;
            while ((line = reader.readLine()) != null) {
                parseLine(line.trim());
            }
        }

        if (parameters.labels == null || parameters.labels.length == 0 || (parameters.labels.length == 1 && parameters.labels[0] == 0)) {
            System.out.println("Warning: The label mass is unknown. The allow_unlabeled parameter will be set to true.");
            parameters.allow_unlabeled = true;
        }

        if (parameters.isTmt35 && parameters.channelNum != 35) {
            throw new IllegalArgumentException("The channel number for TMT 35-plex must be " + 35);
        }
    }

    /**
     * Load a list of input files
     *
     * @param inputFiles List of input files from the command line
     */
    public void loadFileList(List<File> inputFiles) {
        parameters.fileList.addAll(inputFiles);
        // add file names to fNameLi and indMap
        for (File file : inputFiles) {
            String absolutePath = file.getAbsolutePath();
            parameters.fNameLi.add(absolutePath);
            // Check existence of protein.tsv file in the same folder and record if exists
            String proteinPath = absolutePath.replace("psm.tsv", "protein.tsv");
            File proteinFile = new File(proteinPath);
            if (proteinFile.exists()) {
                parameters.proteinFileList.add(proteinFile);
            } else {
                System.out.println("Protein file not found: " + proteinPath);
            }

            Index index = new Index();
            ProteinIndex proteinIndex = new ProteinIndex();
            parameters.indMap.put(absolutePath, index);
            parameters.proteinIndexMap.put(proteinPath, proteinIndex);
        }
    }

    // region helper methods=======================================================================================
    private void parseLine(String line) {

        // Skip comments and empty lines
        if (line.startsWith("#") || line.trim().isEmpty()) {
            return;
        }

        // Split the line into key and value
        String[] parts = line.split(":", 2);
        if (parts.length < 2) {
            System.out.println("Error parsing line: " + line + ". Ignore it.");
            return;
        }

        String key = parts[0].trim();
        String value = parts[1].trim().split("#")[0].trim(); // Remove comments

        try {
            switch (key) {
                case "output":
                    parameters.reportPath = value;
                    break;
                case "combined_protein":
                    parameters.combinedF = new File(value);
                    break;
                case "channel_num":
                    parameters.channelNum = Integer.parseInt(value);
                    break;
                case "ref_tag":
                    parameters.refTag = value;
                    break;
                case "ref_d_tag":
                    parameters.refDTag = value;
                    break;
                case "is_tmt_35":
                    parameters.isTmt35 = Boolean.parseBoolean(value);
                    break;
                case "groupby":
                    parameters.groupBy = Integer.parseInt(value);
                    break;
                case "psm_norm":
                    parameters.psmNorm = Boolean.parseBoolean(value);
                    break;
                case "outlier_removal":
                    parameters.outlierRemoval = Boolean.parseBoolean(value);
                    break;
                case "prot_norm":
                    parameters.protNorm = Integer.parseInt(value);
                    break;
                case "min_pep_prob":
                    parameters.minPepProb = Float.parseFloat(value);
                    break;
                case "min_purity":
                    parameters.minPurity = Float.parseFloat(value);
                    break;
                case "min_percent":
                    parameters.minPercent = Float.parseFloat(value);
                    break;
                case "unique_pep":
                    parameters.uniquePep = Boolean.parseBoolean(value);
                    break;
                case "unique_gene":
                    parameters.uniqueGene = Integer.parseInt(value);
                    break;
                case "best_psm":
                    parameters.bestPsm = Boolean.parseBoolean(value);
                    break;
                case "prot_exclude":
                    parameters.protExcludeAry = value.isEmpty() ? new String[]{"none"} : value.split(",");
                    break;
                case "allow_overlabel":
                    parameters.allow_overlabel = Boolean.parseBoolean(value);
                    break;
                case "allow_unlabeled":
                    parameters.allow_unlabeled = Boolean.parseBoolean(value);
                    break;
                case "mod_tag":
                    handleModTag(value);
                    break;
                case "min_site_prob":
                    parameters.minSiteProb = Float.parseFloat(value);
                    break;
                case "ms1_int":
                    parameters.ms1Int = Boolean.parseBoolean(value);
                    break;
                case "top3_pep":
                    parameters.top3Pep = Boolean.parseBoolean(value);
                    break;
                case "print_RefInt":
                    parameters.print_RefInt = Boolean.parseBoolean(value);
                    break;
                case "add_Ref":
                    parameters.add_Ref = Integer.parseInt(value);
                    break;
                case "max_pep_prob_thres":
                    parameters.max_pep_prob_thres = Double.parseDouble(value);
                    break;
                case "min_ntt":
                    parameters.min_ntt = Integer.parseInt(value);
                    break;
                case "abn_type":
                    parameters.abn_type = Integer.parseInt(value);
                    break;
                case "aggregation_method":
                    parameters.aggregation_method = Integer.parseInt(value);
                    break;
                case "glyco_qval":
                    parameters.glycoQval = Float.parseFloat(value);
                    break;
                case "use_glycan_composition":
                    parameters.useGlycoComposition = Boolean.parseBoolean(value);
                    break;
                case "log2transformed":
                    parameters.log2transformed = Boolean.parseBoolean(value);
                    break;
                case "label_masses":
                    Set<Float> ff = Arrays.stream(value.split(",")).map(Float::parseFloat).collect(Collectors.toSet());
                    parameters.labels = new float[ff.size()];
                    int i = 0;
                    for (Float f : ff) {
                        parameters.labels[i++] = f;
                    }
                    Arrays.sort(parameters.labels);
                    break;
                case "min_resolution":
                    parameters.minResolution = Integer.parseInt(value);
                    break;
                case "min_snr":
                    parameters.minSNR = Float.parseFloat(value);
                    break;
            }
        } catch (Exception e) {
            System.err.println("Error parsing line: " + line);
            throw new RuntimeException("Failed to parse YAML file", e);
        }
    }

    private void handleModTag(String value) {
        if (value.isEmpty() || value.equalsIgnoreCase("none")) {
            parameters.modTagSet.add("none");
            return;
        }

        StringBuilder modifiedAA = new StringBuilder();
        String[] parts = value.split(",");
        String targetMass = "";

        for (int i = 0; i < parts.length; i++) {
            String modTag = parts[i].trim();
            String mass = extractMass(modTag);
            parameters.modTagSet.add(modTag);

            if (i == 0 && !mass.isEmpty()) {
                targetMass = mass;
            }

            if (modTag.equalsIgnoreCase("n-glyco")
                    || modTag.equalsIgnoreCase("o-glyco")) {
                parameters.glycoflag = true;
            }

            extractAA(modTag, modifiedAA);
        }

        // Remove the last pipe character if it exists
        if (modifiedAA.length() > 0) {
            modifiedAA.deleteCharAt(modifiedAA.length() - 1);
        }
        parameters.modAA = modifiedAA.toString();
        parameters.columntag = targetMass.isEmpty() ? "" : modAAWithMass(parameters.modAA, targetMass);
    }

    private void extractAA(String modTag, StringBuilder modifiedAA) {
        Matcher matcher = Constants.MOD_TAG_PATTERN.matcher(modTag); // Match Pattern: A(123.456)
        if (matcher.find()) {
            String aminoAcid = matcher.group(Constants.AA_GROUP);
            if (!modifiedAA.toString().contains(aminoAcid)) {
                modifiedAA.append(aminoAcid).append("|");
            }
        }
    }

    private String extractMass(String modTag) {
        Matcher matcher = Constants.MOD_TAG_PATTERN.matcher(modTag); // Match Pattern: A(123.456)
        if (matcher.find()) {
            return matcher.group(Constants.MASS_GROUP);
        }
        return ""; // for glyco
    }

    private String modAAWithMass(String modAA, String mass) {
        // extract amino acids
        modAA = modAA.replace("|", "");
        // truncate the mass to 2 decimal places
        // FIXME 03: not a good practice to truncate mass number
        // mass = String.format("%.2f", Float.parseFloat(mass));
        mass = mass.substring(0, mass.indexOf(".") + 3);
        return modAA + ":" + mass;

    }
    // endregion helper methods====================================================================================
}
