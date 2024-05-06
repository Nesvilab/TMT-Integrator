package tmtintegrator.config;

import tmtintegrator.pojo.ds_Index;
import tmtintegrator.pojo.ds_Parameters;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

public class ConfigLoader {
    private final ds_Parameters parameters;

    public ConfigLoader() {
        this.parameters = new ds_Parameters();
    }

    public ds_Parameters getParameters() {
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
                parseLine(line);
            }
        }
    }

    /**
     * Load a list of input files
     *
     * @param inputFiles List of input files from the command line
     */
    public void loadFileList(List<File> inputFiles) {
        parameters.FileLi.addAll(inputFiles);
        // add file names to fNameLi and indMap
        for (File file : inputFiles) {
            String absolutePath = file.getAbsolutePath();
            parameters.fNameLi.add(absolutePath);

            ds_Index index = new ds_Index();
            parameters.indMap.put(absolutePath, index);
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
            // TODO: handle invalid line
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
                case "prefix":
                    parameters.prefix = value;
                    break;
                case "log2transformed":
                    parameters.log2transformed = Boolean.parseBoolean(value);
                    break;
            }
        } catch (NumberFormatException e) {
            System.err.println("Error parsing value: " + value);
            throw new RuntimeException("Failed to parse YAML file", e);
        } catch (Exception e) {
            System.err.println("Error parsing line: " + line);
            throw new RuntimeException("Failed to parse YAML file", e);
        }
    }

    private void handleModTag(String value) {
        value = (value.isEmpty()) ? "none" : value;
        if (value.equalsIgnoreCase("none")) {
            parameters.modTagLi.add(value);
            return;
        }

        StringBuilder modifiedAA = new StringBuilder();
        String[] parts = value.split(",");
        String targetMass = "";

        for (int i = 0; i < parts.length; i++) {
            String mass = extractMass(parts[i]);
            String formattedMass = formatValue(parts[i], mass);
            parameters.modTagLi.add(formattedMass.trim());

            if (i == 0 && !mass.isEmpty()) {
                targetMass = mass;
            }

            if (formattedMass.equalsIgnoreCase("n-glyco")
                    || formattedMass.equalsIgnoreCase("o-glyco")) {
                parameters.glycoflag = true;
            }

            if (parts[i].contains("(")) {
                int end = parts[i].indexOf("(");
                String aminoAcid = parts[i].substring(end - 1, end);
                if (!modifiedAA.toString().contains(aminoAcid)) {
                    modifiedAA.append(aminoAcid).append("|");
                }
            }
        }

        // Remove the last pipe character if it exists
        if (modifiedAA.length() > 0) {
            modifiedAA.deleteCharAt(modifiedAA.length() - 1);
        }
        parameters.modAA = modifiedAA.toString();
        parameters.columntag = targetMass.isEmpty() ? "" : modAAWithMass(parameters.modAA, targetMass);
    }

    private String extractMass(String value) {
        if (value.contains("(") && value.contains(")")) {
            int start = value.indexOf("(");
            int end = value.indexOf(")");
            return value.substring(start + 1, end);
        }
        return "";
    }

    private String formatValue(String value, String mass) {
        if (value.contains("(")) {
            return value.substring(0, value.indexOf("(")) + "(" + mass + ")";
        }
        return value;
    }

    private String modAAWithMass(String modAA, String mass) {
        // extract amino acids
        modAA = modAA.replace("|", "");
        // truncate the mass to 2 decimal places
        // TODO: check if this is the correct way to truncate the mass
        // mass = String.format("%.2f", Float.parseFloat(mass));
        mass = mass.substring(0, mass.indexOf(".") + 3);
        return modAA + ":" + mass;

    }
    // endregion helper methods====================================================================================
}
