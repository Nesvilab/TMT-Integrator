package tmtintegrator.config;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Parses command line arguments.
 *
 * @author rogerli on 05/2024
 */
public class ArgumentParser {

    private final File yamlFile;
    private final List<File> inputFiles = new ArrayList<>();
    private final boolean validateParameters;

    public ArgumentParser(String[] args) {
        // validate arguments
        if (args.length < 2) {
            throw new IllegalArgumentException("Missing arguments");
        }
        if (!args[0].endsWith(".yml")) {
            throw new IllegalArgumentException("First argument must be a YAML file");
        }

        // parse arguments
        yamlFile = new File(args[0]);
        if (args.length == 2 && args[1].equalsIgnoreCase("--ValParam")) {
            this.validateParameters = true;
        } else {
            this.validateParameters = false;
            for (int i = 1; i < args.length; i++) {
                inputFiles.add((new File(args[i])).getAbsoluteFile());
            }
        }
    }

    public File getYamlFile() {
        return yamlFile;
    }

    public List<File> getInputFiles() {
        return inputFiles;
    }

    public boolean isValidateOnly() {
        return validateParameters;
    }
}
