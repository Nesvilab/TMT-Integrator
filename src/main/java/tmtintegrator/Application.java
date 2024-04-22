package tmtintegrator;

import tmtintegrator.config.ArgumentParser;

import java.io.IOException;

public class Application {

    private final static String APP_NAME = "TMT Integrator";
    private final static String APP_VERSION = "5.0.9";

    public static void main(String[] args) throws IOException {

        System.out.println(APP_NAME + " " + APP_VERSION);

        long startTime = System.currentTimeMillis();

        // Process command line arguments
        ArgumentParser argumentParser = new ArgumentParser(args);

        // Load parameters from the YAML file
        ConfigLoader configLoader = new ConfigLoader();
        configLoader.loadParameters(argumentParser.getYamlFile());
        if (argumentParser.isValidateOnly()) {
            System.out.println("Validating parameters only");
            return;
        }

        // Load input files
        configLoader.loadFileList(argumentParser.getInputFiles());

        try {
            TMTIntegrator integrator = new TMTIntegrator(configLoader.getParameters());
            integrator.run();
        } catch (Exception e) {
            e.printStackTrace();
        }

        long endTime = System.currentTimeMillis();
        System.out.println("Execution time: " + (endTime - startTime) + " ms");
    }

}
