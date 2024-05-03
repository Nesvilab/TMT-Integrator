package tmtintegrator;

import tmtintegrator.config.ArgumentParser;
import tmtintegrator.config.ConfigLoader;

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
        long startLoadTime = System.currentTimeMillis();
        ConfigLoader configLoader = new ConfigLoader();
        configLoader.loadParameters(argumentParser.getYamlFile());
        if (argumentParser.isValidateOnly()) {
            System.out.println("Validating parameters only");
            return;
        }

        // Load input files
        configLoader.loadFileList(argumentParser.getInputFiles());
        long endLoadTime = System.currentTimeMillis();
        System.out.println("Parameter Loading: " + (endLoadTime - startLoadTime) + " ms");

        try {
            TMTIntegrator integrator = new TMTIntegrator(configLoader.getParameters());
            integrator.run();
        } catch (Exception e) {
            // TODO: handle exception, log error
            e.printStackTrace();
        }

        long endTime = System.currentTimeMillis();
        System.out.println("Execution time: " + (endTime - startTime) + " ms");
    }

}
