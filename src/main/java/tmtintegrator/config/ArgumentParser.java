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

package tmtintegrator.config;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Parses command line arguments.
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
