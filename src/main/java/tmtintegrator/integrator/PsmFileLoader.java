package tmtintegrator.integrator;

import tmtintegrator.constants.GroupBy;
import tmtintegrator.pojo.ProcessedPsmEntry;
import tmtintegrator.pojo.Index;
import tmtintegrator.pojo.Parameters;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;

/**
 * Load preprocessed PSM files.
 *
 * @author rogerli on 05/2024
 */
public class PsmFileLoader {

    private final Parameters parameters;

    public PsmFileLoader(Parameters parameter) {
        this.parameters = parameter;
    }

    /**
     * Load preprocessed PSM files, generate PSM entries, and group them by the specified option.
     *
     * @param groupBy group by option
     * @return map of file path to list of PSM entries
     */
    public Map<String, List<String>> loadPsmFiles(GroupBy groupBy) {
        Map<String, List<String>> fileMap = new HashMap<>();
        for (File file : parameters.fileList) {
            List<String> psmList = processFile(file, groupBy);
            fileMap.put(file.getAbsolutePath(), psmList);
        }
        return fileMap;
    }

    private List<String> processFile(File file, GroupBy groupBy) {
        String fileName = file.getAbsolutePath();
        List<String> psmList = new ArrayList<>();

        // get processed PSM file path
        String psmPath = fileName.replace(".tsv", ".ti");
        try (BufferedReader reader = new BufferedReader(new FileReader(psmPath))) {
            String title = reader.readLine();
            parameters.titleMap.put(fileName, title);
            Index index = parameters.indMap.get(fileName);

            // read PSM entries
            String line;
            while ((line = reader.readLine()) != null) {
                ProcessedPsmEntry psmEntry = new ProcessedPsmEntry(parameters, index, line);
                psmEntry.parsePsmEntry();
                if (psmEntry.isUsed()) {
                    psmEntry.analyzePhosphoSites(psmList, groupBy);
                }
            }
        } catch (Exception e) {
            throw new RuntimeException("Error reading file: " + psmPath, e);
        }
        return psmList;
    }

}
