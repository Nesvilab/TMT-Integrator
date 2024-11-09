package tmtintegrator.pojo.psm;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class PhosphoSiteData {
    int siteLocalCount;
    final StringBuilder siteLocalPos; // site localized position
    final Map<Integer, List<String>> siteLocalMassMap; // Key: position; Value: site localized mass
    final List<String> siteLocalPosList;
    int siteCount;
    final Map<Integer, Double> probMap; // Key: position; Value: site probability
    int startIndex;
    int endIndex;

    PhosphoSiteData() {
        siteLocalCount = 0;
        siteLocalPos = new StringBuilder();
        siteLocalMassMap = new TreeMap<>();
        siteLocalPosList = new ArrayList<>();
        siteCount = 0;
        probMap = new TreeMap<>();
    }
}
