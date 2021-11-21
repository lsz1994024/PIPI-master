package proteomics.Segment;


import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import proteomics.Types.*;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class InferSegment {

    private static final int minTagNum = 200;
    private static final int regionNum = 10;
    private static final int topNumInEachRegion = 20;
    private static final Pattern pattern = Pattern.compile("([nc][0-9a-i])?([A-Z#$].?)");

    private final double ms2Tolerance;
    private TreeMap<Segment, Integer> aaVectorTemplate = new TreeMap<>();
    private TreeMap<Segment, Integer> aaConDB = new TreeMap<>();
    private Map<Double, String> modifiedAAMap = new HashMap<>(35, 1);
    private final Double[] deltaMassArray;
    private Map<String, Double> modifiedAAMassMap = new HashMap<>(35, 1);
    private double[] nTermPossibleMod = null;
    private double[] cTermPossibleMod = null;
    private MassTool massTool;
    private Set<String> theoTag3Str = new HashSet<>();

    public InferSegment(MassTool massTool, Map<String, String> parameterMap, Map<Character, Double> fixModMap) throws Exception {
        this.massTool = massTool;
        this.ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        Map<Character, Double> massTable = massTool.getMassTable();

        char[] standardAaArray = new char[]{'G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W', 'U', 'O'};

        Map<Double, Character> massAaMap = new HashMap<>(25, 1);
        for (char aa : standardAaArray) {
            // # = I/L.
            if (aa == 'I' || aa == 'L') {
                massAaMap.put(massTable.get(aa), '#');
            } else {
                massAaMap.put(massTable.get(aa), aa);
            }
        }

        Character[] aaArray = massAaMap.values().toArray(new Character[0]);

        for (char aa1 : aaArray) {
            for (char aa2 : aaArray) {
                for (char aa3 : aaArray) {
                    theoTag3Str.add(String.format(Locale.US, "%c%c%c", aa1, aa2, aa3));
                }
            }
        }

        for (String tag1 : theoTag3Str) {
            Segment tag1_seg = new Segment(tag1, false);
            if (!aaConDB.containsKey(tag1_seg)) {
                aaConDB.put(tag1_seg, 0);
            }
            for (String tag2 : theoTag3Str) {
                String tag_tag = tag1+tag2;
                Segment tag_tag_seg = new Segment(tag_tag, false);
                if (!aaConDB.containsKey(tag_tag_seg)) {
                    aaConDB.put(tag_tag_seg, 0);
                }
            }
        }

        for (char aa1 : aaArray) {
            for (char aa2 : aaArray) {
                for (char aa3 : aaArray) {
                    aaVectorTemplate.put(new Segment(String.format(Locale.US, "%c%c%c", aa1, aa2, aa3)), 0);
                }
            }
        }

        int idx = 0;
        for (Segment segment : aaVectorTemplate.keySet()) {
            aaVectorTemplate.put(segment, idx);
            ++idx;
        }


        idx = 0;
        for (Segment segment : aaConDB.keySet()) {
            aaConDB.put(segment, idx);
            ++idx;
        }

        // generate a mass aa map containing modified amino acid
        for (double k : massAaMap.keySet()) {
            modifiedAAMap.put(k, massAaMap.get(k).toString());
//            System.out.println(k +" "+ massAaMap.get(k).toString());
        }
        for (String k : parameterMap.keySet()) {
            if (k.startsWith("mod")) {
                String v = parameterMap.get(k);
                if (!v.startsWith("0.0")) {
                    String[] temp = v.split("@");
                    double tempMass = massTable.get(temp[1].charAt(0)) + Double.valueOf(temp[0]);
                    // check if the mass has conflict
                    for (double temp2 : modifiedAAMap.keySet()) {
                        if (Math.abs(temp2 - tempMass) <= ms2Tolerance) {
                            throw new Exception(String.format(Locale.US, "%s and %s have conflict mass values(%f vs %f).", v, modifiedAAMap.get(temp2), tempMass, temp2));
                        }
                    }
                    if (Math.abs(fixModMap.get(temp[1].charAt(0))) < 0.1) {
                        // fix modification and var modification cannot be coexist
                        if ((temp[1].charAt(0) == 'I') || (temp[1].charAt(0) == 'L')) {
                            modifiedAAMap.put(tempMass, temp[1].replace(temp[1].charAt(0), '#'));
                            modifiedAAMassMap.put(temp[1].replace(temp[1].charAt(0), '#'), Double.valueOf(temp[0]));
                        } else {
                            modifiedAAMap.put(tempMass, temp[1]);
                            modifiedAAMassMap.put(temp[1], Double.valueOf(temp[0]));
                        }
                    }
                }
            } else if (k.contentEquals("Nterm")) {
                if (Math.abs(fixModMap.get('n')) < 0.1) {
                    // fix modification and var modification cannot be coexist
                    if (!parameterMap.get(k).startsWith("0.0")) {
                        String[] tempArray = parameterMap.get(k).split(",");
                        nTermPossibleMod = new double[tempArray.length];
                        for (int i = 0; i < tempArray.length; ++i) {
                            nTermPossibleMod[i] = Double.valueOf(tempArray[i].trim());
                        }
                    }
                }
            } else if (k.contentEquals("Cterm")) {
                if (Math.abs(fixModMap.get('c')) < 0.1) {
                    // fix modification and var modification cannot be coexist
                    if (!parameterMap.get(k).startsWith("0.0")) {
                        String[] tempArray = parameterMap.get(k).split(",");
                        cTermPossibleMod = new double[tempArray.length];
                        for (int i = 0; i < tempArray.length; ++i) {
                            cTermPossibleMod[i] = Double.valueOf(tempArray[i].trim());
                        }
                    }
                }
            }
        }
        deltaMassArray = modifiedAAMap.keySet().toArray(new Double[0]);
    }

    public List<ThreeExpAA> inferSegmentLocationFromSpectrum(double precursorMass, TreeMap<Double, Double> plMap) throws Exception {
        return inferThreeAAFromSpectrum(addVirtualPeaks(precursorMass, plMap), precursorMass - massTool.H2O + MassTool.PROTON);
    }

    public SparseVector generateSegmentIntensityVector(List<ThreeExpAA> inputList) {
        SparseVector finalVector = new SparseVector();
        if (inputList.isEmpty()) {
            return finalVector;
        } else {
            for (ThreeExpAA expAaList : inputList) {
                double totalIntensity = expAaList.getTotalIntensity();
                int idx = aaVectorTemplate.get(new Segment(expAaList.getPtmFreeAAString()));
                double value = Math.max(totalIntensity, finalVector.get(idx));
                finalVector.put(idx, value);
            }
            return finalVector;
        }
    }

    public SparseBooleanVector generateSegmentIntensityVectorSBV(List<ThreeExpAA> inputList) {
        SparseBooleanVector finalVector = new SparseBooleanVector();
        if (inputList.isEmpty()) {
            return finalVector;
        } else {
            for (ThreeExpAA expAaList : inputList) {
                double totalIntensity = expAaList.getTotalIntensity();
                int idx = aaVectorTemplate.get(new Segment(expAaList.getPtmFreeAAString()));
                double value = Math.max(totalIntensity, finalVector.get(idx));
                finalVector.put(idx, value);
            }
            return finalVector;
        }
    }

    public SparseBooleanVector generateExpSegMat(List<ThreeExpAA> inputList) {
        SparseBooleanVector finalVector = new SparseBooleanVector();
        if (!inputList.isEmpty()) {
            int idx = 0;
            for (ThreeExpAA tag1 : inputList) {
                idx = aaConDB.get(new Segment(tag1.getPtmFreeAAString(), false));
                float tag1Score = 0;
                boolean isNCtag = false;
                if (Math.abs(tag1.getHeadLocation() - massTool.PROTON) <= ms2Tolerance
                        ||Math.abs(tag1.getHeadLocation() - massTool.PROTON - massTool.H2O) <= ms2Tolerance ){
                    tag1Score = (float) (2.3*tag1.getTotalIntensity());
                    isNCtag = true;
                } else {
                    tag1Score = (float) (1*tag1.getTotalIntensity());
                }

                finalVector.put(idx, Math.max(tag1Score, finalVector.get(idx)));
                for (ThreeExpAA tag2 : inputList) {
                    float connectScore = tag1.connectScore(tag2);
                    if (connectScore == -1) {
                        continue;
                    }
//                    System.out.print(tag1 + "+"+tag2+"+"+connectScore+"+"+"      ");
                    float tag2Score = 0;
                    if (isNCtag) {
                        tag2Score = (3*connectScore);
                    } else {
                        tag2Score = 1*connectScore;
                    }
                    idx = aaConDB.get(new Segment(tag1.join(tag2), false)); // use the idx of tag string to represent tag
                    finalVector.put(idx, (float) (1.5*Math.max(tag2Score, finalVector.get(idx))));
                }
            }
        }
        return finalVector;
    }

    public SparseBooleanVector generateTheoSegMat(String peptide) throws Exception {
        TreeMap<Double, Double> theoPeakMap = generateTheoPeak(peptide);
        List<ThreeExpAA> theoTags = inferThreeAAFromSpectrum(theoPeakMap, calResMass(peptide) + MassTool.PROTON);
        if (peptide.equals("HVGYKPSDEHK")){
            System.out.println("Thoe tags");
        }
        SparseBooleanVector finalVector = new SparseBooleanVector();
        if (!theoTags.isEmpty()) {
            int idx = 0;
            for (ThreeExpAA tag1 : theoTags) {
                idx = aaConDB.get(new Segment(tag1.getPtmFreeAAString(), false));
                finalVector.put(idx, 1);
                for (ThreeExpAA tag2 : theoTags) {
                    float connectScore = tag1.connectScore(tag2);
                    float localScore = 0;
                    if (connectScore == -1) {
                        continue;
                    } else if (connectScore == 3){
                        localScore = (float) 1.0;
                    } else if (connectScore == 2){
                        localScore = (float) 0.8;
                    } else if (connectScore == 1) {
                        localScore = (float) 0.4;
                    }
//                    System.out.print("HERE local");
                    idx = aaConDB.get(new Segment(tag1.join(tag2), false)); // use the idx of tag string to represent tag
                    finalVector.put(idx, (float) (1*Math.max(localScore, finalVector.get(idx))));
                }
            }
        }
        return finalVector;


    }

    public Double calResMass(String seq){ // Caution: nterm mod is not included!
        Double totalMass = 0.0;
        int length = seq.length();
        for (int idx = 0; idx < length; ++idx) {
//            if (! massTool.getMassTable().containsKey(seq.substring(idx, idx + 1)))
//            {
//                System.out.println(seq+" "+seq.substring(idx, idx + 1));
//            }
            totalMass += massTool.getMassTable().get(seq.charAt(idx));
        }

        return totalMass;
    }

    public TreeMap<Double, Double> generateTheoPeak(String pep){
        String peptide = normalizeSequence(pep);
        TreeMap<Double, Double> peakList = new TreeMap<>();
        int pepLen = peptide.length();
        Double intes = 1.0;
        for (int i = 1; i <= pepLen - 1; i++){
            String bIon = peptide.substring(0, i);
            Double bMass = calResMass(bIon);
            peakList.put(bMass, intes);

            String yIon = peptide.substring(i); // to the end
            Double yMass = calResMass(yIon) + massTool.H2O;
            peakList.put(yMass, intes);
        }
        peakList.put(calResMass(peptide), intes);
        peakList.put(massTool.H2O, intes);
        peakList.put(0.0, intes);
//        peakList.put(calResMass(peptide) - massTable.get("H2O"), intes);
        return peakList;
    }

    public SparseBooleanVector generateSegmentBooleanVector(String peptide) {
        String normalizedPeptide = normalizeSequence(peptide);
        Set<Integer> tempSet = new HashSet<>(peptide.length() + 1, 1);
        HashMap<Integer, Double> tempMap = new HashMap<>(peptide.length() + 1, 1);
        int tempKey;
        double tempValue;
        for (int i = 0; i <= normalizedPeptide.length() - 3; ++i) {
            tempSet.add(aaVectorTemplate.get(new Segment(normalizedPeptide.substring(i, i + 3))));
            tempKey = aaVectorTemplate.get(new Segment(normalizedPeptide.substring(i, i + 3)));
            if (tempMap.containsKey(tempKey)){
                tempValue = tempMap.get(tempKey) + 1.0;
                tempMap.put(tempKey, tempValue);
            }
            else{
                tempMap.put(tempKey, 1.0);
            }
        }
        return new SparseBooleanVector(tempMap);
    }

    public static String normalizeSequence(String seq) {
        return seq.replaceAll("[IL]", "#");
    }

    public List<ThreeExpAA> inferThreeAAFromSpectrum(TreeMap<Double, Double> plMap, double cTermMz) throws Exception {
        Double[] mzArray = plMap.keySet().toArray(new Double[0]);
        Double[] intensityArray = plMap.values().toArray(new Double[0]);
        Set<ThreeExpAA> tempSet = new HashSet<>();
        List<ThreeExpAA> outputList = new LinkedList<>();
        for (int i = 0; i < mzArray.length - 3; ++i) {
            double mz1 = mzArray[i];
            double intensity1 = intensityArray[i];
            for (int j = i + 1; j < mzArray.length - 2; ++j) {
                double mz2 = mzArray[j];
                double intensity2 = intensityArray[j];
                String aa1 = inferAA(mz1, mz2, Math.abs(mz1 - MassTool.PROTON) <= ms2Tolerance, false);
                if (aa1 != null) {
                    Matcher matcher = pattern.matcher(aa1);
                    char ptmFreeAA = '\0';
                    double mod = 0;
                    double nTermMod = 0;
                    if (matcher.matches()) {
                        if (modifiedAAMassMap.containsKey(matcher.group(2))) {
                            mod = modifiedAAMassMap.get(matcher.group(2));
                        }
                        ptmFreeAA = matcher.group(2).charAt(0);
                        if (matcher.group(1) != null) {
                            if ((matcher.group(1).charAt(1) - '0' >= 0) && (matcher.group(1).charAt(1) - '0' < 10)) {
                                nTermMod = nTermPossibleMod[matcher.group(1).charAt(1) - '0'];
                            } else {
                                throw new Exception("Something is wrong in inferring tags.");
                            }
                        }
                    } else {
                        throw new NullPointerException(String.format(Locale.US, "Cannot find the PTM free amino acid for %s.", aa1));
                    }
                    ExpAA expAa1 = new ExpAA(aa1, ptmFreeAA, mz1, mz2, intensity1, intensity2, -1, mod, nTermMod, 0);
                    List<List<ExpAA>> tempAasList2 = new LinkedList<>();
                    for (int k = j + 1; k < mzArray.length - 1; ++k) {
                        double mz3 = mzArray[k];
                        double intensity3 = intensityArray[k];
                        String aa2 = inferAA(mz2, mz3, false, false);
                        if (aa2 != null) {
                            mod = 0;
                            if (modifiedAAMassMap.containsKey(aa2)) {
                                mod = modifiedAAMassMap.get(aa2);
                            }
                            ExpAA expAa2 = new ExpAA(aa2, aa2.charAt(0), mz2, mz3, intensity2, intensity3, -1, mod, 0, 0);
                            List<ExpAA> tempAasList3 = new LinkedList<>();
                            for (int l = k + 1; l < mzArray.length; ++l) {
                                double mz4 = mzArray[l];
                                double intensity4 = intensityArray[l];
                                String aa3 = inferAA(mz3, mz4, false, Math.abs(mz4 - cTermMz) <= ms2Tolerance);
                                if (aa3 != null) {
                                    matcher = pattern.matcher(aa3);
                                    ptmFreeAA = '\0';
                                    mod = 0;
                                    double cTermMod = 0;
                                    if (matcher.matches()) {
                                        if (modifiedAAMassMap.containsKey(matcher.group(2))) {
                                            mod = modifiedAAMassMap.get(matcher.group(2));
                                        }
                                        ptmFreeAA = matcher.group(2).charAt(0);
                                        if (matcher.group(1) != null) {
                                            if ((matcher.group(1).charAt(1) - '0' >= 0) && (matcher.group(1).charAt(1) - '0' < 10)) {
                                                cTermMod = cTermPossibleMod[matcher.group(1).charAt(1) - '0'];
                                            } else {
                                                throw new Exception("Something is wrong in inferring tags.");
                                            }
                                        }
                                    } else {
                                        throw new NullPointerException(String.format(Locale.US, "Cannot find the PTM free amino acid for %s.", aa3));
                                    }
                                    ExpAA expAa3 = new ExpAA(aa3, ptmFreeAA, mz3, mz4, intensity3, intensity4, -1, mod, 0, cTermMod);
                                    tempAasList3.add(expAa3);
                                }
                            }
                            for (ExpAA expAas3 : tempAasList3) {
                                List<ExpAA> tempList2 = new LinkedList<>();
                                tempList2.add(expAa2);
                                tempList2.add(expAas3);
                                tempAasList2.add(tempList2);
                            }
                        }
                    }
                    for (List<ExpAA> expAas2 : tempAasList2) {
                        ThreeExpAA threeExpAa = new ThreeExpAA(expAa1, expAas2.get(0), expAas2.get(1));
                        tempSet.add(threeExpAa);
                    }
                }
            }
        }

        // eliminate "overlapped" tags
        ThreeExpAA[] tempArray = tempSet.toArray(new ThreeExpAA[0]);
        List<ThreeExpAA> tempList = new LinkedList<>();
        for (int i = 0; i < tempArray.length; ++i) {
            boolean keep = true;
            for (int j = 0; j < tempArray.length; ++j) {
                if (i != j) {
                    if (tempArray[i].approximateEquals(tempArray[j], 2 * ms2Tolerance)) {
                        if (tempArray[i].getTotalIntensity() < tempArray[j].getTotalIntensity()) {
                            keep = false;
                            break;
                        }
                    }
                }
            }
            if (keep) {
                tempList.add(tempArray[i]);
            }
        }

        if (tempList.size() > minTagNum) {
            double minMz = plMap.firstKey();
            double regionWindow = Math.ceil((plMap.lastKey() - minMz) / regionNum);
            for (ThreeExpAA expAa : tempList) {
                expAa.setRegionIdx((int) Math.floor((expAa.getHeadLocation() - minMz) / regionWindow));
            }
            List<List<Double>> regionIntensityList = new ArrayList<>(20);
            for (int i = 0; i < regionNum; ++i) {
                regionIntensityList.add(new ArrayList<>(100));
            }
            for (ThreeExpAA expAa : tempList) {
                regionIntensityList.get(expAa.getRegionIdx()).add(expAa.getTotalIntensity());
            }
            double[] intensityTArray = new double[regionNum];
            for (int i = 0; i < regionNum; ++i) {
                List<Double> intensityList = regionIntensityList.get(i);
                Collections.sort(intensityList, Collections.reverseOrder());
                if (intensityList.size() > topNumInEachRegion) {
                    intensityTArray[i] = intensityList.get(topNumInEachRegion);
                }
            }
            for (ThreeExpAA expAa : tempList) {
                if (expAa.getTotalIntensity() > intensityTArray[expAa.getRegionIdx()]) {
                    outputList.add(expAa);
                }
            }
            return outputList;
        } else {
            return tempList;
        }
    }

    private String inferAA(double mz1, double mz2, boolean nTerm, boolean cTerm) {
        double mzDiff = mz2 - mz1;
        for (double mass : deltaMassArray) {
            if (Math.abs(mzDiff - mass) <= 2 * ms2Tolerance) {
                return modifiedAAMap.get(mass);
            }
        }

        if (nTerm && (nTermPossibleMod != null)) {
            for (double mass : deltaMassArray) {
                for (int i = 0; i < nTermPossibleMod.length; ++i) {
                    if (Math.abs(mzDiff - mass - nTermPossibleMod[i]) <= 2 * ms2Tolerance) {
                        return "n" + i + modifiedAAMap.get(mass);
                    }
                }
            }
        }

        if (cTerm && (cTermPossibleMod != null)) {
            for (double mass : deltaMassArray) {
                for (int i = 0; i < cTermPossibleMod.length; ++i) {
                    if (Math.abs(mzDiff - mass - cTermPossibleMod[i]) <= 2 * ms2Tolerance) {
                        return "c" + i + modifiedAAMap.get(mass);
                    }
                }
            }
        }

        return null;
    }

    private TreeMap<Double, Double> addVirtualPeaks(double precursorMass, TreeMap<Double, Double> plMap) {
        double totalMass = precursorMass + 2 * MassTool.PROTON;
        TreeMap<Double, Double> finalPlMap = new TreeMap<>();
        for (double mz : plMap.keySet()) {
            finalPlMap.put(mz, plMap.get(mz));
        }
        for (double mz : plMap.keySet()) {
            double anotherMz = totalMass - mz;
            double leftMz = anotherMz - ms2Tolerance;
            double rightMz = anotherMz + ms2Tolerance;
            NavigableMap<Double, Double> temp = null;
            try {
                temp = plMap.subMap(leftMz, true, rightMz, true);
            } catch (IllegalArgumentException ex) {}

            if ((temp == null) || (temp.isEmpty())) {
                finalPlMap.put(anotherMz, plMap.get(mz));
            }
        }

        // Add two virtual peak. Because we have convert all y-ions to b-ions.
        finalPlMap.put(MassTool.PROTON, 1d);
        double cTermMz = precursorMass - massTool.H2O + MassTool.PROTON;
        double leftMz = cTermMz - ms2Tolerance;
        double rightMz = cTermMz + ms2Tolerance;
        NavigableMap<Double, Double> temp = null;
        try {
            temp = plMap.subMap(leftMz, true, rightMz, true);
        } catch (IllegalArgumentException ex) {}
        if ((temp == null) || (temp.isEmpty())) {
            finalPlMap.put(cTermMz, 1d);
        }

        return finalPlMap;
    }
}
