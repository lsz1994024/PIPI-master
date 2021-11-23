package proteomics.Search;

import proteomics.Index.BuildIndex;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import proteomics.Types.*;

import java.util.*;

public class Search {

    private static final int rankNum = 5;

    private List<Peptide> ptmOnlyResult = new LinkedList<>();
    private List<Peptide> ptmFreeResult = new LinkedList<>();
    private List<Peptide> topPepList = new LinkedList<>();
    public int numTieTop1Peps = 0;
    public Set<String> backboneList = new HashSet<>();
    public double topScore = -20d;
    public Search(BuildIndex buildIndex, double precursorMass, int scanNum, SparseBooleanVector scanCode, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double minPtmMass, double maxPtmMass, int localMaxMs2Charge) {
        PriorityQueue<ResultEntry> ptmFreeQueue = new PriorityQueue<>(rankNum * 2);
        PriorityQueue<ResultEntry> ptmOnlyQueue = new PriorityQueue<>(rankNum * 2);
        double scanNormSquare = scanCode.norm2square();
        double leftTol = ms1Tolerance;
        double rightTol = ms1Tolerance;
        if (ms1ToleranceUnit == 1) {
            leftTol = precursorMass - (precursorMass * leftInverseMs1Tolerance);
            rightTol = (precursorMass * rightInverseMs1Tolerance) - precursorMass;
        }
        double leftMass = Math.max(precursorMass + minPtmMass - leftTol, buildIndex.getMinPeptideMass());
        double rightMass = Math.min(precursorMass + maxPtmMass + rightTol, buildIndex.getMaxPeptideMass());

        if (leftMass >= rightMass) {
            return;
        }

        Map<String, Peptide0> peptide0Map = buildIndex.getPeptide0Map();
        TreeMap<Double, Set<String>> massPeptideMap = buildIndex.getMassPeptideMap();

        NavigableMap<Double, Set<String>> subMassPeptideMap = massPeptideMap.subMap(leftMass, true, rightMass, true);

        if (!subMassPeptideMap.isEmpty()) {
            for (double mass : subMassPeptideMap.keySet()) {
                for (String sequence : massPeptideMap.get(mass)) {
                    Peptide0 peptide0 = peptide0Map.get(sequence);
                    double score = 0;
                    double temp1 = Math.sqrt(peptide0.code.norm2square() * scanNormSquare);
                    if (temp1 > 1e-6) {
                        score = peptide0.code.dot(scanCode) / temp1;
//                        if (scanNum == 2132 && sequence.equals("nHPTTPSSTAKc")){
//                            System.out.println("in score 2132 " + score);
//                            System.out.println("pep "+ peptide0.code);
//                            System.out.println("scan "+ scanCode.sparseVector);
//                        }
                    }
                    double deltaMass = mass - precursorMass; // caution: the order matters under ms1ToleranceUnit == 1 situation

                    if (peptide0.isTarget) {
                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                            // PTM-free
                            if (ptmFreeQueue.size() < rankNum) {
                                ptmFreeQueue.add(new ResultEntry(score, sequence, false));
                            } else {
                                if (score > ptmFreeQueue.peek().score) {
                                    ptmFreeQueue.poll();
                                    ptmFreeQueue.add(new ResultEntry(score, sequence, false));
                                }
                            }
                        }

                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                            // PTM-only
                            if (ptmOnlyQueue.size() < rankNum) {
                                ptmOnlyQueue.add(new ResultEntry(score, sequence, false));
                            } else {
                                if (score > ptmOnlyQueue.peek().score) {
                                    ptmOnlyQueue.poll();
                                    ptmOnlyQueue.add(new ResultEntry(score, sequence, false));
                                }
                            }
                        }
                    } else {
                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                            // PTM-free
                            if (ptmFreeQueue.size() < rankNum) {
                                ptmFreeQueue.add(new ResultEntry(score, sequence, true));
                            } else {
                                if (score > ptmFreeQueue.peek().score) {
                                    ptmFreeQueue.poll();
                                    ptmFreeQueue.add(new ResultEntry(score, sequence, true));
                                }
                            }
                        }

                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                            // PTM-only
                            if (ptmOnlyQueue.size() < rankNum) {
                                ptmOnlyQueue.add(new ResultEntry(score, sequence, true));
                            } else {
                                if (score > ptmOnlyQueue.peek().score) {
                                    ptmOnlyQueue.poll();
                                    ptmOnlyQueue.add(new ResultEntry(score, sequence, true));
                                }
                            }
                        }
                    }
                }
            }
        }

        if (!(ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty())) {
            topPepList = mergeResult(ptmFreeQueue, ptmOnlyQueue, massTool, localMaxMs2Charge);
        }


        if (!(ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty())) {
            ptmFreeResult = convertResult(ptmFreeQueue, massTool, localMaxMs2Charge);
            ptmOnlyResult = convertResult(ptmOnlyQueue, massTool, localMaxMs2Charge);
        }
    }

    private List<Peptide> mergeResult(PriorityQueue<ResultEntry> ptmFreeQueue, PriorityQueue<ResultEntry> ptmOnlyQueue, MassTool massTool, int localMaxMs2Charge) {
        List<Peptide> peptideList = new LinkedList<>();
        int freeSize = ptmFreeQueue.size();
        int onlySize = ptmOnlyQueue.size();
        if (!ptmFreeQueue.isEmpty()) {
            for (ResultEntry temp : ptmFreeQueue){
                if (temp.score > topScore){
                    topScore = temp.score;
                    numTieTop1Peps = 1;
                    backboneList.clear();
                    backboneList.add(temp.peptide);
                } else if(temp.score == topScore) {
                    numTieTop1Peps += 1;
                    backboneList.add(temp.peptide);
                }
                peptideList.add(new Peptide(temp.peptide, temp.isDecoy(), massTool, localMaxMs2Charge, temp.score, freeSize));
            }
        }

        if (!ptmOnlyQueue.isEmpty()) {
            for (ResultEntry temp : ptmOnlyQueue){
                if (temp.score > topScore){
                    topScore = temp.score;
                    numTieTop1Peps = 1;
                    backboneList.clear();
                    backboneList.add(temp.peptide);
                } else if(temp.score == topScore) {
                    numTieTop1Peps += 1;
                    backboneList.add(temp.peptide);
                }
                peptideList.add(new Peptide(temp.peptide, temp.isDecoy(), massTool, localMaxMs2Charge, temp.score, onlySize));
            }
        }

        return peptideList;
    }

    private List<Peptide> convertResult(PriorityQueue<ResultEntry> inputQueue, MassTool massTool, int localMaxMs2Charge) {
        List<Peptide> peptideList = new LinkedList<>();
        int globalRank = inputQueue.size();
        while (!inputQueue.isEmpty()) {
            ResultEntry temp = inputQueue.poll();
            peptideList.add(new Peptide(temp.peptide, temp.isDecoy(), massTool, localMaxMs2Charge, temp.score, globalRank));
            --globalRank;
        }

        return peptideList;
    }

    public List<Peptide> getPTMOnlyResult() {
        return ptmOnlyResult;
    }

    public List<Peptide> getPTMFreeResult() {
        return ptmFreeResult;
    }

    public List<Peptide> getAllPepList() {
        return topPepList;
    }
}
