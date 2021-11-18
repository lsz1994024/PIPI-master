package proteomics.Types;


public class ThreeExpAA implements Comparable<ThreeExpAA> {

    private final ExpAA[] threeExpAa;
    private int hashCode;
    private final double totalIntensity;
    private final String ptmFreeAAString;
    private int regionIdx;
    private float[] peakLocations;
    private float[] intensities;

    public ThreeExpAA(ExpAA aa1, ExpAA aa2, ExpAA aa3) {
        threeExpAa = new ExpAA[]{aa1, aa2, aa3};
        peakLocations = new float[]{(float) aa1.getHeadLocation(), (float) aa2.getHeadLocation(), (float) aa3.getHeadLocation(), (float) aa3.getTailLocation()};
        intensities = new float[]{(float) aa1.getHeadIntensity(), (float) aa2.getHeadIntensity(), (float) aa3.getHeadIntensity(), (float) aa3.getTailIntensity()};
        String toString = threeExpAa[0].toString() + "-" + threeExpAa[1].toString() + "-" + threeExpAa[2].toString();
        hashCode = toString.hashCode();

        StringBuilder sb = new StringBuilder(5);
        for (ExpAA aa : threeExpAa) {
            sb.append(aa.getPtmFreeAA());
        }
        ptmFreeAAString = sb.toString();

        double intensity = threeExpAa[0].getHeadIntensity();
        for (ExpAA aa : threeExpAa) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }
    public String join(ThreeExpAA other) { return ptmFreeAAString + other.getPtmFreeAAString();    }
    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        return (other instanceof ThreeExpAA) && (this.hashCode() == other.hashCode());
    }

    public boolean approximateEquals(ThreeExpAA other, double tolerance) {
        for (int i = 0; i < this.size(); ++i) {
            if (!this.get(i).approximateEquals(other.get(i), tolerance)) {
                return false;
            }
        }
        return true;
    }

    public void setTheoLocation(int i, int theoLoc) {
        threeExpAa[i].setTheoLocation(theoLoc);
        // update toString and hashCode
        String toString = threeExpAa[0].toString() + "-" + threeExpAa[1].toString() + "-" + threeExpAa[2].toString();
        hashCode = toString.hashCode();
    }

    private int connects(ThreeExpAA other) {
        int res = -1;
        if (this.threeExpAa[2].getTailLocation() == other.threeExpAa[0].getHeadLocation()) { res = 1;}
        if (this.threeExpAa[1].getTailLocation() == other.threeExpAa[0].getHeadLocation() && this.threeExpAa[2].getTailLocation() == other.threeExpAa[1].getHeadLocation()) { res = 2;}
        if (this.threeExpAa[0].getTailLocation() == other.threeExpAa[0].getHeadLocation() && this.threeExpAa[1].getTailLocation() == other.threeExpAa[1].getHeadLocation() && this.threeExpAa[2].getTailLocation() == other.threeExpAa[2].getHeadLocation()) { res = 3;}
        return res;
    }

    public float connectScore(ThreeExpAA other) {
        float res = 0;
        int connect = this.connects(other);
        if (connect != -1){
            for (int i = 0; i < connect; i++){
                res += other.intensities[i];
            }
        } else {
            res = -1;
        }
        return res;
    }

    public int compareTo(ThreeExpAA other) {
        return Double.compare(threeExpAa[0].getHeadLocation(), other.threeExpAa[0].getHeadLocation());
    }

    public ExpAA[] getExpAAs() {
        return threeExpAa;
    }

    public String getPtmFreeAAString() {
        return ptmFreeAAString;
    }

    public double getTotalIntensity() {
        return totalIntensity;
    }

    public double getHeadLocation() {
        return threeExpAa[0].getHeadLocation();
    }

    public double getTailLocation() {
        return threeExpAa[threeExpAa.length - 1].getTailLocation();
    }

    public ThreeExpAA clone() throws CloneNotSupportedException {
        super.clone();
        return new ThreeExpAA(threeExpAa[0].clone(), threeExpAa[1].clone(), threeExpAa[2].clone());
    }

    public int size() {
        return threeExpAa.length;
    }

    public ExpAA get(int i) {
        return threeExpAa[i];
    }

    public void setRegionIdx(int regionIdx) {
        this.regionIdx = regionIdx;
    }

    public int getRegionIdx() {
        return regionIdx;
    }
}
