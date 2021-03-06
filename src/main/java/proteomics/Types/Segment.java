package proteomics.Types;


public class Segment implements Comparable<Segment> {

    private final String segmentString;
    private final String revertedVersion;
    private final int hashCode;

    public Segment(String segmentString) {
        String temp = reverseString(segmentString);
        // Compare the string and its reversed version. Kept the smaller one.
        int compareResult = segmentString.compareTo(temp);
        if (compareResult > 0) {
            this.segmentString = temp;
            revertedVersion = segmentString;
        } else {
            this.segmentString = segmentString;
            revertedVersion = temp;
        }
        hashCode = segmentString.hashCode();
    }

    public Segment(String segmentString, boolean shouldDiff) {

        if (shouldDiff) {
            this.segmentString = segmentString;
            this.revertedVersion = reverseString(segmentString);
            hashCode = this.segmentString.hashCode();
        } else {
            String temp = reverseString(segmentString);
            // Compare the string and its reversed version. Kept the smaller one.
            int compareResult = segmentString.compareTo(temp);
            if (compareResult > 0) {
                this.segmentString = temp;
                revertedVersion = segmentString;
            } else {
                this.segmentString = segmentString;
                revertedVersion = temp;
            }
            hashCode = this.segmentString.hashCode(); // this or not this has no difference
        }
    }

    public int length() {
        return segmentString.length();
    }

    public  boolean equals(Object other) {
        return (other instanceof Segment) && this.segmentString.contentEquals(((Segment) other).segmentString);
    }

    public String toString() {
        return segmentString;
    }

    public int hashCode() {
        return hashCode;
    }

    public int compareTo(Segment other) {
        return this.segmentString.compareTo(other.segmentString);
    }

    private String reverseString(String seq) {
        return new StringBuilder(seq).reverse().toString();
    }
}
