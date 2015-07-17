/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.core;

import edu.duke.humangenome.sam.MultiPileup;
import edu.duke.humangenome.sam.Pileup;
import edu.duke.humangenome.util.BaseUtils;
import java.util.Map;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.SamLocusIterator;

/**
 *
 * @author Yongzhuang Liu
 */
public class FeatureSelection {

    private final static int FATHER_INDEX = 0;
    private final static int MOTHER_INDEX = 1;
    private final static int OFFSPRING_INDEX = 2;
    private ReferenceSequence referenceSequence;
    private MultiPileup trioPileup;

    public FeatureSelection(ReferenceSequence referenceSequence, MultiPileup trioPileup) {
        this.referenceSequence = referenceSequence;
        this.trioPileup = trioPileup;
    }

    public String extract() {

        String chrom = trioPileup.getReferenceName();
        int pos = trioPileup.getPosition();
        char ref = getRefBase();
        char alt = getMostAltAllele();
        int[] index = new int[]{FATHER_INDEX, MOTHER_INDEX, OFFSPRING_INDEX};

        StringBuilder recordBuilder = new StringBuilder();

        int[] fAlleleCount = new int[2];
        int[] mAlleleCount = new int[2];
        int[] oAlleleCount = new int[2];

        for (int i = 0; i < index.length; i++) {

            Pileup pileup = new Pileup(trioPileup.getLocusInfo(index[i]));

            int readDepth = pileup.getDepth();

            int meanRefBaseQuality = pileup.getMeanBaseQuality(ref);
            int meanAltBaseQuality = pileup.getMeanBaseQuality(alt);

            int meanRefMappingQuality = pileup.getMeanMappingQuality(ref);
            int meanAltMappingQuality = pileup.getMeanMappingQuality(alt);

            int meanRefDistanceToThreePrime = pileup.getMeanDistanceToThreePrime(ref);
            int meanAltDistanceToThreePrime = pileup.getMeanDistanceToThreePrime(alt);

            double refFractionOfMQ0Reads = pileup.getFractionOfMQ0Reads(ref);
            double altFractionOfMQ0Reads = pileup.getFractionOfMQ0Reads(alt);

            double refFractionOfSoftClippedReads = pileup.getFractionOfSoftClippedReads(ref);
            double altFractionOfSoftClippedReads = pileup.getFractionOfSoftClippedReads(alt);

            int forwardRef = pileup.getForwardBaseCount(ref);
            int reverseRef = pileup.getReverseBaseCount(ref);
            int forwardAlt = pileup.getForwardBaseCount(alt);
            int reverseAlt = pileup.getReverseBaseCount(alt);

            int strandBias = getStrandBias(forwardRef, reverseRef, forwardAlt, reverseAlt);

            int refStrandDirection;
            int altStrandDirection;
            if (getStrandDirection(forwardRef, reverseRef)) {
                refStrandDirection = 1;
            } else {
                refStrandDirection = 0;
            }
            if (getStrandDirection(forwardAlt, reverseAlt)) {
                altStrandDirection = 1;
            } else {
                altStrandDirection = 0;
            }

            if (index[i] == FATHER_INDEX) {
                fAlleleCount[0] = forwardRef + reverseRef;
                fAlleleCount[1] = forwardAlt + reverseAlt;
            }
            if (index[i] == MOTHER_INDEX) {
                mAlleleCount[0] = forwardRef + reverseRef;
                mAlleleCount[1] = forwardAlt + reverseAlt;
            }
            if (index[i] == OFFSPRING_INDEX) {
                oAlleleCount[0] = forwardRef + reverseRef;
                oAlleleCount[1] = forwardAlt + reverseAlt;
            }

            double alleleBalance = (double) (forwardAlt + reverseAlt) / readDepth;

/*
            double meanRefNearbyIndels = (double) pileup.getNearbyIndels(ref) / (forwardRef + reverseRef);
            double meanAltNearbyIndels = (double) pileup.getNearbyIndels(alt) /  (forwardAlt+ reverseAlt);
            double meanRefNearbyMismatches = (double) (pileup.getOverallMismatches(referenceSequence, ref)) / (forwardRef + reverseRef);
            double meanAltNearbyMismatches = (double) (pileup.getOverallMismatches(referenceSequence, alt) - forwardAlt - reverseAlt) / (forwardAlt+ reverseAlt);

 * 
 */
            
            
            double meanRefNearbyIndels = divide(pileup.getNearbyIndels(ref),(forwardRef + reverseRef));
            double meanAltNearbyIndels = divide(pileup.getNearbyIndels(alt),(forwardAlt+ reverseAlt));
            double meanRefNearbyMismatches = divide(pileup.getOverallMismatches(referenceSequence, ref),(forwardRef + reverseRef));
            double meanAltNearbyMismatches = divide((pileup.getOverallMismatches(referenceSequence, alt) - forwardAlt - reverseAlt),(forwardAlt+ reverseAlt));
            
            
            
            
            
            recordBuilder.append(alleleBalance + "," + meanRefBaseQuality + "," + meanAltBaseQuality + "," + readDepth + ",");
            recordBuilder.append(meanRefMappingQuality + "," + meanAltMappingQuality + ",");
            recordBuilder.append(meanRefDistanceToThreePrime + "," + meanAltDistanceToThreePrime + ",");
            recordBuilder.append(refFractionOfMQ0Reads + "," + altFractionOfMQ0Reads + ",");
            recordBuilder.append(refFractionOfSoftClippedReads + "," + altFractionOfSoftClippedReads + ",");
            recordBuilder.append(meanRefNearbyMismatches + "," + meanAltNearbyMismatches + ",");
            recordBuilder.append(meanRefNearbyIndels + "," + meanAltNearbyIndels + ",");
            recordBuilder.append(refStrandDirection + "," + altStrandDirection + ",");
            recordBuilder.append(strandBias + ",");
        }
        int pValueOfFatherToOffspring = getPhredPValue(fAlleleCount, oAlleleCount);
        int pValueOfMotherToOffspring = getPhredPValue(mAlleleCount, oAlleleCount);
        recordBuilder.append(pValueOfFatherToOffspring + "," + pValueOfMotherToOffspring);
        return recordBuilder.toString();
    }

    public char getMostAltAllele() {
        int[] base = new int[4];
        Map<Integer, SamLocusIterator.LocusInfo> map = trioPileup.getLocusInfoMap();
        for (Map.Entry<Integer, SamLocusIterator.LocusInfo> entry : map.entrySet()) {
            Pileup pileup = new Pileup(entry.getValue());
            base[BaseUtils.baseToIndex('A')] += pileup.getBothStrandBaseCount('A');
            base[BaseUtils.baseToIndex('C')] += pileup.getBothStrandBaseCount('C');
            base[BaseUtils.baseToIndex('G')] += pileup.getBothStrandBaseCount('G');
            base[BaseUtils.baseToIndex('T')] += pileup.getBothStrandBaseCount('T');
        }
        int maxIndex = -1;
        int maxNumber = 0;
        int refIndex = BaseUtils.baseToIndex(getRefBase());
        for (int i = 0; i < 4; i++) {
            if (i != refIndex) {
                if (base[i] > maxNumber) {
                    maxIndex = i;
                    maxNumber = base[i];
                }
            }
        }
        return BaseUtils.indexToBase(maxIndex);
    }

    public char getRefBase() {
        return Character.toUpperCase((char) referenceSequence.getBases()[trioPileup.getPosition() - 1]);
    }

    public int getStrandBias(int forwardRef, int reverseRef, int forwardAlt, int reverseAlt) {
        FisherExact fe = new FisherExact(forwardRef + reverseRef + forwardAlt + reverseAlt);
        double sb = fe.getCumlativeP(forwardRef, reverseRef, forwardAlt, reverseAlt);
        if (sb > 1) {
            sb = 1.0;
        }
        return (int) (Math.log10(sb) * (-10));
    }

    public boolean getStrandDirection(int forwardAllele, int reverseAllele) {
        if (Math.min(forwardAllele, reverseAllele) == 0) {
            return false;
        }
        return true;
    }

    public int getPhredPValue(int[] pAlleleCount, int[] oAlleleCount) {
        FisherExact fe = new FisherExact(pAlleleCount[0] + pAlleleCount[1] + oAlleleCount[0] + oAlleleCount[1]);
        double p = fe.getCumlativeP(pAlleleCount[0], pAlleleCount[1], oAlleleCount[0], oAlleleCount[1]);
        if (p > 1) {
            p = 1.0;
        }
        return (int) (Math.log10(p) * (-10));
    }
    
    private double divide(int a, int b){
        if(b==0)
            return 0.0;
        else 
            return (double)a/b;
        
    }
}
