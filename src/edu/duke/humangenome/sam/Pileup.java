/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.sam;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.SamLocusIterator.LocusInfo;
import net.sf.picard.util.SamLocusIterator.RecordAndOffset;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.SequenceUtil;

/**
 *
 * @author Yongzhuang Liu
 */
public class Pileup {

    private LocusInfo locusInfo;
    private int[] baseCountTable = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    private char[] bases;
    private int[] baseQualities;
    private int[] mappingQualities;
    
    public Pileup(LocusInfo locusInfo) {
        this.locusInfo = locusInfo;
        initialize();
    }

    public int getDepth() {
        if (locusInfo == null) {
            return 0;
        } else {
            return locusInfo.getRecordAndPositions().size();
        }
    }
    
    private void initialize() {
        if (locusInfo == null) {
            return;
        } else {
            int depth=getDepth();
            bases=new char[depth];
            baseQualities=new int[depth];
            mappingQualities=new int[depth];
            for (int i = 0; i < depth; i++) {
                byte base = locusInfo.getRecordAndPositions().get(i).getReadBase();
                bases[i]=(char)base;
                baseQualities[i] = (int)locusInfo.getRecordAndPositions().get(i).getBaseQuality();
                mappingQualities[i]= locusInfo.getRecordAndPositions().get(i).getRecord().getMappingQuality();
                boolean reverse = locusInfo.getRecordAndPositions().get(i).getRecord().getReadNegativeStrandFlag();
                if (reverse) {
                    base = (byte) (base + 32);
                }
                switch ((char) base) {
                    case 'A':
                        baseCountTable[0]++;
                        break;
                    case 'C':
                        baseCountTable[1]++;
                        break;
                    case 'G':
                        baseCountTable[2]++;
                        break;
                    case 'T':
                        baseCountTable[3]++;
                        break;
                    case 'N':
                        baseCountTable[4]++;
                        break;
                    case 'a':
                        baseCountTable[5]++;
                        break;
                    case 'c':
                        baseCountTable[6]++;
                        break;
                    case 'g':
                        baseCountTable[7]++;
                        break;
                    case 't':
                        baseCountTable[8]++;
                        break;
                    case 'n':
                        baseCountTable[9]++;
                        break;
                    default:
                        break;
                }
            }
            return;
        }
    }

    public int getForwardBaseCount(char base) {
        base = Character.toUpperCase(base);
        switch (base) {
            case 'A':
                return baseCountTable[0];
            case 'C':
                return baseCountTable[1];
            case 'G':
                return baseCountTable[2];
            case 'T':
                return baseCountTable[3];
            case 'N':
                return baseCountTable[4];
            default:
                return 0;
        }
    }

    public int getReverseBaseCount(char base) {
        base = Character.toLowerCase(base);
        switch (base) {
            case 'a':
                return baseCountTable[5];
            case 'c':
                return baseCountTable[6];
            case 'g':
                return baseCountTable[7];
            case 't':
                return baseCountTable[8];
            case 'n':
                return baseCountTable[9];
            default:
                return 0;
        }
    }

    public int getBothStrandBaseCount(char base) {
        base = Character.toUpperCase(base);
        switch (base) {
            case 'A':
                return baseCountTable[0] + baseCountTable[5];
            case 'C':
                return baseCountTable[1] + baseCountTable[6];
            case 'G':
                return baseCountTable[2] + baseCountTable[7];
            case 'T':
                return baseCountTable[3] + baseCountTable[8];
            case 'N':
                return baseCountTable[4] + baseCountTable[9];
            default:
                return 0;
        }
    }

    public LocusInfo getLocusInfo() {
        return locusInfo;
    }
    
    public int getMeanBaseQuality(char base) {
        base = Character.toUpperCase(base);
        int quality = 0;
        int count = 0;
        for (int i = 0; i < getDepth(); i++) {
            if (bases[i] == base) {
                quality += baseQualities[i];
                count++;
            }
        }
        if (count == 0) {
            return 0;
        } else {
            return quality / count;
        }
    }

    public int getMeanMappingQuality(char base) {
        base = Character.toUpperCase(base);
        int quality = 0;
        int count = 0;
        for (int i = 0; i < getDepth(); i++) {
            if (bases[i] == base && mappingQualities[i]!=0) {
                quality += mappingQualities[i];
                count++;
            }
        }
        if (count == 0) {
            return 0;
        } else {
            return quality / count;
        }
    }

    public double getFractionOfMQ0Reads(char base) {

        base = Character.toUpperCase(base);
        int count = 0;
        int numOfMQ0Reads = 0;
        for (int i = 0; i < getDepth(); i++) {
            if (bases[i] == base) {
                count++;
                if (mappingQualities[i] == 0) {
                    numOfMQ0Reads++;
                }
            }
        }
        if (count == 0) {
            return 0;
        } else {
            return (double) numOfMQ0Reads / count;
        }
    }
     
    public int getMeanDistanceToThreePrime(char base) {
        int count = 0;
        int distance = 0;
        base = Character.toUpperCase(base);
        for (int i = 0; i < getDepth(); i++) {
            if (bases[i] == base) {
                count++;
                boolean reverse = locusInfo.getRecordAndPositions().get(i).getRecord().getReadNegativeStrandFlag();
                if (!reverse) {
                    distance += locusInfo.getRecordAndPositions().get(i).getRecord().getReadLength() - locusInfo.getRecordAndPositions().get(i).getOffset() - 1;
                } else {
                    distance += locusInfo.getRecordAndPositions().get(i).getOffset();
                }
            }
        }
        if (count == 0) {
            return 0;
        } else {
            return distance / count;
        }
    }
    
    public double getFractionOfSoftClippedReads(char base) {
        int count = 0;
        int numOfSoftClippedReads = 0;
        base = Character.toUpperCase(base);
        for (int i = 0; i < getDepth(); i++) {
            if (bases[i] == base) {
                count++;
                String cigar = locusInfo.getRecordAndPositions().get(i).getRecord().getCigarString();
                if (cigar.contains("S")) {
                    numOfSoftClippedReads++;
                }
            }
        }
        if (count == 0) {
            return 0;
        } else {
            return (double) numOfSoftClippedReads / count;
        }
    }

    public int getNearbyIndels(char base) {
        base = Character.toUpperCase(base);
        int indelCount = 0;
        for (int i = 0; i < getDepth(); i++) {
            if (bases[i] == base) {
                SAMRecord record = locusInfo.getRecordAndPositions().get(i).getRecord();
                Cigar cigar = record.getCigar();
                for (int j = 0; j < cigar.numCigarElements(); j++) {
                    CigarOperator operator = cigar.getCigarElement(j).getOperator();
                    int length = cigar.getCigarElement(j).getLength();
                    switch (operator) {
                        case M:
                            break;
                        case I:
                            indelCount++;
                            break;
                        case D:
                            indelCount++;
                            break;
                        default:
                            break;
                    }
                }
            }
        }
        return indelCount;
    }
        
    public int getOverallMismatches(ReferenceSequence referenceSequence, char base) {
        base = Character.toUpperCase(base);
        byte[] referenceBases = referenceSequence.getBases();
        int count = 0;
        for (int i = 0; i < getDepth(); i++) {
            if (bases[i] == base) {
                SAMRecord record = locusInfo.getRecordAndPositions().get(i).getRecord();
                count += SequenceUtil.countMismatches(record, referenceBases);
            }
        }
        return count;
    }
    
}
