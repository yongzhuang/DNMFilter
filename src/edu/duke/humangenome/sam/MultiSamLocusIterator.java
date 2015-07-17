/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.sam;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.util.SamLocusIterator.LocusInfo;
import net.sf.samtools.util.CloseableIterator;

/**
 *
 * @author Yongzhuang Liu
 */
public class MultiSamLocusIterator {

    private PriorityQueue<ComparableSamLocusIterator> pq;
    private int numberOfIterators;
    private LocusInfoComparator comparator;

    public MultiSamLocusIterator() {
        this.pq = new PriorityQueue<ComparableSamLocusIterator>();
        this.comparator = getComparator();
        this.numberOfIterators = 0;
    }

    public void add(SamLocusIterator samLocusIterator, int index, int mappingQualityScoreCutoff, int qualityScoreCutoff) {
        samLocusIterator.setEmitUncoveredLoci(false);
        samLocusIterator.setMappingQualityScoreCutoff(mappingQualityScoreCutoff);
        samLocusIterator.setQualityScoreCutoff(qualityScoreCutoff);
        samLocusIterator.iterator();
        pq.add(new ComparableSamLocusIterator(samLocusIterator, comparator, index));
        numberOfIterators++;
    }

    public void add(SamLocusIterator samLocusIterator, int index) {
        samLocusIterator.setEmitUncoveredLoci(false);
        samLocusIterator.iterator();
        pq.add(new ComparableSamLocusIterator(samLocusIterator, comparator, index));
        numberOfIterators++;
    }

    public int getNumberOfIterators() {
        return this.numberOfIterators;
    }

    public Map<Integer, LocusInfo> getLocusInfos() {
        Map<Integer, LocusInfo> map = new HashMap();
        List<ComparableSamLocusIterator> list = new ArrayList();
        if (pq.isEmpty()) {
            return null;
        }
        for (int i = 0; i < numberOfIterators; i++) {
            ComparableSamLocusIterator end = pq.poll();
            list.add(end);
            if (!pq.isEmpty()) {
                ComparableSamLocusIterator next = pq.peek();
                if (end.peek().getPosition() != next.peek().getPosition()) {
                    break;
                }
            } else {
                break;
            }
        }
        for (int i = 0; i < list.size(); i++) {
            map.put(list.get(i).getIndex(), list.get(i).next());
            if (list.get(i).hasNext()) {
                pq.add(list.get(i));
            }
        }
        if (map.size() > 0) {
            return map;
        } else {
            return null;
        }
    }

    private LocusInfoComparator getComparator() {
        return new LocusInfoComparator() {
            public int compare(LocusInfo locusInfo1, LocusInfo locusInfo2) {
                int result = 0;
                if (locusInfo1.getSequenceIndex() < locusInfo2.getSequenceIndex()) {
                    return -1;
                }
                if (locusInfo1.getSequenceIndex() > locusInfo2.getSequenceIndex()) {
                    return 1;
                } else {
                    if (locusInfo1.getPosition() < locusInfo2.getPosition()) {
                        return -1;
                    }
                    if (locusInfo1.getPosition() > locusInfo2.getPosition()) {
                        return 1;
                    }
                }
                return result;
            }
        };
    }

    public void close() {
        for (CloseableIterator<LocusInfo> iterator : pq) {
            iterator.close();
        }
    }
}
