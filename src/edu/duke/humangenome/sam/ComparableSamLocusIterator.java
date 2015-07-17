/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.sam;

import java.util.Comparator;
import java.util.Iterator;
import net.sf.picard.util.PeekableIterator;
import net.sf.picard.util.SamLocusIterator.LocusInfo;
import net.sf.samtools.util.CloseableIterator;

/**
 *
 * @author Yongzhuang Liu
 */
public class ComparableSamLocusIterator extends PeekableIterator<LocusInfo> implements Comparable<ComparableSamLocusIterator> {

    private final Comparator<LocusInfo> comparator;
    private int index;

    public ComparableSamLocusIterator(CloseableIterator<LocusInfo> iterator, Comparator<LocusInfo> comparator, int index) {
        super(iterator);
        this.comparator = comparator;
        this.index = index;
    }

    public int compareTo(final ComparableSamLocusIterator that) {
        LocusInfo locusInfo1 = this.peek();
        LocusInfo locusInfo2 = that.peek();
        return comparator.compare(locusInfo1, locusInfo2);
    }

    public int getIndex() {
        return index;
    }
}

