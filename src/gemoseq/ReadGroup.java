package projects.gemoseq;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/**
 * A group of fragments (read pairs or single reads) that share the exact same
 * alignment structure. All depth-dependent counts are collapsed into {@link #count},
 * so memory and downstream runtime become (nearly) independent of sequencing depth.
 *
 * For a normal group, {@link #frags} has length 1 and frags[0] holds the one or two
 * mates of the fragment structure; {@link #fragW}{@code [0]} is the number of
 * collapsed fragments.
 *
 * For the artificial "dummy" reads used to fill short gaps, all dummy records of a
 * region are accumulated in a single shared group (this reproduces the original
 * behavior where all reads named "dummy" shared one index in the region id map):
 * {@link #frags} then has one entry per dummy record pair and {@link #fragW} is 1
 * per entry.
 */
public class ReadGroup {

	/** Compact representation of one mate (one aligned record). */
	public static final class Mate {
		/** 1-based reference starts of the alignment blocks. */
		public final int[] bRef;
		/** lengths of the alignment blocks. */
		public final int[] bLen;
		/**
		 * per block (index 0 unused): bit0 = the cigar element right before this
		 * block is a skipped (non-indel) region, i.e. the short-read isSplit test;
		 * bit1 = a skipped region has been seen while matching this block (the
		 * long-read isSplit test).
		 */
		public final byte[] bFlag;
		/** long reads: length of a deletion right before the skipped region. */
		public final int[] delBefore;
		/** long reads: length of a deletion right before this block. */
		public final int[] delOff;
		/** number of mismatches in the first/last <=10 bp of each block. */
		public final int nmm;
		public final int alignStart;
		public final int alignEnd;

		public Mate(int[] bRef, int[] bLen, byte[] bFlag, int[] delBefore, int[] delOff, int nmm, int alignStart, int alignEnd) {
			this.bRef = bRef;
			this.bLen = bLen;
			this.bFlag = bFlag;
			this.delBefore = delBefore;
			this.delOff = delOff;
			this.nmm = nmm;
			this.alignStart = alignStart;
			this.alignEnd = alignEnd;
		}
	}

	public static final int FLAG_SPLIT_SHORT = 1;
	public static final int FLAG_SKIPPED_SEEN = 2;

	public final ArrayList<Mate[]> frags = new ArrayList<Mate[]>(1);

	/** weight (number of collapsed fragments) per entry of {@link #frags}; usually a single entry. */
	public int[] fragW = new int[1];

	/** number of collapsed fragments in total (sum of fragW). */
	public int count;

	/** region-local group id (row index in the quantification matrix). */
	public int gid = -1;

	public char strand = '.';

	public boolean dummy = false;

	public int start = Integer.MAX_VALUE;
	public int end = -1;

	public ReadGroup() {
	}

	/** adds a fragment instance with the given mates and weight. */
	public void add(Mate[] mates, int w) {
		int n = frags.size();
		frags.add(mates);
		if(n >= fragW.length) {
			int[] nw = new int[Math.max(2*n, 4)];
			System.arraycopy(fragW, 0, nw, 0, n);
			fragW = nw;
		}
		fragW[n] = w;
		count += w;
		for(Mate m : mates) {
			if(m.alignStart < start) {
				start = m.alignStart;
			}
			if(m.alignEnd > end) {
				end = m.alignEnd;
			}
		}
	}

	public int numFrags() {
		return frags.size();
	}

	/** total number of records (mates) times fragment weight, used for coverage bookkeeping. */
	public int numRecords() {
		int n = 0;
		for(int i=0;i<frags.size();i++) {
			n += fragW[i]*frags.get(i).length;
		}
		return n;
	}

	/**
	 * Converts one SAMRecord into the compact mate representation. The cigar/block
	 * matching logic mirrors the original ReadGraph.addRead / SplicingGraph.addRead
	 * code paths exactly, so that downstream decisions are unchanged.
	 */
	public static Mate convert(SAMRecord read, boolean longReads) {
		List<AlignmentBlock> blockLi = read.getAlignmentBlocks();

		int nb = blockLi.size();
		int[] bRef = new int[nb];
		int[] bLen = new int[nb];
		byte[] bFlag = new byte[nb];
		int[] delBefore = new int[nb];
		int[] delOff = new int[nb];

		int nmm = 0;

		if(nb > 1 && !longReads) {
			// mismatch check of the original ReadGraph.addRead (short reads only)
			char[] chr = Genome.genome.getChromosome(read.getReferenceName());
			byte[] rs = read.getReadBases();
			if(rs != null && rs.length > 0) {
				for(AlignmentBlock block : blockLi) {
					int l = block.getLength();
					int c = block.getReferenceStart()-1;
					int r = block.getReadStart()-1;
					for(int i=0;i<l/2 && i<10;i++) {
						if(chr[c+i] != (char)rs[r+i]) {
							nmm++;
						}
					}
					for(int i=l-1;i>=l/2 && i>l-11;i--) {
						if(chr[c+i] != (char)rs[r+i]) {
							nmm++;
						}
					}
				}
			}
		}

		Cigar cigar = read.getCigar();
		List<CigarElement> els = cigar.getCigarElements();
		int ci = 0;

		for(int k=0;k<nb;k++) {
			AlignmentBlock block = blockLi.get(k);
			bRef[k] = block.getReferenceStart();
			bLen[k] = block.getLength();

			CigarElement el = els.get(ci++);
			CigarElement prev = null;
			boolean hasBeenSkipped = false;
			int db = 0;
			while(el.getLength() != block.getLength() || !el.getOperator().isAlignment()) {
				if(el.getOperator().isIndelOrSkippedRegion() && !el.getOperator().isIndel()) {
					hasBeenSkipped = true;
					db = longReads && prev != null && prev.getOperator()==CigarOperator.DELETION ? prev.getLength() : 0;
				}
				prev = el;
				el = els.get(ci++);
			}
			int dof = longReads && prev != null && prev.getOperator()==CigarOperator.DELETION ? prev.getLength() : 0;

			if(prev != null && prev.getOperator().isIndelOrSkippedRegion() && !prev.getOperator().isIndel()) {
				bFlag[k] |= FLAG_SPLIT_SHORT;
			}
			if(hasBeenSkipped) {
				bFlag[k] |= FLAG_SKIPPED_SEEN;
			}
			delBefore[k] = db;
			delOff[k] = dof;
		}

		return new Mate(bRef, bLen, bFlag, delBefore, delOff, nmm, read.getAlignmentStart(), read.getAlignmentEnd());
	}

	/** Structural signature used for collapsing; everything that can influence downstream decisions. */
	public static final class Signature {
		private final Mate[][] mates;
		private final char strand;
		private final int hash;

		public Signature(Mate[] frag, char strand) {
			this.mates = new Mate[][] {frag};
			this.strand = strand;
			this.hash = computeHash();
		}

		private Signature(Mate[][] mates, char strand) {
			this.mates = mates;
			this.strand = strand;
			this.hash = computeHash();
		}

		private int computeHash() {
			int h = strand;
			for(Mate[] frag : mates) {
				for(Mate m : frag) {
					h = 31*h + Arrays.hashCode(m.bRef);
					h = 31*h + Arrays.hashCode(m.bLen);
					h = 31*h + Arrays.hashCode(m.bFlag);
					h = 31*h + Arrays.hashCode(m.delBefore);
					h = 31*h + Arrays.hashCode(m.delOff);
					h = 31*h + m.nmm;
				}
			}
			return h;
		}

		private static boolean mateEq(Mate a, Mate b) {
			return a.nmm == b.nmm
					&& Arrays.equals(a.bRef, b.bRef)
					&& Arrays.equals(a.bLen, b.bLen)
					&& Arrays.equals(a.bFlag, b.bFlag)
					&& Arrays.equals(a.delBefore, b.delBefore)
					&& Arrays.equals(a.delOff, b.delOff);
		}

		@Override
		public boolean equals(Object o) {
			Signature s = (Signature) o;
			if(s.strand != strand || s.mates.length != mates.length) {
				return false;
			}
			for(int i=0;i<mates.length;i++) {
				if(s.mates[i].length != mates[i].length) {
					return false;
				}
				for(int j=0;j<mates[i].length;j++) {
					if(!mateEq(s.mates[i][j], mates[i][j])) {
						return false;
					}
				}
			}
			return true;
		}

		@Override
		public int hashCode() {
			return hash;
		}
	}
}
