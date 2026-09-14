package projects.gemoma;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class ReadStats {

	private double meanSplit;
	private double sdSplit;
	private double meanReadLen;
	private double factor;

	public ReadStats(int minIntronLength, double factor, String... bams) {
		this(minIntronLength, factor, null, bams);
	}

	/**
	 * Computes read statistics, optionally restricted to a set of references.
	 * With restrictRefs != null, records are fetched by indexed queries per reference
	 * (fast when a valid .bai/.csi index exists) or by streaming the full BAM and
	 * filtering (fallback when no index is present).
	 */
	public ReadStats(int minIntronLength, double factor, String[] restrictRefs, String... bams) {
		this.factor = factor;
		SamReaderFactory srf = SamReaderFactory.makeDefault();
		srf.validationStringency( ValidationStringency.SILENT );

		Accumulator acc = new Accumulator(minIntronLength);

		for(int i=0;i<bams.length;i++) {
			SamReader reader = srf.open(openWithIndex(bams[i]));

			if(restrictRefs == null || restrictRefs.length == 0) {
				Iterator<SAMRecord> recIt = reader.iterator();
				while(recIt.hasNext()) {
					acc.consume(recIt.next());
				}
			}else {
				SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();
				ArrayList<String> refs = new ArrayList<String>();
				for(String ref : restrictRefs) {
					if(dict.getSequence(ref) == null) {
						throw new IllegalArgumentException("ERROR: reference "+ref+" not present in BAM header of "+bams[i]);
					}
					refs.add(ref);
				}
				final ArrayList<String> sorted = new ArrayList<String>(refs);
				sorted.sort((a,b) -> Integer.compare(dict.getSequenceIndex(a), dict.getSequenceIndex(b)));
				if(reader.hasIndex()) {
					for(String ref : sorted) {
						Iterator<SAMRecord> recIt = reader.query(ref, 0, 0, false);
						while(recIt.hasNext()) {
							acc.consume(recIt.next());
						}
					}
				}else {
					System.out.println("WARNING: no index found for "+bams[i]+"; computing read statistics by streaming the full BAM and filtering by reference.");
					HashSet<String> keep = new HashSet<String>(sorted);
					Iterator<SAMRecord> recIt = reader.iterator();
					while(recIt.hasNext()) {
						SAMRecord rec = recIt.next();
						if(rec.getReferenceIndex() != null && keep.contains(rec.getReferenceName())) {
							acc.consume(rec);
						}
					}
				}
			}
			try {
				reader.close();
			} catch(Exception e) {}
		}

		meanSplit = acc.sum/acc.n;
		sdSplit = Math.sqrt( acc.sumSq/acc.n - meanSplit*meanSplit );
		meanReadLen = acc.lenSum/acc.lenN;
		//System.out.println("stats: "+meanSplit+" "+sdSplit);

	}

	/** opens the BAM with an explicitly located index file (.bai or .csi, both naming conventions). */
	private static SamInputResource openWithIndex(String bam) {
		SamInputResource res = SamInputResource.of(new File(bam));
		String[] cands = new String[] {bam+".bai", bam+".csi", bam.replaceAll("\\.bam$", ".bai"), bam.replaceAll("\\.bam$", ".csi")};
		for(String c : cands) {
			File f = new File(c);
			if(f.isFile()) {
				if(f.lastModified() < new File(bam).lastModified()) {
					System.out.println("WARNING: index file "+c+" is older than the BAM file. A stale index silently returns incomplete results; re-create the index before relying on reference-restricted runs.");
				}
				res.index(f);
				break;
			}
		}
		return res;
	}

	private static class Accumulator {
		final int minIntronLength;
		double sum = 0;
		double sumSq = 0;
		double n = 0;
		double lenSum = 0;
		double lenN = 0;

		Accumulator(int minIntronLength) {
			this.minIntronLength = minIntronLength;
		}

		void consume(SAMRecord curr) {
			lenSum += curr.getReadLength();
			lenN++;

			List<AlignmentBlock> li = curr.getAlignmentBlocks();
			if(li.size() > 1) {
				Iterator<AlignmentBlock> it = li.iterator();
				AlignmentBlock al = it.next();
				while(it.hasNext()) {
					AlignmentBlock ne = it.next();
					double diff = Math.abs( ne.getReferenceStart()-al.getReferenceStart()-al.getLength())+1;
					if(diff>minIntronLength) {
						sum += diff;
						sumSq += diff*diff;
						n++;
					}
					al = ne;
				}
			}
		}
	}

	public double getSdSplitLength() {
		return sdSplit;
	}

	public double getMeanSplitLength() {
		return meanSplit;
	}

	public double getMeanReadLength() {
		return meanReadLen;
	}

	public boolean isOK(int len, int num) {
		double min = factor*(len - getMeanSplitLength())/getSdSplitLength();

		return num >= min;
	}

}
