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
		this(minIntronLength, factor, null, false, bams);
	}

	/**
	 * Computes read statistics, optionally restricted to a set of references.
	 * With restrictRefs != null, records are fetched by indexed queries per reference
	 * (fast when a valid .bai/.csi index exists) or by streaming the full BAM and
	 * filtering (fallback when no index is present, or when forceStream is true).
	 */
	public ReadStats(int minIntronLength, double factor, String[] restrictRefs, String... bams) {
		this(minIntronLength, factor, restrictRefs, false, bams);
	}

	public ReadStats(int minIntronLength, double factor, String[] restrictRefs, boolean forceStream, String... bams) {
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
				if(reader.hasIndex() && !forceStream) {
					for(String ref : sorted) {
						java.io.Closeable recIt = (java.io.Closeable) reader.query(ref, 0, 0, false);
						Iterator<SAMRecord> it = (Iterator<SAMRecord>) recIt;
						while(it.hasNext()) {
							acc.consume(it.next());
						}
						// htsjdk allows only one open iterator per reader; close before the next query
						try {
							recIt.close();
						} catch(Exception e) {}
					}
				}else {
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

	private ReadStats(double meanSplit, double sdSplit, double meanReadLen, double factor) {
		this.meanSplit = meanSplit;
		this.sdSplit = sdSplit;
		this.meanReadLen = meanReadLen;
		this.factor = factor;
	}

	/** writes the statistics to a simple text file (see fromFile for the format). */
	public void toFile(String path) throws java.io.IOException {
		toFile(path, null);
	}

	public void toFile(String path, String source) throws java.io.IOException {
		java.io.PrintWriter w = new java.io.PrintWriter(new java.io.FileWriter(path));
		w.println("#ReadStats v1" + (source != null ? "\t" + source : ""));
		w.println("meanSplit\t" + meanSplit);
		w.println("sdSplit\t" + sdSplit);
		w.println("meanReadLen\t" + meanReadLen);
		w.println("factor\t" + factor);
		w.close();
	}

	/** reads statistics previously written with toFile; throws on malformed content. */
	public static ReadStats fromFile(String path) throws java.io.IOException {
		java.io.BufferedReader rd = new java.io.BufferedReader(new java.io.FileReader(path));
		Double meanSplit = null, sdSplit = null, meanReadLen = null;
		double factor = 1.0;
		String line;
		while( (line = rd.readLine()) != null ) {
			line = line.trim();
			if(line.isEmpty() || line.startsWith("#")) {
				continue;
			}
			String[] p = line.split("\\s+", 2);
			if(p.length != 2) {
				continue;
			}
			String k = p[0], v = p[1].trim();
			try {
				if(k.equals("meanSplit")) {
					meanSplit = Double.parseDouble(v);
				}else if(k.equals("sdSplit")) {
					sdSplit = Double.parseDouble(v);
				}else if(k.equals("meanReadLen")) {
					meanReadLen = Double.parseDouble(v);
				}else if(k.equals("factor")) {
					factor = Double.parseDouble(v);
				}
			} catch(NumberFormatException e) {
				rd.close();
				throw new java.io.IOException("malformed value in ReadStats file "+path+": "+line);
			}
		}
		rd.close();
		if(meanSplit == null || sdSplit == null) {
			throw new java.io.IOException("ReadStats file "+path+" lacks meanSplit/sdSplit entries");
		}
		return new ReadStats(meanSplit, sdSplit, meanReadLen == null ? 0.0 : meanReadLen, factor);
	}

	public boolean isOK(int len, int num) {
		double min = factor*(len - getMeanSplitLength())/getSdSplitLength();

		return num >= min;
	}

}
