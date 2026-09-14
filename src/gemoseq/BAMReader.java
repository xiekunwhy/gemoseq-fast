package projects.gemoseq;

import java.io.File;
import java.util.Iterator;
import java.util.LinkedList;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import projects.gemoma.ExtractRNAseqEvidence.Stranded;

public class BAMReader implements Iterator<Region>{
	
	
	private Iterator<SAMRecord> recIt;
	private SamReader reader;
	private SAMRecord curr;
	private int maxIntronLength;
	private double maxcov;
	private double sample;
	
	private Region revRegion;
	private Stranded stranded;
	
	private int maxRegionLength;

	private int maxGapFilled;

	private int minQuality;

	private boolean longReads;

	private boolean collapse;
	private long absCap;
	private int pendingCap;

	private String[] restrictRefs;

	public BAMReader(int maxIntronLength, String bam, double maxcov, double sample, Stranded stranded, int minQuality, int maxRegionLength, int maxGapFilled, boolean longReads) {
		this(maxIntronLength, bam, maxcov, sample, stranded, minQuality, maxRegionLength, maxGapFilled, longReads, true, 4_000_000, 1_000_000, null);
	}

	public BAMReader(int maxIntronLength, String bam, double maxcov, double sample, Stranded stranded, int minQuality, int maxRegionLength, int maxGapFilled, boolean longReads, boolean collapse, long absCap, int pendingCap, String[] restrictRefs) {

		SamReaderFactory srf = SamReaderFactory.makeDefault();
		srf.validationStringency( ValidationStringency.SILENT );
		this.reader = srf.open(openWithIndex(bam));

		this.restrictRefs = restrictRefs;
		if(restrictRefs != null && restrictRefs.length > 0) {
			// validate against the dictionary and sort in dictionary order
			SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();
			java.util.ArrayList<String> ok = new java.util.ArrayList<String>();
			for(String ref : restrictRefs) {
				if(dict.getSequence(ref) == null) {
					throw new IllegalArgumentException("ERROR: reference "+ref+" not present in BAM header of "+bam);
				}
				ok.add(ref);
			}
			ok.sort((a,b) -> Integer.compare(dict.getSequenceIndex(a), dict.getSequenceIndex(b)));
			this.restrictRefs = ok.toArray(new String[0]);
			if(reader.hasIndex()) {
				this.recIt = new RefIterator();
			}else {
				System.out.println("WARNING: no .bai index found for "+bam+" (a .csi index is not supported by this htsjdk version); falling back to streaming the full BAM and filtering by reference. Convert with 'samtools index -b' for fast random access.");
				this.recIt = new FilterIterator();
			}
		}else {
			this.recIt = reader.iterator();
		}

		this.maxIntronLength = maxIntronLength;

		this.maxcov = maxcov;
		this.sample = sample;

		this.revRegion = null;
		this.stranded = stranded;

		this.minQuality = minQuality;

		this.maxRegionLength = maxRegionLength;

		this.maxGapFilled = maxGapFilled;

		this.longReads = longReads;

		this.collapse = collapse;
		this.absCap = absCap;
		this.pendingCap = pendingCap;

	}

	/** opens the BAM with an explicitly located index file (.bai or .csi, both naming conventions). */
	private static SamInputResource openWithIndex(String bam) {
		SamInputResource res = SamInputResource.of(new File(bam));
		String[] cands = new String[] {bam+".bai", bam+".csi", bam.replaceAll("\\.bam$", ".bai"), bam.replaceAll("\\.bam$", ".csi")};
		for(String c : cands) {
			File f = new File(c);
			if(f.isFile()) {
				if(f.lastModified() < new File(bam).lastModified()) {
					System.out.println("WARNING: index file "+c+" is older than the BAM file. A stale index silently returns incomplete results; re-create the index (e.g. 'samtools index "+bam+"') before relying on reference-restricted runs.");
				}
				res.index(f);
				break;
			}
		}
		return res;
	}

	/** iterates the records of the requested references one reference at a time (indexed query). */
	private class RefIterator implements Iterator<SAMRecord> {
		private Iterator<SAMRecord> cur = null;
		private int qi = 0;
		@Override
		public boolean hasNext() {
			while(cur == null || !cur.hasNext()) {
				if(qi >= restrictRefs.length) {
					return false;
				}
				cur = reader.query(restrictRefs[qi++], 0, 0, false);
			}
			return true;
		}
		@Override
		public SAMRecord next() {
			return cur.next();
		}
	}

	/** streams the whole BAM and yields only records on the requested references (no index needed). */
	private class FilterIterator implements Iterator<SAMRecord> {
		private final Iterator<SAMRecord> it = reader.iterator();
		private final java.util.HashSet<String> keep = new java.util.HashSet<String>(java.util.Arrays.asList(restrictRefs));
		private SAMRecord nextRec = null;
		@Override
		public boolean hasNext() {
			while(nextRec == null && it.hasNext()) {
				SAMRecord r = it.next();
				if(r.getReferenceIndex() != null && keep.contains(r.getReferenceName())) {
					nextRec = r;
				}
			}
			return nextRec != null;
		}
		@Override
		public SAMRecord next() {
			SAMRecord r = nextRec;
			nextRec = null;
			return r;
		}
	}
	
	public int getSequenceLength(String chrom) {
		SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();
		
		return dict.getSequence(chrom).getSequenceLength();
	}

	@Override
	public boolean hasNext() {
		return recIt.hasNext() || revRegion != null;
	}

	@Override
	public Region next() {
		
		if(revRegion != null) {
			Region temp = revRegion;
			revRegion = null;
			return temp;
		}
		
		Region region = new Region(maxcov,sample, stranded==Stranded.FR_UNSTRANDED ? '.' : '+', stranded, longReads, collapse, absCap, pendingCap);
		Region revTemp = new Region(maxcov,sample, stranded==Stranded.FR_UNSTRANDED ? '.' : '-', stranded, longReads, collapse, absCap, pendingCap);
		
		if(curr != null) {
			add(region,revTemp,curr);
		}
		
		while(recIt.hasNext()) {
			curr = recIt.next();
			if(curr.getReadLength() + maxIntronLength < curr.getAlignmentEnd()-curr.getAlignmentStart()) {
				curr = null;
				continue;
			}
			if(curr.getMappingQuality() < minQuality) {
				curr = null;
				continue;
			}
			if( 
					(region.getReferenceIndex() != null && !curr.getReferenceIndex().equals(region.getReferenceIndex() ) ) ||
					(revTemp.getReferenceIndex() != null && !curr.getReferenceIndex().equals(revTemp.getReferenceIndex() ) ) 
				) {
				return join(region,revTemp);
			}
			int maxEnd = Math.max(region.getRegionEnd() == null ? -1 : region.getRegionEnd(), revTemp.getRegionEnd() == null ? -1 : revTemp.getRegionEnd());
			if(
					maxEnd > -1 && maxEnd < curr.getAlignmentStart()
				) {
				
				boolean filled = false;
				if(maxEnd - Math.min(region.getRegionStart() == null ? maxEnd : region.getRegionStart(), revTemp.getRegionStart() == null ? maxEnd : revTemp.getRegionStart()) < maxRegionLength) {
					
					if(isFwd(curr) && region.getRegionEnd() != null && region.getRegionEnd()+maxGapFilled>= curr.getAlignmentStart()) {
						SAMRecord[] dummies = getDummyPair(curr, region.getRegionEnd()-1, curr.getAlignmentStart()+1);
						add(region,revTemp,dummies[0]);
						add(region,revTemp,dummies[1]);
						filled = true;
					}
					if(!isFwd(curr) && revTemp.getRegionEnd() != null && revTemp.getRegionEnd()+maxGapFilled>= curr.getAlignmentStart()) {
						SAMRecord[] dummies = getDummyPair(curr, revTemp.getRegionEnd()-1, curr.getAlignmentStart()+1);
						add(region,revTemp,dummies[0]);
						add(region,revTemp,dummies[1]);
						filled = true;
					}
					
				}
				if(!filled) {
					return join(region,revTemp);
				}
			}
			add(region,revTemp,curr);
		}
		return join(region,revTemp);
	}
	
	
	private static final SAMRecord getDummy(SAMRecord curr, int start, int end ) {
		SAMRecord dummy = new SAMRecord(curr.getHeader());
		dummy.setAlignmentStart(start);
		dummy.setMappingQuality(255);
		dummy.setReferenceIndex(curr.getReferenceIndex());
		dummy.setReferenceName(curr.getReferenceName());
		dummy.setReadName("dummy");
		LinkedList<CigarElement> cili = new LinkedList<CigarElement>();
		cili.add(new CigarElement(end-start+1, CigarOperator.M));
		dummy.setCigar(new Cigar(cili));
		dummy.setReadPairedFlag(true);
		dummy.setProperPairFlag(true);
		dummy.setMateReferenceIndex(dummy.getReferenceIndex());
		dummy.setMateReferenceName(dummy.getReferenceName());
		return dummy;
	}

	private SAMRecord[] getDummyPair(SAMRecord curr, int start, int end ) {
		SAMRecord dummy1 = getDummy(curr, start, end);
		SAMRecord dummy2 = getDummy(curr, start, end);
		dummy1.setFirstOfPairFlag(true);
		dummy2.setSecondOfPairFlag(true);
		if(stranded == Stranded.FR_UNSTRANDED) {
			dummy1.setReadNegativeStrandFlag(false);
			dummy2.setReadNegativeStrandFlag(true);
		}else {
			if(curr.getReadPairedFlag()) {
				if(curr.getFirstOfPairFlag()) {
					dummy1.setReadNegativeStrandFlag(curr.getReadNegativeStrandFlag());
					dummy1.setMateNegativeStrandFlag(curr.getMateNegativeStrandFlag());
					dummy2.setReadNegativeStrandFlag(curr.getMateNegativeStrandFlag());
					dummy2.setMateNegativeStrandFlag(curr.getReadNegativeStrandFlag());
				}else {
					dummy1.setReadNegativeStrandFlag(curr.getMateNegativeStrandFlag());
					dummy1.setMateNegativeStrandFlag(curr.getReadNegativeStrandFlag());
					dummy2.setReadNegativeStrandFlag(curr.getReadNegativeStrandFlag());
					dummy2.setMateNegativeStrandFlag(curr.getMateNegativeStrandFlag());
				}
			}
		}
		return new SAMRecord[] {dummy1,dummy2};
	}
	
	
	private Region join(Region region, Region revTemp) {
		
		if(region.getReferenceIndex() == null) {
			return revTemp;
		}
		if(revTemp.getReferenceIndex() == null) {
			return region;
		}
		
				
		if(!longReads) {
			if(revTemp.getTheoreticalNumberOfReads() < 10) {
				return region;
			}
			if(region.getTheoreticalNumberOfReads() < 10) {
				return revTemp;
			}
		}
		if(revTemp.getTheoreticalNumberOfReads()*50 < region.getTheoreticalNumberOfReads()) {//TODO 50
			double jacc = jaccard(region, revTemp);
			if(jacc > 0.5) {
				region.join(revTemp);//TODO discard instead of join?
				return region;
			}
		}

		if(region.getTheoreticalNumberOfReads()*50 < revTemp.getTheoreticalNumberOfReads()) {
			double jacc = jaccard(region,revTemp);
			if(jacc > 0.5) {
				revTemp.join(region);
				return revTemp;
			}
		}
		
		revRegion = revTemp;
		
		region.setRevRegion(revRegion);
		revRegion.setRevRegion(region);
		
		return region;
	}
	
	private double cor(Region region, Region revTemp) {
		
		int[] covFwd = region.getCoverageByBlocks();
		int[] covRev = revTemp.getCoverageByBlocks();
		
		int regionStart = region.getRegionStart();
		int regionStartRev = revTemp.getRegionStart();
		
		int globalStart = Math.min(regionStart, regionStartRev);
		int globalEnd = Math.max(regionStart+covFwd.length, regionStartRev+covRev.length);
		
		
		double ex = 0.0;
		double ey = 0.0;
		double exy = 0.0;
		double ex2 = 0.0;
		double ey2 = 0.0;
		double n = globalEnd-globalStart;
		for(int i=globalStart;i<globalEnd;i++) {
			double v1 = i<regionStart || i>= regionStart+covFwd.length ? 0.0 : covFwd[i-regionStart];
			double v2 = i<regionStartRev || i >= regionStartRev+covRev.length ? 0.0: covRev[i-regionStartRev];
			
			ex += v1;
			ey += v2;
			exy += v1*v2;
			ex2 += v1*v1;
			ey2 += v2*v2;
			
		}
		
		ex /= n;
		ey /= n;
		exy /= n;
		ex2 /= n;
		ey2 /= n;
		
		double cov = exy - ex*ey;
		double v1 = ex2 - ex*ex;
		double v2 = ey2 - ey*ey;
		
		if(v1 == 0 || v2 == 0) {
			return 0;
		}else {
			return cov/Math.sqrt(v1*v2);	
		}
	}
	
	
	private static double getThreshold(int[] vals) {
		double sum = 0.0;
		double n = 0.0;
		for(int i=0;i<vals.length;i++) {
			if(vals[i] > 0) {
				sum += vals[i];
				n ++;
			}
		}
		return sum/n*0.1;
	}
	
	private double jaccard(Region region, Region revTemp) {
		
		int[] covFwd = region.getCoverageByBlocks();
		int[] covRev = revTemp.getCoverageByBlocks();
		
		double t1 = getThreshold(covFwd);
		double t2 = getThreshold(covRev);
		
		int regionStart = region.getRegionStart();
		int regionStartRev = revTemp.getRegionStart();
		
		
		int globalStart = Math.min(regionStart, regionStartRev);
		int globalEnd = Math.max(regionStart+covFwd.length, regionStartRev+covRev.length);
		
		
		double inter = 0.0; 
		double union1 = 0.0;
		double union2 = 0.0;
		double sum1 = 0.0;
		double sum2 = 0.0;
		for(int i=globalStart;i<globalEnd;i++) {
			double v1 = i<regionStart || i>= regionStart+covFwd.length ? 0.0 : covFwd[i-regionStart];
			double v2 = i<regionStartRev || i >= regionStartRev+covRev.length ? 0.0: covRev[i-regionStartRev];
			
			sum1 += v1;
			sum2 += v2;
			
			if(v1 > t1 && v2 > t2) {
				inter += 1;
			}
			if(v1 > t1) {
				union1 ++;
			}
			if(v2 > t2) {
				union2 ++;
			}
		}
		double union = sum1 > sum2 ? union2 : union1;
		
		return inter/union;
		
	}
	

	
	private boolean isFwd(SAMRecord sr) {
		if(stranded == Stranded.FR_UNSTRANDED) {
			return true;
		}else if(stranded == Stranded.FR_FIRST_STRAND){
			if( 
					(sr.getReadPairedFlag() && sr.getFirstOfPairFlag() && sr.getReadNegativeStrandFlag() ) ||
					(sr.getReadPairedFlag() && sr.getSecondOfPairFlag() && sr.getMateNegativeStrandFlag() ) ||
					( (!sr.getReadPairedFlag() || !sr.getProperPairFlag()) && sr.getReadNegativeStrandFlag() )
				) {
				return true;
			}else {
				return false;
			}
		}else {
			if( 
					(sr.getReadPairedFlag() && sr.getFirstOfPairFlag() && sr.getReadNegativeStrandFlag() ) ||
					(sr.getReadPairedFlag() && sr.getSecondOfPairFlag() && sr.getMateNegativeStrandFlag() ) ||
					( (!sr.getReadPairedFlag() || !sr.getProperPairFlag()) && sr.getReadNegativeStrandFlag() )
				) {
				return false;
			}else {
				return true;
			}
		}
	}
	

	private void add(Region region, Region revTemp, SAMRecord sr) {
		if(isFwd(sr)) {
			region.addRead(sr);
		}else {
			revTemp.addRead(sr);
		}
	}
	
	
	
}
