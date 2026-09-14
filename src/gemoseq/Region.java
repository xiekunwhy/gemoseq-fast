package projects.gemoseq;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Random;

import de.jstacs.utils.IntList;
import de.jstacs.utils.ToolBox;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import projects.gemoma.ExtractRNAseqEvidence.Stranded;
import projects.gemoseq.ReadGroup.Mate;

/**
 * A region of the genome holding the mapped reads in collapsed form:
 * fragments (read pairs or single reads) with identical alignment structure are
 * stored once with a weight ({@link ReadGroup}), which makes memory consumption
 * and downstream runtime nearly independent of sequencing depth.
 *
 * Down-sampling semantics of the original implementation are kept: reads are
 * accepted with probability sampleProb (halved every time the region grows beyond
 * its coverage limit), and the limit is the minimum of maxcov*regionLength (as
 * before) and an absolute cap maxReadsPerRegion (new, guards against extreme
 * coverage spikes and gap-free mega regions).
 */
public class Region {

	private static Random r = new Random(113);
	private double maxcov;
	private double sample;

	private LinkedHashMap<ReadGroup.Signature, ReadGroup> groups;
	private ArrayList<ReadGroup> rawGroups; // used instead of the map when collapse=false
	private ArrayList<ReadGroup> groupList;
	private LinkedHashMap<String, Object[]> pending; // read name -> {mate, strand}
	private ReadGroup dummyGroup;

	private long totalRecords;
	private long nOut;

	private Integer regionStart;
	private Integer regionEnd;
	private String chrom;
	private Integer lastRefIdx;
	private double sampleProb;

	private int[] coverageByBlocks; // guarded by this; volatile double-check

	private Region revRegion;

	private char strand;

	private Stranded stranded;
	private boolean longReads;
	private boolean collapse;
	private long absCap;
	private int pendingCap;

	private int[] groupWeights;
	private char[] groupStrands;
	private volatile boolean finalized;

	public Region(double maxcov, double sample, char strand) {
		this(maxcov, sample, strand, Stranded.FR_UNSTRANDED, false, true, 4_000_000, 1_000_000);
	}

	public Region(double maxcov, double sample, char strand, Stranded stranded, boolean longReads, boolean collapse, long absCap, int pendingCap) {
		this.groups = collapse ? new LinkedHashMap<ReadGroup.Signature, ReadGroup>() : null;
		this.rawGroups = collapse ? null : new ArrayList<ReadGroup>();
		this.groupList = new ArrayList<ReadGroup>();
		this.pending = new LinkedHashMap<String, Object[]>();
		this.chrom = null;
		this.regionStart = null;
		this.regionEnd = null;
		this.sampleProb = 1.0;
		this.sample = sample;
		this.maxcov = maxcov;
		this.nOut = 0;
		this.totalRecords = 0;
		this.strand = strand;
		this.stranded = stranded;
		this.longReads = longReads;
		this.collapse = collapse;
		this.absCap = absCap;
		this.pendingCap = pendingCap;
		this.finalized = false;
	}

	public char getStrand() {
		return strand;
	}

	public void setRevRegion(Region region) {
		this.revRegion = region;
	}

	public Region getRevRegion() {
		return revRegion;
	}

	public String getChrom() {
		return chrom;
	}

	private void finalizeGroups() {
		if(finalized) {
			return;
		}
		synchronized(this) {
			if(finalized) {
				return;
			}
			// flush incomplete pairs as singleton fragments
			for(Map.Entry<String, Object[]> e : pending.entrySet()) {
				addFrag(new Mate[] {(Mate) e.getValue()[0]}, (Character) e.getValue()[1]);
			}
			pending.clear();
			if(dummyGroup != null) {
				groupList.add(dummyGroup);
			}
			if(groups != null) {
				groupList.addAll(groups.values());
				groups = null;
			}
			if(rawGroups != null) {
				groupList.addAll(rawGroups);
				rawGroups = null;
			}
			groupWeights = new int[groupList.size()];
			groupStrands = new char[groupList.size()];
			for(int i=0;i<groupList.size();i++) {
				ReadGroup g = groupList.get(i);
				g.gid = i;
				// dummy gap-fill reads all shared a single read index in the original
				// implementation, i.e. they count as one observation in quantification
				groupWeights[i] = g.dummy ? 1 : g.count;
				groupStrands[i] = g.strand;
			}
			finalized = true;
		}
	}

	public ArrayList<ReadGroup> getGroups() {
		finalizeGroups();
		return groupList;
	}

	/** fragment weight (number of collapsed fragments) per group id. */
	public int[] getGroupWeights() {
		finalizeGroups();
		return groupWeights;
	}

	/** strand per group id. */
	public char[] getGroupStrands() {
		finalizeGroups();
		return groupStrands;
	}

	public int getNumGroups() {
		finalizeGroups();
		return groupList.size();
	}

	public int getTheoreticalNumberOfReads() {
		synchronized(this) {
			return (int)Math.round( (double)totalRecords/sampleProb);
		}
	}

	public int[] getCoverageByBlocks() {
		int[] c = coverageByBlocks;
		if(c == null) {
			synchronized(this) {
				c = coverageByBlocks;
				if(c == null) {
					this.computeCoverageByBlocks();
					c = coverageByBlocks;
				}
			}
		}
		return c;
	}

	private ArrayList<ReadGroup> currentGroups() {
		// caller must hold this monitor (synchronized)
		ArrayList<ReadGroup> all = new ArrayList<ReadGroup>();
		if(finalized) {
			all.addAll(groupList);
			return all;
		}
		if(dummyGroup != null) {
			all.add(dummyGroup);
		}
		if(groups != null) {
			all.addAll(groups.values());
		}
		if(rawGroups != null) {
			all.addAll(rawGroups);
		}
		for(Object[] e : pending.values()) {
			ReadGroup pg = new ReadGroup();
			pg.strand = (Character) e[1];
			pg.add(new Mate[] {(Mate) e[0]}, 1);
			all.add(pg);
		}
		return all;
	}

	private int[] computeCoverageByReads() {
		synchronized(this) {
			int[] cov = new int[regionEnd-regionStart+1];
			for(ReadGroup g : currentGroups()) {
				for(int f=0;f<g.frags.size();f++) {
					for(Mate m : g.frags.get(f)) {
						for(int i=m.alignStart;i<m.alignEnd;i++) {
							cov[i-regionStart] += g.fragW[f];
						}
					}
				}
			}
			return cov;
		}
	}

	private void computeCoverageByBlocks() {
		// caller holds this monitor (via getCoverageByBlocks)
		int[] cov = new int[regionEnd-regionStart+1];
		for(ReadGroup g : currentGroups()) {
			for(int f=0;f<g.frags.size();f++) {
				for(Mate m : g.frags.get(f)) {
					for(int b=0;b<m.bRef.length;b++) {
						for(int i=m.bRef[b];i<m.bRef[b]+m.bLen[b];i++) {
							cov[i-regionStart] += g.fragW[f];
						}
					}
				}
			}
		}
		this.coverageByBlocks = cov;
	}


	public Region[] splitByCov(int maxLen) {
		finalizeGroups();
		int[] cov = computeCoverageByReads();

		IntList splitPoints = new IntList();
		splitPoints(splitPoints, 0, cov.length, cov, maxLen);

		splitPoints.add(0);
		splitPoints.add(cov.length);

		splitPoints.sort();

		int[] sar = splitPoints.toArray();

		Region[] regions = new Region[splitPoints.length()-1];
		for(int i=0;i<regions.length;i++) {
			regions[i] = new Region(maxcov, sample, strand, stranded, longReads, collapse, absCap, pendingCap);
			regions[i].chrom = this.chrom;
			regions[i].lastRefIdx = this.lastRefIdx;
		}

		for(ReadGroup g : groupList) {
			for(int f=0;f<g.frags.size();f++) {
				Mate[] mates = g.frags.get(f);
				// find the sub-region of each mate; a fragment is kept whole iff all
				// its mates fall into the same sub-region (as in the original code,
				// records not fully contained in any sub-region are dropped)
				int[] js = new int[mates.length];
				for(int k=0;k<mates.length;k++) {
					js[k] = -1;
					int astart = mates[k].alignStart-regionStart;
					int aend = mates[k].alignEnd-regionStart;
					for(int j=1;j<sar.length;j++) {
						if(astart >= sar[j-1] && aend < sar[j]) {
							js[k] = j-1;
							break;
						}
					}
				}
				boolean same = js[0] >= 0;
				for(int k=1;k<mates.length && same;k++) {
					same = js[k] == js[0];
				}
				if(same) {
					regions[js[0]].addGroupDirect(mates, g.fragW[f], g.strand, g.dummy);
				}else {
					for(int k=0;k<mates.length;k++) {
						if(js[k] >= 0) {
							regions[js[k]].addGroupDirect(new Mate[] {mates[k]}, g.fragW[f], g.strand, g.dummy);
						}
					}
				}
			}
		}

		// drop sub-regions that received no reads (e.g. intervals fully inside a long
		// intron, where every read spans the split boundaries); returning them would
		// crash downstream with a NullPointerException (upstream issue #75)
		int nNonEmpty = 0;
		for(Region reg : regions) {
			if(reg.getRegionStart() != null) {
				nNonEmpty++;
			}
		}
		if(nNonEmpty < regions.length) {
			Region[] nonEmpty = new Region[nNonEmpty];
			int j = 0;
			for(Region reg : regions) {
				if(reg.getRegionStart() != null) {
					nonEmpty[j++] = reg;
				}
			}
			return nonEmpty;
		}

		return regions;
	}

	private void splitPoints(IntList splitPoints, int start, int end, int[] cov, int maxLen) {
		int len = end-start+1;

		int idx = ToolBox.getMinIndex(start+(int)(len*0.2), start+(int)(len*0.8), cov);
		splitPoints.add(idx);

		if(idx - start > maxLen) {
			splitPoints(splitPoints, start, idx, cov, maxLen);
		}
		if(end - idx > maxLen) {
			splitPoints(splitPoints, idx, end, cov, maxLen);
		}

	}


	public void addRead(SAMRecord read) {
		synchronized(this) {
			addReadInternal(read);
		}
	}

	private void addReadInternal(SAMRecord read) {
		if(sampleProb != 1.0 && r.nextDouble() >= sampleProb) {
			nOut++;
			return;
		}

		if(chrom == null) {
			chrom = read.getReferenceName();
		}else {
			if(!chrom.equals(read.getReferenceName())) {
				throw new RuntimeException();
			}
		}
		if(regionStart == null || read.getAlignmentStart()< regionStart) {
			regionStart = read.getAlignmentStart();
		}
		if(regionEnd == null || read.getAlignmentEnd()> regionEnd) {
			regionEnd = read.getAlignmentEnd();
		}
		lastRefIdx = read.getReferenceIndex();
		totalRecords++;

		Mate mate = ReadGroup.convert(read, longReads);
		char st = computeStrand(read);

		if("dummy".equals(read.getReadName())) {
			if(dummyGroup == null) {
				dummyGroup = new ReadGroup();
				dummyGroup.dummy = true;
				dummyGroup.strand = st;
			}
			dummyGroup.add(new Mate[] {mate}, 1);
		}else if(!read.getReadPairedFlag()) {
			addFrag(new Mate[] {mate}, st);
		}else {
			Object[] first = pending.remove(read.getReadName());
			if(first == null) {
				pending.put(read.getReadName(), new Object[] {mate, st});
				if(pending.size() > pendingCap) {
					// pathological depth: force-complete the eldest entry as singleton
					Iterator<Map.Entry<String, Object[]>> it = pending.entrySet().iterator();
					Map.Entry<String, Object[]> eldest = it.next();
					it.remove();
					addFrag(new Mate[] {(Mate) eldest.getValue()[0]}, (Character) eldest.getValue()[1]);
				}
			}else {
				// last record wins for the strand, as in the original id map overwrite
				addFrag(new Mate[] {(Mate) first[0], mate}, st);
			}
		}

		long limit = Math.min((long)(maxcov*(regionEnd - regionStart+1)), absCap);
		if(totalRecords > limit) {
			downsample();
		}

	}

	private void addFrag(Mate[] mates, char st) {
		ReadGroup g;
		if(collapse) {
			ReadGroup.Signature sig = new ReadGroup.Signature(mates, st);
			g = groups.get(sig);
			if(g == null) {
				g = new ReadGroup();
				g.strand = st;
				groups.put(sig, g);
			}
		}else {
			g = new ReadGroup();
			g.strand = st;
			rawGroups.add(g);
		}
		g.add(mates, 1);
	}

	private void addGroupDirect(Mate[] mates, int w, char st, boolean dummy) {
		if(dummy) {
			if(dummyGroup == null) {
				dummyGroup = new ReadGroup();
				dummyGroup.dummy = true;
				dummyGroup.strand = st;
			}
			dummyGroup.add(mates, w);
			totalRecords += (long) w*mates.length;
			updateBounds(mates);
			return;
		}
		if(groups == null) {
			// splitByCov/merge path: sub-regions always merge by signature
			groups = new LinkedHashMap<ReadGroup.Signature, ReadGroup>();
		}
		ReadGroup.Signature sig = new ReadGroup.Signature(mates, st);
		ReadGroup g = groups.get(sig);
		if(g == null) {
			g = new ReadGroup();
			g.strand = st;
			groups.put(sig, g);
		}
		g.add(mates, w);
		totalRecords += (long) w*mates.length;
		updateBounds(mates);
	}

	private void updateBounds(Mate[] mates) {
		for(Mate m : mates) {
			if(regionStart == null || m.alignStart < regionStart) {
				regionStart = m.alignStart;
			}
			if(regionEnd == null || m.alignEnd > regionEnd) {
				regionEnd = m.alignEnd;
			}
		}
	}

	private static int stochasticRound(double x) {
		int f = (int) Math.floor(x);
		return f + (r.nextDouble() < x-f ? 1 : 0);
	}

	private void downsample() {
		sampleProb *= sample;

		long removed = 0;
		if(groups != null) {
			Iterator<Map.Entry<ReadGroup.Signature, ReadGroup>> it = groups.entrySet().iterator();
			while(it.hasNext()) {
				ReadGroup g = it.next().getValue();
				removed += downsampleGroup(g);
				if(g.count == 0) {
					it.remove();
				}
			}
		}else {
			Iterator<ReadGroup> it = rawGroups.iterator();
			while(it.hasNext()) {
				ReadGroup g = it.next();
				removed += downsampleGroup(g);
				if(g.count == 0) {
					it.remove();
				}
			}
		}
		if(dummyGroup != null) {
			removed += downsampleGroup(dummyGroup);
			if(dummyGroup.count == 0) {
				dummyGroup = null;
			}
		}
		totalRecords -= removed;
		nOut -= removed;
	}

	private long downsampleGroup(ReadGroup g) {
		long removed = 0;
		int nf = 0;
		for(int i=0;i<g.frags.size();i++) {
			int w = g.fragW[i];
			int nw = stochasticRound(w*sample);
			removed += (long)(w-nw)*g.frags.get(i).length;
			if(nw > 0) {
				g.fragW[nf] = nw;
				g.frags.set(nf, g.frags.get(i));
				nf++;
			}
		}
		while(g.frags.size() > nf) {
			g.frags.remove(g.frags.size()-1);
		}
		g.count = 0;
		for(int i=0;i<nf;i++) {
			g.count += g.fragW[i];
		}
		return removed;
	}

	private char computeStrand(SAMRecord sr) {
		char strand = '.';
		if(stranded == Stranded.FR_FIRST_STRAND) {
			strand = '-';
			if(sr.getReadPairedFlag() && sr.getFirstOfPairFlag() && sr.getReadNegativeStrandFlag()) {
				strand = '+';
			}else if(sr.getReadPairedFlag() && sr.getSecondOfPairFlag() && sr.getMateNegativeStrandFlag()) {
				strand = '+';
			}else if(!sr.getReadPairedFlag()) {
				if(sr.getReadNegativeStrandFlag()) {
					strand = '+';
				}else {
					strand = '-';
				}
			}
		}else if(stranded == Stranded.FR_SECOND_STRAND) {
			strand = '+';
			if(sr.getReadPairedFlag() && sr.getFirstOfPairFlag() && sr.getReadNegativeStrandFlag()) {
				strand = '-';
			}else if(sr.getReadPairedFlag() && sr.getSecondOfPairFlag() && sr.getMateNegativeStrandFlag()) {
				strand = '-';
			}else if(!sr.getReadPairedFlag()) {
				if(sr.getReadNegativeStrandFlag()) {
					strand = '-';
				}else {
					strand = '+';
				}
			}
		}
		return strand;
	}

	/**
	 * Rescale factor for abundances: 1/sampling probability of this region, i.e. the
	 * factor that restores read counts reduced by down-sampling back to (approximately)
	 * the original scale. 1.0 when no down-sampling happened.
	 */
	public double getScale() {
		return 1.0/sampleProb;
	}


	public Integer getReferenceIndex() {
		return lastRefIdx;
	}

	public Integer getRegionStart() {
		return regionStart;
	}

	public Integer getRegionEnd() {
		return regionEnd;
	}

	public String toString() {
		return chrom+" "+regionStart+"-"+regionEnd+": "+totalRecords;
	}

	public ReadGraph buildGraph(int minIntronLength, int maxGapFilled, int maxMM, boolean longReads) {
		finalizeGroups();
		ReadGraph sg = new ReadGraph(chrom, regionStart, regionEnd, minIntronLength);
		for(ReadGroup g : groupList) {
			sg.addRead(g, maxMM, longReads);
		}
		sg.finalizeGaps(maxGapFilled);
		return sg;
	}

	public void join(Region revTemp) {
		synchronized(this) {
			synchronized(revTemp) {
				joinInternal(revTemp);
			}
		}
	}

	private void joinInternal(Region revTemp) {
		if(! this.chrom.equals(revTemp.chrom ) ) {
			throw new RuntimeException("chromosomes do not match");
		}

		if(revTemp.totalRecords == 0) {
			return;
		}

		if(this.finalized || revTemp.finalized) {
			throw new RuntimeException("too late");
		}

		if(this.sampleProb == revTemp.sampleProb) {
			this.regionStart = Math.min(this.regionStart, revTemp.regionStart);
			this.regionEnd = Math.max(this.regionEnd, revTemp.regionEnd);
			mergeFrom(revTemp, 1.0);
		}else if(this.sampleProb <= revTemp.sampleProb) {
			double rel = this.sampleProb / revTemp.sampleProb;
			mergeFrom(revTemp, rel);
			this.regionStart = Math.min(this.regionStart, revTemp.regionStart);
			this.regionEnd = Math.max(this.regionEnd, revTemp.regionEnd);
		}else {
			double rel = revTemp.sampleProb / this.sampleProb;
			// re-home this region's content with reduced probability
			ArrayList<ReadGroup> old = new ArrayList<ReadGroup>();
			if(this.groups != null) {
				old.addAll(this.groups.values());
			}
			if(this.rawGroups != null) {
				old.addAll(this.rawGroups);
			}
			ReadGroup oldDummy = this.dummyGroup;
			this.groups = new LinkedHashMap<ReadGroup.Signature, ReadGroup>();
			this.rawGroups = null;
			this.dummyGroup = null;
			this.totalRecords = 0;
			rehome(old, oldDummy, rel);
			mergeFrom(revTemp, 1.0);
			this.regionStart = Math.min(this.regionStart, revTemp.regionStart);
			this.regionEnd = Math.max(this.regionEnd, revTemp.regionEnd);
			this.sampleProb = revTemp.sampleProb;
		}
	}

	private void rehome(ArrayList<ReadGroup> old, ReadGroup oldDummy, double keepProb) {
		for(ReadGroup g : old) {
			rehomeGroup(g, keepProb);
		}
		if(oldDummy != null) {
			ReadGroup ng = new ReadGroup();
			ng.dummy = true;
			ng.strand = oldDummy.strand;
			this.dummyGroup = ng;
			rehomeGroupInto(oldDummy, keepProb, ng);
		}
	}

	private void rehomeGroup(ReadGroup g, double keepProb) {
		ReadGroup ng = new ReadGroup();
		ng.strand = g.strand;
		rehomeGroupInto(g, keepProb, ng);
		if(ng.count > 0) {
			// signature of merged content is equivalent to the original group's
			groups.put(new ReadGroup.Signature(ng.frags.get(0), ng.strand), ng);
		}
	}

	private void rehomeGroupInto(ReadGroup g, double keepProb, ReadGroup ng) {
		for(int i=0;i<g.frags.size();i++) {
			int keep = stochasticRound(g.fragW[i]*keepProb);
			if(keep > 0) {
				ng.add(g.frags.get(i), keep);
				totalRecords += (long) keep*g.frags.get(i).length;
			}
		}
	}

	private void mergeFrom(Region o, double keepProb) {
		// flush the other's pending pairs as singletons (their mates are not in this region)
		o.finalizeGroups();
		for(ReadGroup g : o.groupList) {
			if(g.dummy) {
				for(int i=0;i<g.frags.size();i++) {
					int keep = keepProb >= 1.0 ? g.fragW[i] : stochasticRound(g.fragW[i]*keepProb);
					if(keep > 0) {
						addGroupDirect(g.frags.get(i), keep, g.strand, true);
					}
				}
			}else {
				for(int i=0;i<g.frags.size();i++) {
					int keep = keepProb >= 1.0 ? g.fragW[i] : stochasticRound(g.fragW[i]*keepProb);
					if(keep > 0) {
						addGroupDirect(g.frags.get(i), keep, g.strand, false);
					}
				}
			}
		}
		this.nOut += o.nOut;
	}

}
