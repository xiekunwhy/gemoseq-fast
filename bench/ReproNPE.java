import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

/** one gene with a 30kb intron; reads tile the exons (and splice over the intron) at ~100x. */
public class ReproNPE {
	public static void main(String[] args) throws Exception {
		String genome = args[0];
		String bam = args[1];
		Random r = new Random(7);
		// genome 100kb, gene: exon1 [20000,21000) intron 30kb exon2 [51000,52000)
		StringBuilder sb = new StringBuilder();
		Random rg = new Random(42);
		for(int i=0;i<100000;i++) sb.append("ACGT".charAt(rg.nextInt(4)));
		String seq = sb.toString();
		// splice sites GT..AG
		seq = seq.substring(0,21000) + "GT" + seq.substring(21002,50998) + "AG" + seq.substring(51000);
		java.io.FileWriter fw = new java.io.FileWriter(genome);
		fw.write(">chr1\n");
		for(int i=0;i<seq.length();i+=60) fw.write(seq.substring(i, Math.min(i+60, seq.length()))+"\n");
		fw.close();

		SAMFileHeader header = new SAMFileHeader();
		header.addSequence(new SAMSequenceRecord("chr1", 100000));
		header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
		SAMFileWriter w = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, new File(bam));
		int id = 0;
		ArrayList<SAMRecord> all = new ArrayList<SAMRecord>();
		for(int k=0;k<30000;k++) {
			int start = 19500 + r.nextInt(32000);
			int fragLen = 400;
			int e1s = 20000, e1e = 21000, e2s = 51000, e2e = 52000;
			int r1s = start, r2s = start+fragLen-100;
			all.add(rec(header, "f"+(id)+"/1", seq, r1s, e1s, e1e, e2s, e2e, r2s, true));
			all.add(rec(header, "f"+(id)+"/2", seq, r2s, e1s, e1e, e2s, e2e, r1s, false));
			id++;
		}
		all.sort((a,b) -> Integer.compare(a.getAlignmentStart(), b.getAlignmentStart()));
		for(SAMRecord rec : all) {
			w.addAlignment(rec);
		}
		w.close();
		System.out.println("wrote "+id+" pairs");
	}

	static SAMRecord rec(SAMFileHeader header, String name, String seq, int rs, int e1s, int e1e, int e2s, int e2e, int mateStart, boolean first) {
		SAMRecord rec = new SAMRecord(header);
		rec.setReadName(name);
		rec.setReferenceName("chr1");
		List<CigarElement> cig = new ArrayList<CigarElement>();
		int pos = rs;
		if(rs < e1e) {
			int l = Math.min(rs+100, e1e)-rs;
			cig.add(new CigarElement(l, CigarOperator.M));
			pos = rs+l;
		}
		if(pos < e2s && rs+100 > e2s) {
			// spliced over the intron
			cig.add(new CigarElement(e2s-e1e, CigarOperator.N));
			cig.add(new CigarElement(rs+100-e2s, CigarOperator.M));
		}else if(pos < rs+100) {
			cig.add(new CigarElement(rs+100-pos, CigarOperator.M));
		}
		rec.setAlignmentStart(rs);
		rec.setCigarString(cig.toString().replaceAll("[\\[\\] ]","").replaceAll("(\\d+)([A-Z])", "$1$2,").replaceAll(",$",""));
		// build cigar string manually instead
		StringBuilder cs = new StringBuilder();
		for(CigarElement ce : cig) cs.append(ce.getLength()).append(ce.getOperator().name().charAt(0));
		rec.setCigarString(cs.toString());
		rec.setReadBases(seq.substring(rs-1, Math.min(rs+99, seq.length())).getBytes());
		byte[] qual = new byte[100];
		java.util.Arrays.fill(qual, (byte)'I');
		rec.setBaseQualities(qual);
		rec.setMappingQuality(60);
		rec.setReadPairedFlag(true);
		rec.setProperPairFlag(true);
		rec.setMateReferenceName("chr1");
		rec.setMateAlignmentStart(mateStart);
		if(first) rec.setFirstOfPairFlag(true); else rec.setSecondOfPairFlag(true);
		return rec;
	}
}
