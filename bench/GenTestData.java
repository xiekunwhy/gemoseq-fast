import htsjdk.samtools.*;
import java.io.*;
import java.util.*;

public class GenTestData {
    static final String CHR = "chr1";
    static final int CHR_LEN = 2000000;
    static final int READ_LEN = 100;
    static final int HOT_START = 1000000, HOT_END = 1100000;
    static final Random rng = new Random(42);

    static byte[] genome = new byte[CHR_LEN];
    static final ArrayList<Gene> genes = new ArrayList<Gene>();
    static final ArrayList<int[]> reserved = new ArrayList<int[]>(); // splice-site 2bp windows [s,e)
    static long fragCounter = 0;
    static final int[] cnt = new int[CHR_LEN]; // reused counting-sort buckets
    static final byte[] QUALS = new byte[READ_LEN];
    static { Arrays.fill(QUALS, (byte) 40); }

    static class Gene {
        String id;
        boolean plus, hotspot;
        ArrayList<int[]> exons = new ArrayList<int[]>(); // genomic order, 0-based [s,e)
        int start() { return exons.get(0)[0]; }
        int end() { return exons.get(exons.size() - 1)[1]; }
        int splicedLen() { int s = 0; for (int[] e : exons) s += e[1] - e[0]; return s; }
    }

    public static void main(String[] args) throws Exception {
        byte[] b = { 'A', 'C', 'G', 'T' };
        for (int i = 0; i < CHR_LEN; i++) genome[i] = b[rng.nextInt(4)];

        // 300 normal genes, starts ~6kb apart from 100,000, strand alternating
        boolean plus = true;
        int placed = 0;
        for (int i = 0; i < 300; i++) {
            int start = 100000 + i * 6000 + (rng.nextInt(2001) - 1000);
            Gene g = makeGene(start, plus, CHR_LEN - 1000);
            if (g == null) { System.err.println("WARN: normal gene " + i + " not placed"); continue; }
            g.id = String.format("g%04d", ++placed); g.hotspot = false; genes.add(g);
            plus = !plus;
        }
        System.err.println("normal genes placed: " + placed);

        // 30 dense hotspot genes in [1,000,000, 1,100,000], overlap allowed
        plus = true;
        int hplaced = 0;
        for (int i = 0; i < 30; i++) {
            Gene g = null;
            for (int t = 0; t < 500 && g == null; t++) {
                int start = HOT_START + rng.nextInt((HOT_END - HOT_START) - 18000);
                g = makeGene(start, plus, HOT_END);
            }
            if (g == null) { System.err.println("WARN: hotspot gene " + i + " not placed"); continue; }
            g.id = String.format("g%04d", ++placed); g.hotspot = true; genes.add(g); hplaced++;
            plus = !plus;
        }
        System.err.println("hotspot genes placed: " + hplaced);

        // splice-site dinucleotides: plus GT..AG, minus CT..AC (genomic left..right)
        for (Gene g : genes) {
            for (int k = 0; k + 1 < g.exons.size(); k++) {
                int is = g.exons.get(k)[1], ie = g.exons.get(k + 1)[0];
                if (g.plus) { genome[is] = 'G'; genome[is + 1] = 'T'; genome[ie - 2] = 'A'; genome[ie - 1] = 'G'; }
                else        { genome[is] = 'C'; genome[is + 1] = 'T'; genome[ie - 2] = 'A'; genome[ie - 1] = 'C'; }
            }
        }

        writeFasta(new File("genome.fa"));
        writeFai(new File("genome.fa.fai"));
        writeGtf(new File("truth.gtf"));

        SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord(CHR, CHR_LEN));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        buildBam(header, 30.0, 3000.0, new File("low.bam"));
        buildBam(header, 2000.0, 100000.0, new File("hi.bam"));
        System.err.println("ALL DONE");
    }

    static Gene makeGene(int start, boolean plus, int maxEnd) {
        for (int t = 0; t < 200; t++) {
            int nEx = 1 + rng.nextInt(8);
            ArrayList<int[]> ex = new ArrayList<int[]>();
            int p = start;
            for (int k = 0; k < nEx; k++) {
                int el = 80 + rng.nextInt(321); // exon 80..400
                ex.add(new int[] { p, p + el });
                p += el;
                if (k + 1 < nEx) p += 50 + rng.nextInt(1951); // intron 50..2000
            }
            if (ex.get(nEx - 1)[1] > maxEnd) continue;
            boolean bad = false;
            for (int k = 0; k + 1 < nEx && !bad; k++) {
                int is = ex.get(k)[1], ie = ex.get(k + 1)[0];
                if (conflict(is) || conflict(ie - 2)) bad = true;
            }
            if (bad) continue;
            for (int k = 0; k + 1 < nEx; k++) {
                int is = ex.get(k)[1], ie = ex.get(k + 1)[0];
                reserved.add(new int[] { is, is + 2 });
                reserved.add(new int[] { ie - 2, ie });
            }
            Gene g = new Gene(); g.plus = plus; g.exons = ex;
            return g;
        }
        return null;
    }

    static boolean conflict(int s) {
        for (int[] r : reserved) if (s < r[1] && s + 2 > r[0]) return true;
        return false;
    }

    static void writeFasta(File f) throws Exception {
        BufferedWriter w = new BufferedWriter(new FileWriter(f));
        w.write(">" + CHR); w.newLine();
        for (int i = 0; i < CHR_LEN; i += 80) {
            w.write(new String(genome, i, Math.min(80, CHR_LEN - i), "US-ASCII"));
            w.newLine();
        }
        w.close();
    }

    static void writeFai(File f) throws Exception {
        // header ">chr1\n" = 6 bytes; 80 bases + '\n' per line
        PrintWriter pw = new PrintWriter(new FileWriter(f));
        pw.println(CHR + "\t" + CHR_LEN + "\t6\t80\t81");
        pw.close();
    }

    static void writeGtf(File f) throws Exception {
        ArrayList<Gene> sorted = new ArrayList<Gene>(genes);
        Collections.sort(sorted, new Comparator<Gene>() { public int compare(Gene a, Gene b) { return a.start() - b.start(); } });
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)));
        for (Gene g : sorted) {
            String strand = g.plus ? "+" : "-";
            String tid = g.id + ".t1";
            pw.println(CHR + "\tGenTestData\ttranscript\t" + (g.start() + 1) + "\t" + g.end() + "\t.\t" + strand + "\t.\tgene_id \"" + g.id + "\"; transcript_id \"" + tid + "\";");
            for (int[] e : g.exons)
                pw.println(CHR + "\tGenTestData\texon\t" + (e[0] + 1) + "\t" + e[1] + "\t.\t" + strand + "\t.\tgene_id \"" + g.id + "\"; transcript_id \"" + tid + "\";");
        }
        pw.close();
    }

    // transcript offset -> genomic position (0-based). minus strand: offset 0 = highest genomic coord.
    static int[] tx2g(Gene g) {
        int[] m = new int[g.splicedLen()];
        int o = 0;
        if (g.plus) {
            for (int[] e : g.exons) for (int p = e[0]; p < e[1]; p++) m[o++] = p;
        } else {
            for (int i = g.exons.size() - 1; i >= 0; i--) {
                int[] e = g.exons.get(i);
                for (int p = e[1] - 1; p >= e[0]; p--) m[o++] = p;
            }
        }
        return m;
    }

    // genomic blocks (sorted, 0-based [s,e)) covered by transcript offsets [ro, ro+len)
    static int[][] blocksArr(int[] m, int ro, int len) {
        ArrayList<int[]> r = new ArrayList<int[]>(4);
        int s = m[ro], prev = m[ro];
        for (int i = 1; i < len; i++) {
            int p = m[ro + i];
            if (Math.abs(p - prev) == 1) { prev = p; continue; }
            r.add(new int[] { Math.min(s, prev), Math.max(s, prev) + 1 });
            s = p; prev = p;
        }
        r.add(new int[] { Math.min(s, prev), Math.max(s, prev) + 1 });
        Collections.sort(r, new Comparator<int[]>() { public int compare(int[] a, int[] b) { return a[0] - b[0]; } });
        return r.toArray(new int[0][]);
    }

    static SAMRecord makeRecord(Gene g, int[] m, int o, int f, boolean left, String name, SAMFileHeader header) {
        int roA = o, roB = o + f - READ_LEN;
        int sA = Math.min(m[roA], m[roA + READ_LEN - 1]);
        int sB = Math.min(m[roB], m[roB + READ_LEN - 1]);
        int roLeft = sA <= sB ? roA : roB;
        int roRight = roLeft == roA ? roB : roA;
        int ro = left ? roLeft : roRight;
        int roMate = left ? roRight : roLeft;
        int[][] bl = blocksArr(m, ro, READ_LEN);
        int start0 = bl[0][0];
        int endEx = bl[bl.length - 1][1];
        int mateStart0 = Math.min(m[roMate], m[roMate + READ_LEN - 1]);
        int mateEndEx = Math.max(m[roMate], m[roMate + READ_LEN - 1]) + 1;
        int span = Math.max(endEx, mateEndEx) - Math.min(start0, mateStart0);

        ArrayList<CigarElement> ce = new ArrayList<CigarElement>(4);
        byte[] seq = new byte[READ_LEN];
        int off = 0, prevEnd = -1;
        for (int i = 0; i < bl.length; i++) {
            if (i > 0) ce.add(new CigarElement(bl[i][0] - prevEnd, CigarOperator.N));
            int l = bl[i][1] - bl[i][0];
            ce.add(new CigarElement(l, CigarOperator.M));
            System.arraycopy(genome, bl[i][0], seq, off, l);
            off += l; prevEnd = bl[i][1];
        }

        SAMRecord r = new SAMRecord(header);
        r.setReadName(name);
        r.setReferenceName(CHR);
        r.setAlignmentStart(start0 + 1);
        r.setCigar(new Cigar(ce));
        r.setReadBases(seq);
        r.setBaseQualities(QUALS);
        r.setMappingQuality(60);
        r.setReadPairedFlag(true);
        r.setProperPairFlag(true);
        r.setReadUnmappedFlag(false);
        r.setMateUnmappedFlag(false);
        r.setMateNegativeStrandFlag(false);
        if (left) { r.setFirstOfPairFlag(true);  r.setInferredInsertSize(span); }
        else      { r.setSecondOfPairFlag(true); r.setInferredInsertSize(-span); }
        r.setMateReferenceName(CHR);
        r.setMateAlignmentStart(mateStart0 + 1);
        r.setAttribute("NM", 0);
        return r;
    }

    // writes one coordinate-sorted BAM for a single gene; false if gene too short for reads
    static boolean writeGeneBam(Gene g, double depth, SAMFileHeader header, File out) {
        int L = g.splicedLen();
        if (L < READ_LEN) return false;
        int[] m = tx2g(g);
        int nEx = g.exons.size();
        double[] mult = new double[nEx];
        for (int i = 0; i < nEx; i++) mult[i] = 0.2 + 2.8 * rng.nextDouble();

        // exons in transcript order, with per-exon abundance multipliers (non-uniform coverage)
        int[] txOff = new int[nEx], txLen = new int[nEx];
        double[] w = new double[nEx];
        double tot = 0; int o = 0;
        for (int i = 0; i < nEx; i++) {
            int gi = g.plus ? i : nEx - 1 - i;
            txOff[i] = o; txLen[i] = g.exons.get(gi)[1] - g.exons.get(gi)[0]; o += txLen[i];
            w[i] = txLen[i] * mult[gi]; tot += w[i];
        }

        int nfrag = Math.max(1, (int) Math.round(depth * L / 200.0));
        int[] fo = new int[nfrag], ff = new int[nfrag];
        int[] rstart = new int[2 * nfrag]; // [2i]=left read start, [2i+1]=right read start
        String[] names = new String[nfrag];
        for (int i = 0; i < nfrag; i++) {
            double x = rng.nextDouble() * tot;
            int k = 0;
            while (k < nEx - 1) { x -= w[k]; if (x <= 0) break; k++; }
            int oo = txOff[k] + rng.nextInt(txLen[k]);
            int fl = 250 + (int) Math.round(50 * rng.nextGaussian());
            if (fl < READ_LEN) fl = READ_LEN;
            if (fl > 500) fl = 500;
            if (fl > L) fl = L;
            if (oo + fl > L) oo = L - fl;
            fo[i] = oo; ff[i] = fl;
            names[i] = String.format("r%08d", ++fragCounter);
            int sA = Math.min(m[oo], m[oo + READ_LEN - 1]);
            int sB = Math.min(m[oo + fl - READ_LEN], m[oo + fl - 1]);
            if (sA <= sB) { rstart[2 * i] = sA; rstart[2 * i + 1] = sB; }
            else          { rstart[2 * i] = sB; rstart[2 * i + 1] = sA; }
        }

        // counting sort of all reads of this gene by alignment start
        Arrays.fill(cnt, 0);
        for (int j = 0; j < 2 * nfrag; j++) cnt[rstart[j]]++;
        int[] off = new int[CHR_LEN];
        int sum = 0;
        for (int p = 0; p < CHR_LEN; p++) { off[p] = sum; sum += cnt[p]; }
        int[] ord = new int[2 * nfrag];
        for (int j = 0; j < 2 * nfrag; j++) ord[off[rstart[j]]++] = j;
        off = null;

        SAMFileWriter wr = new SAMFileWriterFactory().makeBAMWriter(header, true, out);
        for (int j = 0; j < 2 * nfrag; j++) {
            int idx = ord[j];
            wr.addAlignment(makeRecord(g, m, fo[idx >> 1], ff[idx >> 1], (idx & 1) == 0, names[idx >> 1], header));
        }
        wr.close();
        return true;
    }

    static void mergeBams(ArrayList<File> ins, SAMFileHeader header, File out) throws Exception {
        int n = ins.size();
        SamReader[] rd = new SamReader[n];
        SAMRecordIterator[] it = new SAMRecordIterator[n];
        final SAMRecord[] cur = new SAMRecord[n];
        PriorityQueue<Integer> pq = new PriorityQueue<Integer>(new Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                int d = cur[a].getAlignmentStart() - cur[b].getAlignmentStart();
                return d != 0 ? d : a - b;
            }
        });
        for (int i = 0; i < n; i++) {
            rd[i] = SamReaderFactory.makeDefault().open(ins.get(i));
            it[i] = rd[i].iterator();
            if (it[i].hasNext()) { cur[i] = it[i].next(); pq.add(i); }
        }
        SAMFileWriter w = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, out);
        while (!pq.isEmpty()) {
            int i = pq.poll();
            w.addAlignment(cur[i]);
            if (it[i].hasNext()) { cur[i] = it[i].next(); pq.add(i); }
            else cur[i] = null;
        }
        w.close();
        for (int i = 0; i < n; i++) { it[i].close(); rd[i].close(); }
    }

    static void buildBam(SAMFileHeader header, double dNorm, double dHot, File out) throws Exception {
        File tmpDir = new File("tmp_genes");
        tmpDir.mkdirs();
        ArrayList<File> temps = new ArrayList<File>();
        ArrayList<Gene> sorted = new ArrayList<Gene>(genes);
        Collections.sort(sorted, new Comparator<Gene>() { public int compare(Gene a, Gene b) { return a.start() - b.start(); } });
        int skipped = 0;
        long t0 = System.currentTimeMillis();
        for (Gene g : sorted) {
            File tmp = new File(tmpDir, g.id + ".bam");
            if (writeGeneBam(g, g.hotspot ? dHot : dNorm, header, tmp)) temps.add(tmp);
            else { skipped++; tmp.delete(); }
            if (g.hotspot) System.err.println(out.getName() + ": " + g.id + " (" + (System.currentTimeMillis() - t0) / 1000 + "s)");
        }
        System.err.println(out.getName() + ": per-gene BAMs done, skipped=" + skipped + ", merging " + temps.size() + " files (" + (System.currentTimeMillis() - t0) / 1000 + "s)");
        long t1 = System.currentTimeMillis();
        mergeBams(temps, header, out);
        System.err.println(out.getName() + ": merged (" + (System.currentTimeMillis() - t1) / 1000 + "s)");
        for (File f : temps) f.delete();
        tmpDir.delete();
        validate(out);
    }

    static void validate(File bam) throws Exception {
        SamReader r = SamReaderFactory.makeDefault().open(bam);
        SAMRecordIterator it = r.iterator();
        long n = 0; int prev = -1, bad = 0, maxEnd = 0;
        while (it.hasNext()) {
            SAMRecord rec = it.next();
            n++;
            if (rec.getAlignmentStart() < prev) bad++;
            prev = rec.getAlignmentStart();
            if (rec.getReadLength() != READ_LEN) bad++;
            int cs = 0;
            for (CigarElement e : rec.getCigar().getCigarElements())
                if (e.getOperator() == CigarOperator.M) cs += e.getLength();
            if (cs != READ_LEN) bad++;
            maxEnd = Math.max(maxEnd, rec.getAlignmentEnd());
        }
        it.close(); r.close();
        System.err.println(bam + ": records=" + n + " bad=" + bad + " maxEnd=" + maxEnd + " sizeMB=" + bam.length() / 1048576);
    }
}
