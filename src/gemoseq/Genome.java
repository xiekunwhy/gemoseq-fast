package projects.gemoseq;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.HashSet;

import projects.gemoma.Tools;

public class Genome {

	public static Genome genome;
	public static HashMap<String,Character> code;
	public static boolean[][][] isStart;
	public static boolean[][][] isStop;

	private HashMap<String,char[]> chromosomes;

	public static void init(String path) throws Exception {
		init(path, null);
	}

	/**
	 * Initializes the genome, optionally loading only the given references.
	 * Loading is done with a streaming parser when a filter is given, so that
	 * unneeded chromosomes never occupy memory.
	 */
	public static void init(String path, HashSet<String> onlyRefs) throws Exception {
		genome = onlyRefs == null ? new Genome(path) : new Genome(path, onlyRefs);
		code = Tools.getCode(Genome.class.getClassLoader().getResourceAsStream( "projects/gemoseq/genetic_code.txt"));
		code.put("NNN", 'X');



		int maxIdx = 256;

		isStart = new boolean[maxIdx][maxIdx][maxIdx];
		isStop = new boolean[maxIdx][maxIdx][maxIdx];




		for(String key : code.keySet()) {
			Character val = code.get(key);
			if(val == '*') {
				char[] keys = key.toCharArray();
				isStop[keys[0]][keys[1]][keys[2]] = true;
			}else if(val == 'M') {
				char[] keys = key.toCharArray();
				isStart[keys[0]][keys[1]][keys[2]] = true;
			}
		}

	}

	public String[] getChromosomeNames() {
		if(chromosomes.isEmpty() && fai != null) {
			return fai.keySet().toArray(new String[0]);
		}
		return this.chromosomes.keySet().toArray(new String[0]);
	}

	private Genome(String path) throws Exception {
		this.chromosomes = new HashMap<String,char[]>();
		initFai(path);
		if(fai != null) {
			// lazy loading: chromosomes are read from the fasta on demand via the .fai index
			return;
		}
		HashMap<String,String> temp = Tools.getFasta(path, 5, ".*");
		for(String key : temp.keySet()) {
			this.chromosomes.put(key, temp.get(key).toCharArray());
		}
	}

	/** streaming fasta parser keeping only the requested references. */
	private Genome(String path, HashSet<String> onlyRefs) throws Exception {
		this.chromosomes = new HashMap<String,char[]>();
		BufferedReader rd = new BufferedReader(new FileReader(path));
		String line;
		StringBuilder sb = null;
		boolean keep = false;
		while( (line = rd.readLine()) != null ) {
			if(line.startsWith(">")) {
				if(keep && sb != null) {
					String name = currentName;
					chromosomes.put(name, sb.toString().toCharArray());
				}
				currentName = null;
				keep = false;
				String h = line.substring(1).trim();
				int sp = h.indexOf(' ');
				if(sp > 0) {
					h = h.substring(0, sp);
				}
				if(onlyRefs.contains(h)) {
					keep = true;
					currentName = h;
					sb = new StringBuilder();
				}
			}else if(keep) {
				sb.append(line.trim().toUpperCase());
			}
		}
		if(keep && sb != null && currentName != null) {
			chromosomes.put(currentName, sb.toString().toCharArray());
		}
		rd.close();
		for(String ref : onlyRefs) {
			if(!chromosomes.containsKey(ref)) {
				throw new Exception("reference "+ref+" not found in genome file "+path);
			}
		}
	}

	private String currentName = null;

	public char[] getChromosome(String id) {
		char[] c = chromosomes.get(id);
		if(c == null) {
			synchronized(this) {
				c = chromosomes.get(id);
				if(c == null) {
					c = lazyLoad(id);
					if(c != null) {
						chromosomes.put(id, c);
					}
				}
			}
		}
		return c;
	}

	/** lazily loads one chromosome from the genome fasta using the .fai index (saves the full-genome memory footprint). */
	private char[] lazyLoad(String id) {
		if(fai == null) {
			return null;
		}
		long[] rec = fai.get(id);
		if(rec == null) {
			return null;
		}
		try {
			java.io.RandomAccessFile raf = new java.io.RandomAccessFile(fastaPath, "r");
			long offset = rec[0];
			long len = rec[1];
			long lineBases = rec[2];
			long lineWidth = rec[3];
			StringBuilder sb = new StringBuilder((int)Math.min(len, Integer.MAX_VALUE));
			byte[] buf = new byte[8192];
			long remaining = len;
			long pos = offset;
			while(remaining > 0) {
				raf.seek(pos);
				int n = raf.read(buf, 0, (int)Math.min(buf.length, remaining + (remaining/lineBases+1)*(lineWidth-lineBases) + lineWidth));
				if(n < 0) {
					break;
				}
				for(int i=0;i<n && remaining>0;i++) {
					byte b = buf[i];
					if(b != '\n' && b != '\r') {
						sb.append(Character.toUpperCase((char)b));
						remaining--;
					}
				}
				pos += n;
			}
			raf.close();
			return sb.toString().toCharArray();
		} catch(Exception e) {
			throw new RuntimeException("lazy loading of chromosome "+id+" from "+fastaPath+" failed: "+e.getMessage());
		}
	}

	private void initFai(String path) {
		this.fastaPath = path;
		java.io.File f = new java.io.File(path + ".fai");
		if(!f.isFile()) {
			return;
		}
		try {
			fai = new HashMap<String, long[]>();
			java.io.BufferedReader rd = new java.io.BufferedReader(new java.io.FileReader(f));
			String line;
			while( (line = rd.readLine()) != null ) {
				String[] p = line.split("\t");
				if(p.length >= 5) {
					fai.put(p[0], new long[] {Long.parseLong(p[2]), Long.parseLong(p[1]), Long.parseLong(p[3]), Long.parseLong(p[4])});
				}
			}
			rd.close();
		} catch(Exception e) {
			fai = null;
		}
	}

	private HashMap<String, long[]> fai;
	private String fastaPath;

}
