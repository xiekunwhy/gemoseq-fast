package projects.gemoseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import projects.gemoma.ReadStats;

/**
 * Standalone computation of the read statistics (intron length mean/sd, mean read
 * length) used by GeMoSeq's splice-graph pruning. Writes them to a small text file
 * that can be fed back into GeMoSeq ("Read statistics" parameter) so that many
 * restricted runs can share one whole-genome statistic instead of recomputing it
 * per chromosome.
 */
public class ReadStatsTool implements JstacsTool {

	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<Parameter>();
		try {
			pars.add(new FileParameter("Mapped reads","Mapped Reads in BAM format, coordinate sorted","bam",true));
			pars.add(new SimpleParameter(DataType.STRING,"Output file","Output text file for the statistics",true,"readstats.stats"));
			pars.add(new SimpleParameter(DataType.INT,"Shortest intron length","Length of the shortest intron considered",true,new NumberValidator<Integer>(0, Integer.MAX_VALUE),10));
			pars.add(new FileParameter("Reference list","Optional: text file with one reference name per line; compute statistics only from these references","txt,list,tsv,csv,bed",false));
		} catch(ParameterException ex) {
			ex.printStackTrace();
		}
		return new ToolParameterSet(this.getToolName(),pars);
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {

		String bamFile = ((FileParameter)parameters.getParameterForName("Mapped reads")).getFileContents().getFilename();
		String out = (String) parameters.getParameterForName("Output file").getValue();
		int minIntronLength = (int) parameters.getParameterForName("Shortest intron length").getValue();

		String[] restrictRefs = null;
		Object listVal = parameters.getParameterForName("Reference list").getValue();
		if(listVal != null) {
			String listFile = ((FileParameter)parameters.getParameterForName("Reference list")).getFileContents().getFilename();
			java.util.LinkedList<String> refs = new java.util.LinkedList<String>();
			BufferedReader rd = new BufferedReader(new FileReader(listFile));
			String line;
			while( (line = rd.readLine()) != null ) {
				line = line.trim();
				if(line.isEmpty() || line.startsWith("#")) {
					continue;
				}
				String[] parts = line.split("\\s+");
				if(!parts[0].isEmpty() && !refs.contains(parts[0])) {
					refs.add(parts[0]);
				}
			}
			rd.close();
			if(refs.isEmpty()) {
				throw new Exception("Reference list file "+listFile+" contains no reference names");
			}
			restrictRefs = refs.toArray(new String[0]);
		}

		protocol.append("Computing read statistics from "+bamFile+(restrictRefs == null ? " (whole BAM)" : " ("+restrictRefs.length+" reference(s))")+"\n");
		ReadStats stats = restrictRefs == null ? new ReadStats(minIntronLength, 1.0, bamFile) : new ReadStats(minIntronLength, 1.0, restrictRefs, bamFile);

		// write directly to the final location: <outdir>/<o> (outdir defaults to ".")
		File dir = TranscriptPrediction.outdirForTemp != null ? TranscriptPrediction.outdirForTemp : new File(".");
		File target = new File(dir, out);
		File parent = target.getParentFile();
		if(parent != null && !parent.isDirectory()) {
			parent.mkdirs();
		}
		stats.toFile(target.getPath(), bamFile);
		protocol.append("meanSplit="+stats.getMeanSplitLength()+" sdSplit="+stats.getSdSplitLength()+" meanReadLen="+stats.getMeanReadLength()+"\n");
		protocol.append("Statistics written to "+target.getPath()+"\n");

		TextResult tr = new TextResult("Read Statistics", "Read statistics for GeMoSeq", new FileParameter.FileRepresentation(target.getAbsolutePath()), false, "txt", getToolName(), null, true);
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
	}

	@Override
	public String getToolName() {
		return "Read Statistics";
	}

	@Override
	public String getToolVersion() {
		return "1.2.3-fast";
	}

	@Override
	public String getShortName() {
		return "readstats";
	}

	@Override
	public String getDescription() {
		return "compute read statistics (intron length distribution) for GeMoSeq pruning";
	}

	@Override
	public String getHelpText() {
		return "";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

	@Override
	public ToolResult[] getTestCases(String path) {
		return null;
	}

	@Override
	public void clear() {
	}

	@Override
	public String[] getReferences() {
		return null;
	}
}
