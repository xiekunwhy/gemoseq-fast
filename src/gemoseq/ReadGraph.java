package projects.gemoseq;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

import de.jstacs.utils.IntList;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import projects.gemoma.ReadStats;

public class ReadGraph {

	String chrom;
	int regionStart;
	private int minIntronLength;
	private int lastIdx;
	
	Node[] nodes;
	LinkedList<Edge> edges;
	
	private static final Node EMPTY_NODE = new Node();

	public ReadGraph(String chrom, int regionStart, int regionEnd, int minIntronLength) {
		this.chrom = chrom;
		this.regionStart = regionStart;
		this.nodes = new Node[regionEnd-regionStart+1];
		for(int i=0;i<nodes.length;i++) {
			nodes[i] = new Node();
		}
		this.edges = new LinkedList<Edge>();
		this.minIntronLength = minIntronLength;

	}

	/** split-off graphs are read-only; positions outside the component share one empty node. */
	private ReadGraph(String chrom, int regionStart, int length, int minIntronLength, boolean spanMode) {
		this.chrom = chrom;
		this.regionStart = regionStart;
		this.nodes = new Node[length];
		Arrays.fill(this.nodes, EMPTY_NODE);
		this.edges = new LinkedList<Edge>();
		this.minIntronLength = minIntronLength;
	}
	
	public void finalizeGaps(int maximumGapFilled) {
		int nzero = -1;
		for(int i=0;i<nodes.length;i++) {
			if(nodes[i].getNumberOfReads()>0 || nodes[i].getNumberOfIncomingEdges()>0 || nodes[i].getNumberOfOutgoingEdges()>0) {
				if(nzero > 0 && nzero < maximumGapFilled) {
					for(int j=i-nzero;j<i;j++) {
						nodes[j].nReads++;
						nodes[j].nDummy++;
						Edge e = nodes[j].addOutgoing(j, j+1, edges, nodes);
						e.nReads++;
					}
					Edge e = nodes[i-nzero-1].addOutgoing(i-nzero-1, i-nzero, edges, nodes);
					e.nReads++;
				}
				nzero = 0;
			}else if(nzero >= 0){
				nzero++;
			}
		}
	}
	
	public void condense() {
		int off = 0;
		while(nodes[off].nOut==0 && nodes[off].nIn==0) {
			off++;
		}
		int end = nodes.length-1;
		while(nodes[end].nIn == 0 && nodes[end].nOut==0) {
			end--;
		}
		if(off > nodes.length*0.5 || end < nodes.length*0.5) {
			Node[] temp = new Node[end-off+1];
			System.arraycopy(nodes, off, temp, 0, temp.length);
			nodes = temp;
			regionStart += off;
			for(Edge e : edges) {
				e.relStart -= off;
				e.relEnd -= off;
			}
		}
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer();
		Iterator<Edge> edgeIt = edges.iterator();
		while(edgeIt.hasNext()) {
			Edge e = edgeIt.next();
			if(e.getLength()>1) {
				sb.append(e.toString(regionStart)+"\n");
			}
		}
		
		for(int i=0;i<nodes.length;i++) {
			sb.append("node"+i+": "+nodes[i].nReads+" ("+nodes[i].nDummy+") {");
			for(Edge edge : nodes[i].getOutgoingEdges()) {sb.append(edge.toString(regionStart)+", ");}
			sb.append("}\n");
		}
		
		return sb.toString();
	}
	
	public void addRead(SAMRecord read, int maxMM, boolean longReads) {
		throw new UnsupportedOperationException("use addRead(ReadGroup,...) — SAMRecord-based ingestion was replaced by compact fragment groups");
	}

	public void addRead(ReadGroup g, int maxMM, boolean longReads) {
		for(int f=0;f<g.frags.size();f++) {
			int w = g.fragW[f];
			for(ReadGroup.Mate m : g.frags.get(f)) {
				addMate(m, w, g.dummy, maxMM, longReads);
			}
		}
	}

	private void addMate(ReadGroup.Mate m, int w, boolean isDummy, int maxMM, boolean longReads) {
		int nb = m.bRef.length;
		if(nb > 1 && m.nmm > maxMM) {
			return;
		}

		int lastRelEnd = -1;

		for(int k=0;k<nb;k++) {
			int relStart = m.bRef[k]-regionStart;

			if(lastRelEnd > -1) {
				boolean isSplit = true;
				if(relStart != lastRelEnd+1) {
					if(longReads) {
						isSplit = (m.bFlag[k] & ReadGroup.FLAG_SKIPPED_SEEN) != 0;
					}else {
						isSplit = (m.bFlag[k] & ReadGroup.FLAG_SPLIT_SHORT) != 0;
					}
				}
				if(relStart-lastRelEnd < minIntronLength || !isSplit) {
					for(int i=lastRelEnd+1;i<=relStart;i++) {
						Edge e = nodes[i-1].addOutgoing(i-1, i, edges, nodes);
						e.nReads += w;
					}
				}else {
					Edge e = nodes[lastRelEnd+m.delBefore[k]].addOutgoing(lastRelEnd+m.delBefore[k], relStart-m.delOff[k], edges, nodes);
					e.nReads += w;
				}
			}

			for(int i=0;i<m.bLen[k];i++) {
				nodes[relStart+i].nReads += w;
				if(isDummy) {
					nodes[relStart+i].nDummy += w;
				}
				if(i>0) {
					Edge e = nodes[relStart+i-1].addOutgoing(relStart+i-1, relStart+i, edges, nodes);
					e.nReads += w;
				}
			}
			lastRelEnd = relStart+m.bLen[k]-1;
		}

	}
	
	
	private void remove(Collection<Edge> toRemove) {
		// 'edges' is a LinkedList: Collection.removeAll on it is O(n*m); do a single pass instead
		HashSet<Edge> set = new HashSet<Edge>(toRemove);
		for(Edge e : toRemove) {
			nodes[e.relStart].removeOutgoing(e.relEnd);
			nodes[e.relEnd].removeIncoming(e.relStart);
		}
		edges.removeIf(set::contains);
	}
	
	public void pruneByAbsoluteNumberOfReads(double minNum, boolean onlyIntrons) {
		Iterator<Edge> it = edges.iterator();
		LinkedList<Edge> toRemove = new LinkedList<ReadGraph.Edge>();
		while(it.hasNext()) {
			Edge e = it.next();
			if(( !onlyIntrons || e.getLength() >= minIntronLength) && e.getNumberOfReads() < minNum) {
				toRemove.add(e);
			}
		}
		remove(toRemove);
	}
	
	public void pruneBySplitLength(ReadStats stats) {
		Iterator<Edge> it = edges.iterator();
		LinkedList<Edge> toRemove = new LinkedList<ReadGraph.Edge>();
		while(it.hasNext()) {
			Edge e = it.next();
			int start = e.relStart;
			int end = e.relEnd;
			int num = e.getNumberOfReads();
			//double min = 2.0*(end-start+1 - stats.getMeanSplitLength())/stats.getSdSplitLength();
			if(!stats.isOK(end-start+1, num)) {
				toRemove.add(e);
			}
		}
		remove(toRemove);		
	}
	
	public void pruneByRelativeNumberOfReads(double percent, boolean onlyIntrons) {
		Iterator<Edge> it = edges.iterator();
		LinkedList<Edge> toRemove = new LinkedList<ReadGraph.Edge>();
		while(it.hasNext()) {
			Edge e = it.next();
			int start = e.relStart;
			int end = e.relEnd;
			double maxNum = Math.max(nodes[start].getNumberOfReads(),nodes[end].getNumberOfReads());
			if(( !onlyIntrons || e.getLength() >= minIntronLength) && e.getNumberOfReads() < maxNum * percent) {
				toRemove.add(e);
			}
		}
		remove(toRemove);
	}
	
	
	public void pruneAlternativeIntronsLong(int maxDist) {
		Iterator<Edge> it = edges.iterator();
		LinkedList<Edge> toRemove = new LinkedList<ReadGraph.Edge>();
		
		while(it.hasNext()) {
			Edge e = it.next();
			if(e.getLength() >= minIntronLength) {
				int max = 0;
				int start = e.relStart;
				for(int i=Math.max(0, start-maxDist);i<Math.min(nodes.length, start+maxDist);i++) {
					for(Edge second : nodes[i].getOutgoingEdges()) {
						if(second.relEnd == e.relEnd) {
							if(second.nReads > max) {
								max = second.nReads;
							}
						}
					}
				}
				int end = e.relEnd;
				for(int i=Math.max(0, end-maxDist);i<Math.min(nodes.length, end+maxDist);i++) {
					for(Edge second : nodes[i].getIncomingEdges()) {
						if(second.relStart == e.relStart) {
							if(second.nReads > max) {
								max = second.nReads;
							}
						}
					}
				}
				if(e.nReads*5 < max) {
					toRemove.add(e);
				}
			}
		}
		remove(toRemove);
	}
	
	
	public int[] proposeSplitByCoverageOnly(int minTranscriptLength, int minExonSize) {
		double[] covBefore = new double[nodes.length];
		double[] covAfter = new double[nodes.length];
		double[] nBefore = new double[nodes.length];
		double[] nAfter = new double[nodes.length];
		for(int i=0;i<nodes.length;i++) {
			double n = 0;
			for(Edge e : nodes[i].getIncomingEdges()) {
				int s = e.relStart;
				n += e.getNumberOfReads();
				covBefore[i] += ( covBefore[s] + nodes[s].getNumberOfReads() ) * e.getNumberOfReads();
				
				nBefore[i] += (nBefore[s]+1) * e.getNumberOfReads();
			}
			if(n>0) {
				covBefore[i] /= n;
				nBefore[i] /= n;
			}
		}
		for(int i=nodes.length-1;i>=0;i--) {
			int n = 0;
			for(Edge e : nodes[i].getOutgoingEdges()) {
				int s = e.relEnd;
				n += e.getNumberOfReads();
				covAfter[i] += ( covAfter[s] + nodes[s].getNumberOfReads() ) * e.getNumberOfReads();
				
				nAfter[i] += (nAfter[s]+1) * e.getNumberOfReads();
				
			}
			if(n>0) {
				covAfter[i] /= n;
				nAfter[i] /= n;
			}
		}
		
		for(int i=0;i<nodes.length;i++) {
			if(nBefore[i]>0) {
				covBefore[i] /= nBefore[i];
			}
			if(nAfter[i]>0) {
				covAfter[i] /= nAfter[i];
			}
		}
		
		double[] sliwin = new double[nodes.length];
		int l=50;
		sliwin[0] = 0;
		for(int i=0;i<Math.min(nodes.length, l);i++) {
			sliwin[0] += nodes[i].getNumberOfReads();
		}
		for(int i=1;i<nodes.length;i++) {
			int ls = Math.max(0, i-l);
			int le = Math.min(nodes.length-1,i+l-1);
			sliwin[i] = sliwin[i-1];
			if(ls-1>=0) {
				sliwin[i] -= nodes[ls-1].getNumberOfReads();
			}
			if(le<=nodes.length-1) {
				sliwin[i] += nodes[le].getNumberOfReads();
			}
		}
		for(int i=0;i<nodes.length;i++) {
			int ls = Math.max(0, i-l);
			int le = Math.min(nodes.length-1,i+l-1);
			sliwin[i] /= (le-ls+1);
			sliwin[i] -= i*(nodes.length-i)/(double)(nodes.length*nodes.length);
//			System.out.println((i+regionStart)+" "+sliwin[i]);
		}
		
		
		int start = 0;
		int end = 0;
		
		IntList locs = new IntList();
		
		for(int i=0;i<nodes.length;i++) {
			if(nodes[i].getNumberOfOutgoingEdges() == 1 &&
					nodes[i].getOutgoingEdges().iterator().next().getRelEnd() == i+1 &&
					nodes[i+1].getNumberOfIncomingEdges() == 1) {
				end = i;
			}else {
				
				if(start > minTranscriptLength && nodes.length-end > minTranscriptLength 
						&& nodes[i].getNumberOfIncomingEdges()>0 && nodes[i].getNumberOfOutgoingEdges()>0) {
					//System.out.println((regionStart+start)+"-"+(regionStart+end)+": "+start+" "+end+" "+nodes.length);
					int[] regCov = new int[end-start+1];
					int cov = 0;
					for(int j=start;j<=end;j++) {
						cov += nodes[j].getNumberOfReads();
						regCov[j-start] = cov;
					}

					double locMinCoverage = Double.POSITIVE_INFINITY;

					int locMinCoverageIdx = -1;

					for(int j=start+minExonSize;j<end-minExonSize;j++) {
						double covBeforeLoc = covBefore[j];
						double covAfterLoc = covAfter[j];
						double covLocal = sliwin[j];
//						System.out.println((regionStart+j)+": "+covBefore+" "+covAfter+" "+covLocal);
						if(covBeforeLoc > 0 && covAfterLoc > 0) {
							if( ( Math.abs( Math.log(covBeforeLoc) - Math.log(covAfterLoc) )/Math.log(2) > 2.0 && covLocal*2.0 < Math.min(covBeforeLoc, covAfterLoc)) || 
									covLocal/Math.min(covBeforeLoc, covAfterLoc) < 0.2 ) {
								if( covLocal/Math.min(covBeforeLoc, covAfterLoc) < locMinCoverage ) {
									locMinCoverageIdx = j;
									locMinCoverage = covLocal/Math.min(covBeforeLoc, covAfterLoc);
								}
							}
						}
					}

					
					
					if(locMinCoverageIdx > -1) {
						locs.add(locMinCoverageIdx+regionStart);
					}
				}
				start = i;
				end = i;
			}
		}
		
		
		
		
		int[] temp = locs.toArray();
		Arrays.sort(temp);
		return temp;
	}
	
	
	private static double[] sliwin(int l, int[] cov) {
		double[] sliwin = new double[cov.length];

		sliwin[0] = 0;
		for(int i=0;i<Math.min(cov.length, l);i++) {
			sliwin[0] += cov[i];
		}
		for(int i=1;i<cov.length;i++) {
			int ls = Math.max(0, i-l);
			int le = Math.min(cov.length-1,i+l-1);
			sliwin[i] = sliwin[i-1];
			if(ls-1>=0) {
				sliwin[i] -= cov[ls-1];
			}
			if(le<=cov.length-1) {
				sliwin[i] += cov[le];
			}
		}
		for(int i=0;i<cov.length;i++) {
			int ls = Math.max(0, i-l);
			int le = Math.min(cov.length-1,i+l-1);
			sliwin[i] /= (le-ls+1);
			sliwin[i] -= i*(cov.length-i)/(double)(cov.length*cov.length);
		}
		return sliwin;
	}
	
	public int[] proposeSplitByReadCoverage(int minTranscriptLength, int minExonSize) {
		int[] cov = new int[nodes.length];
		
		for(int i=0;i<nodes.length;i++) {
			cov[i] += nodes[i].getNumberOfReads()-nodes[i].getNumberOfDummyReads();
		}
		
		for(Edge e : edges) {
			int start = e.relStart;
			int end = e.relEnd;
			int ec = e.nReads;
			for(int i=start+1;i<end;i++) {
				cov[i] += ec;
			}
		}
		
		
		double[] sliwin = sliwin(50,cov);
		int l2 = 1000;
		double[] sliwin2 = sliwin(l2,cov);
		
		
		int lastNode = nodes.length;
		while( lastNode > 0 && nodes[lastNode-1].getNumberOfOutgoingEdges()==0 ) {
			lastNode--;
		}
		
		int start = 0;
		int end = 0;
		
		IntList locs = new IntList();
		
		for(int i=0;i<nodes.length;i++) {
//			System.out.print((regionStart+i)+" ");
			if(nodes[i].getNumberOfOutgoingEdges() == 1 &&
					nodes[i].getOutgoingEdges().iterator().next().getRelEnd() == i+1 &&
					nodes[i+1].getNumberOfIncomingEdges() == 1) {
				end = i;
			}else {
				if(
						nodes[i].getNumberOfIncomingEdges()>0 && nodes[i].getNumberOfOutgoingEdges()>0
						|| i==lastNode) {//TODO changed here
				

					LinkedList<Integer> locMinCoverageIdx = new LinkedList<>();
					LinkedList<Double> locMinCoverage = new LinkedList<>();
					LinkedList<Integer> locMaxFcIdx = new LinkedList<>();
					LinkedList<Double> locMaxFc = new LinkedList<>();

					for(int j=Math.max(start+minExonSize,minTranscriptLength);j<Math.min(end-minExonSize,nodes.length-minTranscriptLength);j++) {
						double covBeforeLoc = sliwin2[ Math.max(0,j-l2) ];
						double covAfterLoc = sliwin2[Math.min(nodes.length-1, j+l2)];
						double covLocal = sliwin[j];

						if(covBeforeLoc > 0 && covAfterLoc > 0) {
							
							double fc = Math.abs( Math.log(covBeforeLoc+1.0) - Math.log(covAfterLoc+1.0) )/Math.log(2);
							if( fc > 2.0 ) {
								
								if(locMaxFcIdx.size() == 0 || j-minTranscriptLength*2 > locMaxFcIdx.getLast()) {
									locMaxFcIdx.add(j);
									locMaxFc.add(fc);
								}else if(fc > locMaxFc.getLast()) {
									locMaxFc.set(locMaxFc.size()-1, fc);
									locMaxFcIdx.set(locMaxFcIdx.size()-1, j);
								}
								
							}
							
							if( (covLocal+1.0)/Math.min(covBeforeLoc+1.0, covAfterLoc+1.0) < 0.2) {
								if(locMinCoverageIdx.size() == 0 || j-minTranscriptLength*2 > locMinCoverageIdx.getLast()) {
									locMinCoverageIdx.add(j);
									locMinCoverage.add(covLocal);
								}else if( covLocal < locMinCoverage.getLast() ) {
									locMinCoverage.set(locMinCoverage.size()-1, covLocal);
									locMinCoverageIdx.set(locMinCoverageIdx.size()-1, j);
								}
							}
						}
					}

					
					
					if(locMinCoverageIdx.size() > 0) {
						for(int k=0;k<locMinCoverageIdx.size();k++) {
							locs.addConditional(locMinCoverageIdx.get(k)+regionStart);
						}
						
					}
					if(locMaxFcIdx.size() > 0) {
						double t = Double.NEGATIVE_INFINITY;
						if(locMaxFcIdx.size() > 3) {
							Double[] temp = locMaxFc.toArray(new Double[0]);
							Arrays.sort(temp);
							t = temp[temp.length-3];
							
						}
						Iterator<Integer> it = locMaxFcIdx.iterator();
						Iterator<Double> it2 = locMaxFc.iterator();
						while(it.hasNext()) {
							int locIdx = it.next();
							double val = it2.next();
							if(val > t) {
								if(locIdx > -1) {
									while(locIdx > minTranscriptLength && sliwin[locIdx-1] < sliwin[locIdx]) {
										double covBeforeLoc = sliwin2[ Math.max(0,locIdx-1-l2) ];
										double covAfterLoc = sliwin2[Math.min(nodes.length-1, locIdx-1+l2)];
										double fc = Math.abs( Math.log(covBeforeLoc+1.0) - Math.log(covAfterLoc+1.0) )/Math.log(2);
										if(fc > 2.0) {
											locIdx--;
										}else {
											break;
										}
									}
									while(locIdx < nodes.length-minTranscriptLength && sliwin[locIdx+1] < sliwin[locIdx]) {
										double covBeforeLoc = sliwin2[ Math.max(0,locIdx+1-l2) ];
										double covAfterLoc = sliwin2[Math.min(nodes.length-1, locIdx+1+l2)];
										double fc = Math.abs( Math.log(covBeforeLoc+1.0) - Math.log(covAfterLoc+1.0) )/Math.log(2);
										if(fc > 2.0) {
											locIdx++;
										}else {
											break;
										}
									}

									locs.addConditional(locIdx+regionStart);
								}
							}
						}
					}
				}
				start = i;
				end = i;
			}
		}
		
		
		int[] temp = locs.toArray();
		Arrays.sort(temp);
		return temp;
		
		
		
		
		
	}
	
	
	
	
	
	public int[] proposeSplitsByLocalCoverage(int minExonSize, int minTranscriptLength, double logCoverageFC, double percentLocalCoverageDrop, double minCoverage) {
		int start = 0;
		int end = 0;
		
		IntList locs = new IntList();
		
		for(int i=0;i<nodes.length;i++) {
			if(nodes[i].getNumberOfOutgoingEdges() == 1 &&
					nodes[i].getOutgoingEdges().iterator().next().getRelEnd() == i+1 &&
					nodes[i+1].getNumberOfIncomingEdges() == 1) {
				end = i;
			}else {
				
				if(start > minTranscriptLength && nodes.length-end > minTranscriptLength && end-minExonSize > start+minExonSize) {
					int[] regCov = new int[end-start+1];
					int cov = 0;
					for(int j=start;j<=end;j++) {
						cov += nodes[j].getNumberOfReads();
						regCov[j-start] = cov;
					}

					double locMinCoverage = Double.POSITIVE_INFINITY;
					double minBefore = 0;
					double minAfter = 0;
					int locMinCoverageIdx = -1;

					for(int j=start+minExonSize;j<end-minExonSize;j++) {
						double covBefore = regCov[j-start]/(double)(j-start);
						double covAfter = (cov - regCov[j-start])/(double)(end-j-1);
						double covLocal = nodes[j].getNumberOfReads();
						if(covBefore > minCoverage && covAfter > minCoverage) {
							if( ( Math.abs( Math.log(covBefore) - Math.log(covAfter) )/Math.log(2) > logCoverageFC && covLocal*2.0 < Math.min(covBefore, covAfter)) || 
									covLocal/Math.min(covBefore, covAfter) < percentLocalCoverageDrop ) {
								if( covLocal/Math.min(covBefore, covAfter) < locMinCoverage ) {
									locMinCoverageIdx = j;
									locMinCoverage = covLocal/Math.min(covBefore, covAfter);
									minBefore = covBefore;
									minAfter = covAfter;
								}
							}
						}
					}

					
					
					if(locMinCoverageIdx > -1) {
						locs.add(locMinCoverageIdx+regionStart);
					}
				}
				start = i;
				end = i;
			}
		}
		return locs.toArray();
	}
	
	
	public void pruneByLocalCoverage(int minExonSize, int minTranscriptLength, double logCoverageFC, double percentLocalCoverageDrop, double minCoverage) {
		
		int start = 0;
		int end = 0;
		for(int i=0;i<nodes.length;i++) {
			if(nodes[i].getNumberOfOutgoingEdges() == 1 &&
					nodes[i].getOutgoingEdges().iterator().next().getRelEnd() == i+1 &&
					nodes[i+1].getNumberOfIncomingEdges() == 1) {
				end = i;
			}else {
				
				if(start > minTranscriptLength && nodes.length-end > minTranscriptLength && end-minExonSize > start+minExonSize) {
					//System.out.println((regionStart+start)+"-"+(regionStart+end)+": "+start+" "+end+" "+nodes.length);
					int[] regCov = new int[end-start+1];
					int cov = 0;
					for(int j=start;j<=end;j++) {
						cov += nodes[j].getNumberOfReads();
						regCov[j-start] = cov;
					}

					double locMinCoverage = Double.POSITIVE_INFINITY;
					double minBefore = 0;
					double minAfter = 0;
					int locMinCoverageIdx = -1;

					for(int j=start+minExonSize;j<end-minExonSize;j++) {
						double covBefore = regCov[j-start]/(double)(j-start);
						double covAfter = (cov - regCov[j-start])/(double)(end-j-1);
						double covLocal = nodes[j].getNumberOfReads();

						if(covBefore > minCoverage && covAfter > minCoverage) {
							if( ( Math.abs( Math.log(covBefore) - Math.log(covAfter) )/Math.log(2) > logCoverageFC && covLocal*2.0 < Math.min(covBefore, covAfter)) || 
									covLocal/Math.min(covBefore, covAfter) < percentLocalCoverageDrop ) {
								if( covLocal/Math.min(covBefore, covAfter) < locMinCoverage ) {
									locMinCoverageIdx = j;
									locMinCoverage = covLocal/Math.min(covBefore, covAfter);
									minBefore = covBefore;
									minAfter = covAfter;
								}
							}
						}
					}

					if(locMinCoverageIdx > -1) {
						Collection<Edge> es = nodes[locMinCoverageIdx].getOutgoingEdges();
						if(es.size() > 1) {
							throw new RuntimeException("should not happen");
						}
						remove(es);
					}
				}
				start = i;
				end = i;
			}
		}
		
	}
	
	
	public ReadGraph nextSplit() {
		int idx = -1;
		for(int i=lastIdx;i<nodes.length;i++) {
			if(nodes[i] != null && nodes[i].getNumberOfOutgoingEdges() > 0) {
				idx = i;
				break;
			}
		}
		if(idx == -1) {
			return null;
		}

		lastIdx = idx;

		LinkedList<Integer> stack = new LinkedList<Integer>();
		stack.add(idx);

		ArrayList<Integer> compIdx = new ArrayList<Integer>();
		ArrayList<Node> compNodes = new ArrayList<Node>();
		LinkedList<Edge> compEdges = new LinkedList<Edge>();

		int min = nodes.length-1;
		int max = 0;

		while(stack.size() > 0) {
			int curr = stack.removeFirst();
			if(nodes[curr] != null) {
				Node n = nodes[curr];
				if(curr < min) {
					min = curr;
				}
				if(curr > max) {
					max = curr;
				}
				Collection<Edge> out = n.getOutgoingEdges();
				Collection<Edge> in = n.getIncomingEdges();

				for(Edge edge : out) {
					if(nodes[edge.relEnd] != null) {
						stack.add(edge.relEnd);
					}
					compEdges.add(edge);
				}
				for(Edge edge : in) {
					if(nodes[edge.relStart] != null) {
						stack.add(edge.relStart);
					}
				}

				compIdx.add(curr);
				compNodes.add(n);
				nodes[curr] = null;
			}
		}

		ReadGraph g = new ReadGraph(chrom, regionStart+min, max-min+1, minIntronLength, true);
		for(int i=0;i<compIdx.size();i++) {
			g.nodes[compIdx.get(i)-min] = compNodes.get(i);
		}
		for(Edge e : compEdges) {
			e.relStart -= min;
			e.relEnd -= min;
		}
		g.edges = compEdges;

		return g;

	}
	
	
	public static class Node{
		
		private HashMap<Integer,Edge> outgoing;
		private HashMap<Integer,Edge> incoming;
		private int nOut;
		private int nIn;
		
		private int nReads;
		private int nDummy;
		
		public Node() {
			// maps are created lazily on first edge to keep per-base memory low
			this.nReads = 0;
			this.nDummy = 0;
			this.nOut = 0;
			this.nIn = 0;
		}
		
		public int getNumberOfReads() {
			return nReads;
		}
		
		public int getNumberOfDummyReads() {
			return nDummy;
		}

		public void removeOutgoing(int relEnd) {
			if(outgoing != null && outgoing.remove(relEnd) != null) {
				nOut--;
			}

		}
		
		public void removeIncoming(int relStart) {
			if(incoming != null && incoming.remove(relStart) != null) {
				nIn--;
			}
		}
		
		public Edge addOutgoing(int relStart, int relEnd, LinkedList<Edge> list, Node[] nodes) {

			if(outgoing == null) {
				outgoing = new HashMap<Integer, ReadGraph.Edge>();
			}
			Edge e = outgoing.get(relEnd);
			if(e != null) {
				return e;
			}
			Edge newEdge = new Edge(relStart,relEnd);
			this.outgoing.put(relEnd,newEdge);
			nOut++;
			if(nodes[relEnd].incoming == null) {
				nodes[relEnd].incoming = new HashMap<Integer, ReadGraph.Edge>();
			}
			nodes[relEnd].incoming.put(relStart, newEdge);
			nodes[relEnd].nIn++;
			list.add(newEdge);
			return newEdge;
		}
		
		public int getNumberOfOutgoingEdges() {
			return nOut;//outgoing.size();
		}
		
		public int getNumberOfIncomingEdges() {
			return nIn;//incoming.size();
		}

		public Collection<Edge> getOutgoingEdges() {
			return outgoing == null ? java.util.Collections.<Edge>emptyList() : outgoing.values();
		}

		public Collection<Edge> getIncomingEdges() {
			return incoming == null ? java.util.Collections.<Edge>emptyList() : incoming.values();
		}
		
		
	}
	
	public static class Edge{
		
		private int relStart;
		private int relEnd;
		
		private int nReads;
		
		public Edge(int relStart, int relEnd) {
			this.relStart = relStart;
			this.relEnd = relEnd;
			this.nReads = 0;
		}
		
		public Edge(int relStart, int relEnd, int nReads) {
			this(relStart,relEnd);
			this.nReads = nReads;
		}
		
		public Edge clone() throws CloneNotSupportedException {
			Edge clone = (Edge) super.clone();
			clone.relStart = relStart;
			clone.relEnd = relEnd;
			clone.nReads = nReads;
			
			return clone;
		}
		
		
		public int getNumberOfReads() {
			return nReads;
		}
		
		public String toString(int regionStart) {
			return (regionStart+relStart)+" <-> "+(regionStart+relEnd)+": "+nReads;
		}
		
		
		public int getRelStart() {
			return relStart;
		}
		
		public int getRelEnd() {
			return relEnd;
		}
		
		public int getLength() {
			return relEnd - relStart - 1;
		}		
		
		
	}

	public LinkedList<Edge> getIntronEdges() {
		return edges.stream().filter(e -> e.getRelEnd()-e.getRelStart()>=minIntronLength).collect(Collectors.toCollection(LinkedList<Edge>::new));
	}
	
}
