import java.io.*;
import java.util.*;
import java.text.SimpleDateFormat;

public class calDis {
	public int[][] calculate (String bc) throws IOException{
		String divLine = "\r\n";
		String divItem = "\t";
		
		//read in FASTA file into a 'list'
		File fastaFile = new File("./bc_groups/"+bc+".fasta");
		readFasta fastaReader = new readFasta();
		fastaReader.readFastaFast(fastaFile);
		List<faUnit> faList = fastaReader.getSeqs();
		//System.out.println("reading sequences finished");
		
		//for filtering to unique sequences
                System.out.println("filtering start");
                printCurrentTime();
		dataProcessor dataProcessor = new dataProcessor();
		List<faUnit> uniqueList = dataProcessor.findUniqueWithHash(faList);
		//System.out.println("filtering to unique sequences finished");
		//printCurrentTime();		
	
		Map<String,Integer> seqIndex = new HashMap<String,Integer>();	
		seqIndex = dataProcessor.seqIndex;
		//for calculating Hamming distance
		int uniqueNum = uniqueList.size();
		int[][] dist_matrix = new int[uniqueNum][uniqueNum]; //distance matrix
		int clusThresh = 23; //clustering distance threshold
		for (int i = 0; i < uniqueNum;i++) {
			for (int j = i+1; j < uniqueNum; j++) {
				String seqA = uniqueList.get(i).getSeq();
				String seqB = uniqueList.get(j).getSeq();
				int hamDist = dataProcessor.calHam(seqA, seqB);
				dist_matrix[i][j] = hamDist;
				dist_matrix[j][i] = hamDist;
			}
			dist_matrix[i][i] = 0;
		}
		
		//for calculating levenstein distance
		int[][] dist_matrix_lev = new int[dist_matrix.length][];
		int[][] dist_mix = new int[dist_matrix.length][];
		for (int i=0;i<dist_matrix.length;i++) {
			dist_matrix_lev[i] = dist_matrix[i].clone();
			dist_mix[i] = dist_matrix[i].clone();
		}
		
		for (int i = 0; i <uniqueNum;i++) {
			for (int j = i+1; j< uniqueNum;j++) {
				if (dist_matrix[i][j] >clusThresh) {
					String seqA = uniqueList.get(i).getSeq();
					String seqB = uniqueList.get(j).getSeq();
					int levDist = dataProcessor.calLev(seqA, seqB);
					dist_matrix_lev[i][j] = levDist;
					dist_matrix_lev[j][i] = levDist;
					dist_mix[i][j] = levDist;
					dist_mix[j][i] = levDist;
				}
				else {
					dist_matrix_lev[i][j] = -1;
					dist_matrix_lev[j][i] = -1;
				}
			}
		}
		//System.out.println("calculating leven distance finished");
		//printCurrentTime();



		//write out distance matrix of original sequences - write out from uniqueList & dist_mix
		Iterator<faUnit> it2 = faList.iterator();
		String[] orgSeqList = new String[faList.size()];
		int count = 0;
		while(it2.hasNext()) {
			String seqName = it2.next().getName();
			orgSeqList[count] = seqName;
			count++;
		}
		int[][] dist_org = new int[orgSeqList.length][orgSeqList.length];
		
                for (int i = 0;i < orgSeqList.length;i++) {
                        for (int j = i+1;j < orgSeqList.length;j++) {
                                String iName = orgSeqList[i];
                                String jName = orgSeqList[j];
                                int iIndex = seqIndex.get(iName)-1;
                                int jIndex = seqIndex.get(jName)-1;
                                dist_org[i][j] = dist_mix[iIndex][jIndex];
				dist_org[j][i] = dist_mix[iIndex][jIndex];
                        }
                }
		//System.out.println("generating orgin matrix finished");
		//printCurrentTime();		


		return(dist_org);

		//System.out.println("origin matrix finished");
		//printCurrentTime();
	}

	private static void writeRe (File outputFile, String content) {
		try {
			FileWriter fileWriter = new FileWriter(outputFile);
			fileWriter.write(content);
			fileWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	private static void printCurrentTime() {
		SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		System.out.println(df.format(new Date()));
	}
}

