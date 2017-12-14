import java.io.*;
import java.util.*;
import java.text.SimpleDateFormat;

public class allClust {
	public static void main(String[] args) throws IOException{
		String listFile = args[0];
		int thr = 23;	
		try {
			BufferedReader reader = new BufferedReader(new FileReader(new File(listFile)));
			
			List<String> tempList = new ArrayList<String>();
			String temp;
			while ((temp = reader.readLine()) != null) {
				tempList.add(temp);
			}
			reader.close();
			
			for (int i = 0;i < tempList.size();i++) {
				String tempLine = tempList.get(i);
				int status = 0;		
				String[] cut = tempLine.split("\t");
				String bc = cut[0];
				int len = bc.length();
				String bc1 = bc;
                                bc = bc.substring(0,len-6);
		
				//calculate distance matrix
				calDis caler = new calDis();
				int[][] disMat = caler.calculate(bc);

				int[][] disMat2 = new int[disMat.length][];
				for (int a = 0;a<disMat2.length;a++) {
					disMat2[a] = disMat[a].clone();
				}

				//run RCM
				for (int a = 0;a < disMat2.length;a++) {
					for (int b = 0; b < disMat2.length;b++) {
						int cutInt = disMat2[a][b];
						if (cutInt <= thr) {
							disMat2[a][b] = 1;
						}
						else {
							disMat2[a][b] = 0;
						}
					}
				}

				int[] rcmRe = new int[disMat2.length];
				if (disMat2.length > 1) { 
					RcmOrder rcm = new RcmOrder();
					rcmRe = rcm.ordering(disMat2);
				
				}
				else if (disMat2.length == 1){
					rcmRe[0] = 1;
				}

				//run group
				HashMap<Integer,Integer> rcmRc = new HashMap<Integer,Integer>();
				
				int[][] afRcm = new int[rcmRe.length][rcmRe.length];
				for (int j = 0;j < rcmRe.length;j++) {
					for (int k = 0;k <rcmRe.length;k++) {
						afRcm[j][k] = disMat[rcmRe[j]-1][rcmRe[k]-1];
					}
					
					rcmRc.put(j,rcmRe[j]-1);
				}


				RcmGroup grouper = new RcmGroup();
				ArrayList<int[]> subClust = grouper.clust(afRcm);

				String fasFile = "./bc_groups/"+bc+".fasta";
				String grpFile = "./levGroup/"+bc+"_group.res";			
				//read in sequence,name and index info
				HashMap<Integer,String> faInd = new HashMap<Integer,String>();
				HashMap<String,String> faSeq = new HashMap<String,String>();
				int ind = 0;
				String curName = "";

				BufferedReader seqReader = new BufferedReader(new FileReader(fasFile));

				String seqTemp;
				while ((seqTemp = seqReader.readLine()) != null && seqTemp.length() >0 ) {
					seqTemp = seqTemp.trim();
					String head = seqTemp.substring(0,1);
					if (head.equals(">")) {
						curName = seqTemp.substring(1);
						faInd.put(ind,curName);
						ind++;
					}
					else {
						faSeq.put(curName,seqTemp);
					}
				}
				seqReader.close();
				//transform group results
				FileWriter writer = new FileWriter(grpFile);
				BufferedWriter bw = new BufferedWriter(writer); 
				int gpID = 0;
				for (int[] re:subClust) {
					gpID++;
					for (int j = 0;j < re.length;j++) {
						int p = re[j];
						int orId = rcmRc.get(p);
						String orName = faInd.get(orId);
						String orSeq = faSeq.get(orName);
						bw.write(orName+"\t"+gpID);
						bw.newLine();
					}
				}
				bw.flush();
				bw.close();
				writer.close();
				
					
			}
		}catch(Exception e){
			e.printStackTrace();
		}
        }


	public static void updateList(String listFile, List<String> list) {
		try{
			FileWriter writer = new FileWriter(new File(listFile));
			BufferedWriter bw = new BufferedWriter(writer);
			Iterator<String> it = list.iterator();
			while(it.hasNext()) {
				String tempLine = it.next();
				bw.write(tempLine);
				bw.newLine();
			}
			bw.flush();
			bw.close();
			writer.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}
}

