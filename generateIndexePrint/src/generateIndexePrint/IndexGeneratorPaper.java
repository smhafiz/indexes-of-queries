package generateIndexePrint;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

public class IndexGeneratorPaper {

	public static void main(String[] args) throws IOException {
		//		try {
		//			System.setOut(new PrintStream(new File("test_sanity.txt")));
		//		} catch (Exception e) {
		//			e.printStackTrace();
		//		}
		Integer max_file_size = 4800*1024;//System parameter
		int max_P_limit = 100;//system parameter
		int max_number_of_topics_per_paper = 20;//system parameter
		int most_frequent_topics = 1000;//66 has support 50, 190 has support 20,  395 has support 10
		long total_time_elapsed = 0;
		int top_k_rank = 4;
		int number_of_criteria = 1;


		String fileName = "complete_full_with_cites.txt";
		BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
		long startTime, estimatedTime;

		startTime = System.nanoTime();
		String line;
		Paper paper;
		ArrayList<Paper> paper_dataset = new ArrayList<Paper>();
		ArrayList<String> all_topics = new ArrayList<String>();
		int counter = 0;
		int number_of_true_citation = 0;
		BufferedWriter bw_row = new BufferedWriter(new FileWriter(new File("paper_file_sizes.txt")));
		while(true)
		{
			line = br.readLine();
			if(line == null)
				break;
			String[] parts = line.split("\t\t");

			paper = Paper.getPaper(parts,max_file_size,max_number_of_topics_per_paper);
			if(paper!=null) {
				bw_row.write(paper.file_size+"\n");
				counter++;
				paper.id = counter;
				paper_dataset.add(paper);
				if(paper.citation>0) number_of_true_citation++;
				else {
					Random rand = new Random();
					paper.citation = rand.nextInt(100);
				} 
				all_topics.addAll(paper.keywords);
			}
		}
		br.close();
		bw_row.close();
		estimatedTime = System.nanoTime() - startTime;
		total_time_elapsed+=estimatedTime;
		int R_limit = paper_dataset.size();
		System.out.println("Number of paper files (<="+max_file_size/(1024*1024)+"MB): " + R_limit + ", creates " + (float)max_file_size*R_limit/(1024*1024*1024) + " GB.");
		System.out.println("Metadata to Java Dataset: "+(double)estimatedTime/1000000 + " milliseconds");
		System.out.println("Number of crawled citation: "+ number_of_true_citation);
		
		startTime = System.nanoTime();
		Set<String> unique_topics = new HashSet<String>(all_topics);
		ArrayList<WordCountPair> sorted_topics = new ArrayList<WordCountPair>();
		for(String keyword:unique_topics){
			WordCountPair wcp = new WordCountPair();
			wcp.word = keyword;
			wcp.count = Collections.frequency(all_topics, keyword);
			sorted_topics.add(wcp);
		}
		Collections.sort(sorted_topics, (o1, o2) -> o2.count.compareTo(o1.count));
		estimatedTime = System.nanoTime() - startTime;
		total_time_elapsed+=estimatedTime;
		System.out.println("Find most frequent topic words: "+(double)estimatedTime/1000000 + " milliseconds");

		float average_elapsed_time = 0;
		int number_of_trials = 100;
		long[] all_elapsed = new long[number_of_trials];
		long sum_elapsed = 0;
		for(int f=0;f<number_of_trials;f++){
			startTime = System.nanoTime();
			int[][][][] sparse_matrices = new int[top_k_rank][number_of_criteria][most_frequent_topics][R_limit];
			for(int i=0;i<most_frequent_topics;i++){
				WordCountPair wcp = sorted_topics.get(i);
//				System.out.println(wcp.word+"->"+wcp.count);
				String topic_word = wcp.word;
				ArrayList<Paper> vid_ds_per_topic = prunePaperDS(paper_dataset,topic_word);
//								System.out.println("size: "+ vid_ds_per_topic.size());
				int[][][] topKRowsByFourCriteraForSingleKeyword = getTopKRowsByFourCriteraForSingleKeyword(vid_ds_per_topic,top_k_rank,R_limit,topic_word,number_of_criteria);
				for(int j=0;j<top_k_rank;j++){
					for(int k=0;k<number_of_criteria;k++){
						sparse_matrices[j][k][i] = topKRowsByFourCriteraForSingleKeyword[k][j];
					}
				}
			}
						writeCCSFilesForBatchTopK(sparse_matrices);
//						writeCCSFilesForTopKInEachCriterion(sparse_matrices);
			//			writeCCSFIlesForBatchQueries(sparse_matrices);
			all_elapsed[f] = System.nanoTime() - startTime;
			sum_elapsed+=all_elapsed[f];
		}
		average_elapsed_time = sum_elapsed/number_of_trials;
		double sum_power = 0;
		for(int lll=0;lll<number_of_trials;lll++){
			sum_power +=(all_elapsed[lll]-average_elapsed_time)*(all_elapsed[lll]-average_elapsed_time);
		}
		double s_d = Math.sqrt(sum_power/number_of_trials);
		total_time_elapsed+=average_elapsed_time;
		System.out.println("Generate indexes of queries files: "+(double)average_elapsed_time/1000000 + " milliseconds, with s.d.: " + s_d/1000000 );
		System.out.println("Summation of all time elapsed: "+(double)total_time_elapsed/1000000000 + " seconds");
	}
	//1_date, 1_citation, ...
	private static void writeCCSFilesForTopKInEachCriterion(int[][][][] sparse_matrices) throws IOException {
		String[] endings = {"date","citation"};

		for(int i=0;i<sparse_matrices.length;i++) {
			String initial = String.valueOf(i+1);
			for(int j=0;j<sparse_matrices[0].length;j++){
				int p = sparse_matrices[0][0].length;
				int r = sparse_matrices[0][0][0].length;
				int[] rowIndex = new int[p];
				int[] columnPtr = new int[r+1];
				BufferedWriter bw_row = new BufferedWriter(new FileWriter(new File(initial+"_"+endings[j]+"_"+String.valueOf(p)+"_"+String.valueOf(r)+".row")));
				bw_row.write(p+" ");
				BufferedWriter bw_col = new BufferedWriter(new FileWriter(new File(initial+"_"+endings[j]+"_"+String.valueOf(p)+"_"+String.valueOf(r)+".col")));
				bw_col.write(r+" "+columnPtr[0]+" ");
				int rowIndex_counter = 0;
				for(int k=0;k<r;k++){
					columnPtr[k+1]=columnPtr[k];
					for(int l=0;l<p;l++){
						if(sparse_matrices[i][j][l][k]==1){
							rowIndex[rowIndex_counter]=l;
							columnPtr[k+1]++;
							bw_row.write(rowIndex[rowIndex_counter]+" ");
							rowIndex_counter++;
						}
					}
					bw_col.write(columnPtr[k+1]+" ");

				}
				bw_row.close();
				bw_col.close();
			}
		}
	}
	//date_1, date_2, ...
	private static void writeCCSFilesForBatchTopK(int[][][][] sparse_matrices) throws IOException {
		String[] initials = {"citation","date"};
		//int[top_k_rank][number_of_criteria][most_frequent_topics][R_limit];
		for(int i=0;i<sparse_matrices[0].length;i++) {
			for(int j=0;j<sparse_matrices.length;j++){
				int p = sparse_matrices[0][0].length;
				int r = sparse_matrices[0][0][0].length;
				int[] rowIndex = new int[p];
				int[] columnPtr = new int[r+1];
				BufferedWriter bw_row = new BufferedWriter(new FileWriter(new File(initials[i]+"_"+String.valueOf(j+1)+"_"+String.valueOf(p)+"_"+String.valueOf(r)+".row")));
				bw_row.write(p+" ");
				BufferedWriter bw_col = new BufferedWriter(new FileWriter(new File(initials[i]+"_"+String.valueOf(j+1)+"_"+String.valueOf(p)+"_"+String.valueOf(r)+".col")));
				bw_col.write(r+" "+columnPtr[0]+" ");
				int rowIndex_counter = 0;
				for(int k=0;k<r;k++){
					columnPtr[k+1]=columnPtr[k];
					for(int l=0;l<p;l++){
						if(sparse_matrices[j][i][l][k]==1){
							rowIndex[rowIndex_counter]=l;
							columnPtr[k+1]++;
							bw_row.write(rowIndex[rowIndex_counter]+" ");
							rowIndex_counter++;
						}
					}
					bw_col.write(columnPtr[k+1]+" ");

				}
				bw_row.close();
				bw_col.close();
			}
		}

	}
	private static int[][][] getTopKRowsByFourCriteraForSingleKeyword(
			ArrayList<Paper> paper_dataset, int top_k, int R_limit,
			String topic_word, int number_of_criteria) {
		int[][][] topKRowsByFourCriteraForSingleKeyword = new int[number_of_criteria][top_k][R_limit];//first:criterion;second:sorting rank;third:database_row
		Collections.sort(paper_dataset, 
				(o1, o2) -> o2.citation.compareTo(o1.citation));
		for(int i=0;i<top_k;i++){
			topKRowsByFourCriteraForSingleKeyword[0][i][paper_dataset.get(i).id-1]=1;
		}
		if(number_of_criteria==1) return topKRowsByFourCriteraForSingleKeyword;
		Collections.sort(paper_dataset, 
				(o1, o2) -> o2.publish_date.compareTo(o1.publish_date));
		for(int i=0;i<top_k;i++){
			topKRowsByFourCriteraForSingleKeyword[1][i][paper_dataset.get(i).id-1]=1;
		}

		if(number_of_criteria==2) return topKRowsByFourCriteraForSingleKeyword;
		return topKRowsByFourCriteraForSingleKeyword;
	}
	private static ArrayList<Paper> prunePaperDS(
			ArrayList<Paper> paper_dataset, String word) throws IOException {
		ArrayList<Paper> prunedVidDS = new ArrayList<Paper>();
		for(Paper single_vid:paper_dataset){
			if(single_vid.keywords.contains(word) ) {// || single_vid.title.toLowerCase().contains(word)){
				prunedVidDS.add(single_vid);
				
			}
		}
		return prunedVidDS;
	}
	public static void generateCCSArrays(ArrayList<Paper> paper_dataset,int P_limit, int R_limit, int[] rowIndex, int[] columnPtr){
		int counter = 0;
		int[] columnOccupancy = new int[R_limit];
		for(Paper single_vid:paper_dataset){
			columnOccupancy[single_vid.id-1]=1;
			counter++;
			if(counter==P_limit) break;
		}
		columnPtr[0] = 0;
		counter = 0;
		for(int i=0;i<R_limit;i++){
			if(columnOccupancy[i]==1){
				for(int j=0;j<P_limit;j++) {
					if(paper_dataset.get(j).id-1==i){
						rowIndex[counter]=j;
						counter++;
						break;
					}
				}
			}
			columnPtr[i+1]=columnPtr[i]+columnOccupancy[i];
		}
	}
//	private static void writeCCSfilesForSingleKeyword(ArrayList<Paper> paper_dataset, int P_limit, int R_limit, long total_time_elapsed, String topic_word) throws IOException{
//		long startTime, estimatedTime = 0;
//
//
//		startTime = System.nanoTime();
//		Collections.sort(paper_dataset, 
//				(o1, o2) -> o2.publish_date.compareTo(o1.publish_date));
//		estimatedTime = System.nanoTime() - startTime;
//		total_time_elapsed+=estimatedTime;
//		System.out.println("Sort \'" + topic_word + "\'  based dataset (on date): "+(double)estimatedTime/1000000 + " microseconds");
//
//		/****************** By date *****************/
//		startTime = System.nanoTime();
//		int[] rowIndex = new int[P_limit];
//		int[] columnPtr = new int[R_limit+1];
//		generateCCSArrays(paper_dataset, P_limit, R_limit, rowIndex, columnPtr);
//		estimatedTime = System.nanoTime() - startTime;
//		total_time_elapsed+=estimatedTime;
//		System.out.println("CCS arrays creation (on date): "+(double)estimatedTime/1000000 + " microseconds");
//
//
//		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(topic_word+"_date.row")));
//		bw.write(P_limit+" ");
//
//		for(int i=0;i<rowIndex.length;i++)
//		{
//			bw.write(rowIndex[i]+" ");
//		}
//		bw.close();
//
//
//		bw = new BufferedWriter(new FileWriter(new File(topic_word+"_date.col")));
//		bw.write(R_limit+" ");
//
//		for(int i=0;i<columnPtr.length;i++)
//		{
//			bw.write(columnPtr[i]+" ");
//		}
//		bw.close();
//
//		/****************** By views *****************/
//		startTime = System.nanoTime();
//		Collections.sort(paper_dataset, 
//				(o1, o2) -> o2.citation.compareTo(o1.citation));
//		estimatedTime = System.nanoTime() - startTime;
//		System.out.println("Sort \'" + topic_word + "\'  based dataset (on views): "+(double)estimatedTime/1000000 + " microseconds");
//		startTime = System.nanoTime();
//		generateCCSArrays(paper_dataset, P_limit, R_limit, rowIndex, columnPtr);
//		estimatedTime = System.nanoTime() - startTime;
//		total_time_elapsed+=estimatedTime;
//		System.out.println("CCS arrays creation (on views): "+(double)estimatedTime/1000000 + " microseconds");
//
//		bw = new BufferedWriter(new FileWriter(new File(topic_word+"_views.row")));
//		bw.write(P_limit+" ");
//
//		for(int i=0;i<rowIndex.length;i++)
//		{
//			bw.write(rowIndex[i]+" ");
//		}
//		bw.close();
//
//		bw = new BufferedWriter(new FileWriter(new File(topic_word+"_views.col")));
//		bw.write(R_limit+" ");
//
//		for(int i=0;i<columnPtr.length;i++)
//		{
//			bw.write(columnPtr[i]+" ");
//		}
//		bw.close();
//
//		/****************** By file size *****************/
//		//		startTime = System.nanoTime();
//		//		Collections.sort(paper_dataset, 
//		//				(o1, o2) -> o2.file_size.compareTo(o1.file_size));
//		//		estimatedTime = System.nanoTime() - startTime;
//		//		System.out.println("Sort \'" + topic_word + "\'  based dataseton file size: "+(double)estimatedTime/1000000 + " microseconds");
//		//
//		//		startTime = System.nanoTime();
//		//		generateCCSArrays(paper_dataset, P_limit, R_limit, rowIndex, columnPtr);
//		//		estimatedTime = System.nanoTime() - startTime;
//		//		total_time_elapsed+=estimatedTime;
//		//		System.out.println("CCS arrays creation (on file size): "+(double)estimatedTime/1000000 + " microseconds");
//		//		
//		//		bw = new BufferedWriter(new FileWriter(new File(topic_word+"_size.row")));
//		//		bw.write(P_limit+" ");
//		//		
//		//		for(int i=0;i<rowIndex.length;i++)
//		//		{
//		//			bw.write(rowIndex[i]+" ");
//		//		}
//		//		bw.close();
//		//		
//		//		bw = new BufferedWriter(new FileWriter(new File(topic_word+"_size.col")));
//		//		bw.write(R_limit+" ");
//		//		
//		//		for(int i=0;i<columnPtr.length;i++)
//		//		{
//		//			bw.write(columnPtr[i]+" ");
//		//		}
//		//		bw.close();
//
//
//		System.out.println("Summation of all time elapsed: "+(double)total_time_elapsed/1000000000 + " seconds");
//
//	}
	


}
