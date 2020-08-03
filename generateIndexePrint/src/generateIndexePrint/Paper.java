package generateIndexePrint;

import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;

public class Paper {
	int id;
	int serial_no;
	String title;
	ArrayList<String> authors;
	Date publish_date;
	ArrayList<String> keywords;
	Integer file_size;
	Integer citation;

	static DateFormat df = new SimpleDateFormat("d MMM yyyy");

	public static Paper getPaper(String[] parts, Integer max_file_size, int max_number_of_topics) {
		Paper paper = new Paper();
		paper.serial_no = Integer.parseInt(parts[0].trim());
		paper.title = parts[1].replaceAll("\"", "").trim();
		String[] author_name = parts[2].trim().split(" and ");
		paper.authors = new ArrayList<String>();
		if(author_name.length > 0)
		{	
			for(int i=0;i<author_name.length;i++){
				if(author_name[i].contains(",")) {
					String[] sublist = author_name[i].split(",");
					for(int j=0;j<sublist.length;j++){
						if(sublist[j].trim().length()>1) paper.authors.add(sublist[j].trim().toLowerCase());
					}

				} else if(author_name[i].trim().length()>1) paper.authors.add(author_name[i].trim().toLowerCase()); 
			}

		} else return null;
		try {
//			System.out.println(parts[3].trim());
			paper.publish_date = df.parse(parts[3].trim());
//			System.out.println(paper.serial_no);
//			System.out.println(paper.publish_date);
		} catch (ParseException e) {
			System.out.println(paper.serial_no);
			e.printStackTrace();
		}
		String[] keywords_list = parts[4].replaceAll("\"", "").trim().split("[/,]+");
		paper.keywords = new ArrayList<String>();

		if(keywords_list.length>0) {
			for(int i=0;i<keywords_list.length;i++){
				if(keywords_list[i].trim().length()>1) paper.keywords.add(keywords_list[i].trim().toLowerCase());
			}

		} else paper.keywords.add("Not mentioned");
		try {
			paper.file_size = Integer.parseInt(parts[5].trim());
			if(paper.file_size>max_file_size) {
				//System.out.println("s_n: " + paper.serial_no + ", size: " + paper.file_size);
				return null;
				}
			paper.citation = Integer.parseInt(parts[6].trim());
		} catch (Exception e) {
			System.out.println(paper.serial_no);
			e.printStackTrace();
		}
		
		/*if(paper.citation==0) { 
			Random rand = new Random(); 
			paper.citation = rand.nextInt(20)+1;
		}*/
		return paper;

	}

}
